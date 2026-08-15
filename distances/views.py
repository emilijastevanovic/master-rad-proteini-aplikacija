import hashlib
import logging
import textwrap

from django.conf import settings
from django.http import HttpResponse, JsonResponse
from django.core.cache import cache
from neo4j.exceptions import Neo4jError

from .neo4j_client import run_query
from .query_builder import build_graph3d_protein
from .stats_queries import stat_aa_composition, stat_aa_seq, stat_aa_heatmap_df, stat_ss_pairs, stat_segment_lengths, stat_segment_lengths_by_ss, stat_distance_summary, stat_outlier_ss, stat_ss_distribution, stat_range_summary
from .protein_analysis import make_heat_maps, make_segment_hist_png, make_segment_hist_by_ss_png, make_aa_type_heatmap_png

logger = logging.getLogger(__name__)

LATEX_HEATMAP_PROTEINS = {"1K3I", "1A17", "1G2Y", "1FEW"}
LATEX_HEATMAP_DIR = settings.BASE_DIR / "latex_heatmaps"
LATEX_SEGMENT_HIST_DIR = settings.BASE_DIR / "latex_segment_histograms"

STAT_HANDLERS = {
    "aa_composition": stat_aa_composition,
    "ss_distribution": stat_ss_distribution,
    "sequence": stat_aa_seq,
    "ss_pairs": stat_ss_pairs,
    "distance_summary": stat_distance_summary,
    "outlier_ss": stat_outlier_ss,
    "range_summary": stat_range_summary,
}

# Statistike koje pored osnovnih parametara primaju i filtere grafa
# (sekundarna struktura, k, izabrane aminokiseline).
GRAPH_FILTER_STATS = {"range_summary"}

STATS_CACHE_TTL = 60 * 60 * 6

# Parametri koje svaki handler STVARNO koristi u svom Cypher upitu — oni i
# čine ključ keša. Namerno se razlikuju po statistici: npr. sastav
# aminokiselina zavisi samo od proteina i lanca, pa promena tipa rastojanja
# ne sme da pravi novi unos u kešu.
#
# PAŽNJA: osim `range_summary`, nijedan handler ne koristi
# max_distance/min_distance (vidi stats_queries.py), zato ih kod ostalih nema
# ni u ključu. Ako se neki handler izmeni da poštuje te granice, MORA se dodati
# i ovde — inače će keš vraćati rezultat izračunat za drugi opseg. Statistika
# koja nije navedena u ovoj mapi se jednostavno ne kešira.
STAT_CACHE_PARAMS = {
    "aa_composition":   ("protein", "chain"),
    "ss_distribution":  ("protein", "chain"),
    "sequence":         ("protein", "chain"),
    "ss_pairs":         ("protein", "chain", "type"),
    "distance_summary": ("protein", "chain", "type"),
    "outlier_ss":       ("protein", "chain", "type"),
    # "type_raw", a ne "type": za range_summary prazan tip znači „svi tipovi",
    # što je drugačiji rezultat od izričito izabranog „caca"
    "range_summary":    ("protein", "chain", "type_raw", "max_distance", "min_distance",
                         "ss", "k", "aminonames"),
}


def _norm_stat_type(value):
    """Ista normalizacija kao u stats_queries — prazan tip znači 'caca'."""
    return value if (value and value != "prazno") else "caca"


# Uđe u ključ keša da bi izmena formule odbacila ranije zapamćene rezultate.
# Povećaj kad god se promeni način računanja neke statistike.
# v2: stat_outlier_ss više ne odseca partnere sa nižim indeksom.
# v3: stat_outlier_ss vraća i closest_residues.
STATS_FORMULA_VERSION = "v3"


def _stats_cache_key(name, key_params):
    raw = "|".join([STATS_FORMULA_VERSION, name] + [str(key_params[p]) for p in STAT_CACHE_PARAMS[name]])
    return "stats:" + hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _should_cache_stat(result):
    """Prazan rezultat se ne kešira — da uvoz novog dump-a ne bi ostao
    zaklonjen zapamćenim 'nema podataka' do isteka TTL-a."""
    if not isinstance(result, dict) or result.get("empty"):
        return False
    if result.get("total") == 0:
        return False

    lists = [v for v in result.values() if isinstance(v, list)]
    if lists and all(not v for v in lists):
        return False

    return True


def _clean_distance(value, label):
    """Validira granicu rastojanja i vraća normalizovan string (ili None).

    Prihvata i zarez kao decimalni znak, jer je na našoj tastaturi prirodno
    kucnuti „8,5". Vraća string, a ne float, da bi `0` i dalje bilo tretirano
    kao zadata granica — kao i pre ove provere.
    """
    if value is None:
        return None

    text = str(value).strip().replace(",", ".")
    if not text or text == "prazno":
        return None

    try:
        float(text)
    except ValueError:
        raise ValueError(f"{label}: očekivan je broj, uneto „{value}”")

    return text


def _clean_k(value):
    """Validira sekvencijalno rastojanje (ceo broj)."""
    if value is None:
        return None

    text = str(value).strip()
    if not text or text == "prazno":
        return None

    try:
        int(text)
    except ValueError:
        raise ValueError(f"sekvencijalno rastojanje: očekivan je ceo broj, uneto „{value}”")

    return text


def _safe_filename_part(value, default="all"):
    text = str(value or default).strip()
    if not text or text == "prazno":
        text = default

    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in text)


def save_latex_heatmap_png(png, protein, chain, type_, max_distance, min_distance, aminonames=None,
                           suffix="distance_heatmap.png"):
    protein = protein.strip().upper()
    if protein not in LATEX_HEATMAP_PROTEINS:
        return None

    LATEX_HEATMAP_DIR.mkdir(parents=True, exist_ok=True)

    # izbor aminokiselina ulazi u ime da mapa podskupa ne bi pregazila mapu
    # celog proteina; kod dužih izbora ime se skraćuje i dopunjuje hešom da
    # dva različita izbora ne bi završila u istoj datoteci
    if aminonames:
        aa_part = _safe_filename_part("-".join(aminonames))
        if len(aa_part) > 40:
            aa_part = aa_part[:40] + "-" + hashlib.sha256(aa_part.encode()).hexdigest()[:8]
        aa_part = "aa-" + aa_part
    else:
        aa_part = "aa-sve"

    filename = "_".join([
        protein,
        _safe_filename_part(chain),
        _safe_filename_part(type_, "caca"),
        f"max-{_safe_filename_part(max_distance)}",
        f"min-{_safe_filename_part(min_distance)}",
        aa_part,
        suffix,
    ])
    path = LATEX_HEATMAP_DIR / filename
    path.write_bytes(png)
    logger.info("saved LaTeX heatmap: %s", path)
    return path


def save_latex_segment_hist_png(png, protein, chain, by_ss=False):
    protein = protein.strip().upper()
    if protein not in LATEX_HEATMAP_PROTEINS:
        return None

    LATEX_SEGMENT_HIST_DIR.mkdir(parents=True, exist_ok=True)

    filename = "_".join([
        protein,
        _safe_filename_part(chain),
        "segment_lengths_by_ss_histogram.png" if by_ss else "segment_lengths_histogram.png",
    ])
    path = LATEX_SEGMENT_HIST_DIR / filename
    path.write_bytes(png)
    logger.info("saved LaTeX segment histogram: %s", path)
    return path

CHAINS_CACHE_TTL = 60 * 60 * 6


def chains(request):
    """Lanci dostupni za zadati protein, sa brojem aminokiselina po lancu.

    Koristi se da se padajuća lista lanaca popuni tek kad korisnik izabere
    protein — umesto ranije fiksne liste koja je nudila i lance kojih u tom
    proteinu nema.
    """
    protein = request.GET.get("protein")
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    protein = protein.strip().upper()
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    cache_key = "chains:" + hashlib.sha256(protein.encode("utf-8")).hexdigest()
    cached = cache.get(cache_key)
    if cached is not None:
        logger.debug("chains cache hit: %s", protein)
        return JsonResponse(cached)

    try:
        rows = run_query("""
            MATCH (a:AminoAcid)
            WHERE a.protein = $protein AND a.chain IS NOT NULL
            RETURN a.chain AS chain, count(a) AS n
            ORDER BY chain
        """, protein=protein)
    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)

    data = {
        "protein": protein,
        "chains": [{"chain": r["chain"], "n": r["n"]} for r in rows],
    }

    # prazan rezultat se ne kešira — nepostojeći protein ili tek uvezen dump ne
    # sme da ostane zapamćen kao „nema lanaca" do isteka TTL-a
    if data["chains"]:
        cache.set(cache_key, data, CHAINS_CACHE_TTL)

    return JsonResponse(data)


def graph3d_protein(request):
    protein = request.GET.get("protein")
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)
    protein = protein.strip().upper()
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    aminoname_raw = request.GET.get("aminoname", "")
    try:
        params = {
            "protein": protein,
            "aminonames": [a for a in aminoname_raw.split(",") if a],
            "ss": request.GET.get("ss", "/"),
            "k": _clean_k(request.GET.get("k")),
            "type": request.GET.get("type", "prazno"),
            "max_distance": _clean_distance(request.GET.get("max_distance"), "gornja granica"),
            "min_distance": _clean_distance(request.GET.get("min_distance"), "donja granica"),
            "chain": request.GET.get("chain"),
        }
    except ValueError as e:
        return JsonResponse({"error": "invalid_parameter", "detail": str(e)}, status=400)

    nodes = {}
    edges = []

    q, cy_params = build_graph3d_protein(params)

    logger.debug("cypher: %s | params: %s", q, cy_params)

    try:
        rows = run_query(q, **cy_params)
    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)
    logger.debug("row count: %d", len(rows))

    for r in rows:

        if "idx" in r:  
            nid = f"{r['protein']}:{r['chain']}:{r['idx']}"
            if nid not in nodes:
                nodes[nid] = {
                    "id": nid,
                    "protein": r["protein"],
                    "chain": r["chain"],
                    "index": r["idx"],
                    "aa": r["aa"],
                    "ss": r.get("ss"),
                    "x": r.get("x"),
                    "y": r.get("y"),
                    "z": r.get("z"),
                }

        else:
            n1 = f"{r['protein']}:{r['chain']}:{r['index']}"
            n2 = f"{r['protein2']}:{r['chain2']}:{r['index2']}"

            if n1 not in nodes:
                nodes[n1] = {
                    "id": n1,
                    "protein": r["protein"],
                    "chain": r["chain"],
                    "index": r["index"],
                    "aa": r["aa"],
                    "ss": r.get("ss1"),
                    "x": r.get("x"),
                    "y": r.get("y"),
                    "z": r.get("z"),
                }

            if n2 not in nodes:
                nodes[n2] = {
                    "id": n2,
                    "protein": r["protein2"],
                    "chain": r["chain2"],
                    "index": r["index2"],
                    "aa": r["aa2"],
                    "ss": r.get("ss2"),
                    "x": r.get("x2"),
                    "y": r.get("y2"),
                    "z": r.get("z2"),
                }

            edge_id = f"{n1}--{n2}"
            edges.append({
                "id": edge_id,
                "source": n1,
                "target": n2,
                "type": r.get("type"),
                "value": r.get("value"),
            })

    return JsonResponse({
        "nodes": list(nodes.values()),
        "edges": edges
    })


def stats(request):
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))
    type_ = request.GET.get("type")
    max_distance = request.GET.get("max_distance")
    min_distance = request.GET.get("min_distance")

    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    include = request.GET.get("include")
    if include:
        requested = [x.strip() for x in include.split(",") if x.strip()]
    else:
        requested = ["aa_composition"]

    # filteri grafa — koristi ih samo statistika izabranog opsega
    try:
        graph_filters = {
            "ss": request.GET.get("ss"),
            "k": _clean_k(request.GET.get("k")),
            "max_distance": _clean_distance(request.GET.get("max_distance"), "gornja granica"),
            "min_distance": _clean_distance(request.GET.get("min_distance"), "donja granica"),
            "aminonames": [a for a in request.GET.get("aminoname", "").split(",") if a],
        }
    except ValueError as e:
        return JsonResponse({"error": "invalid_parameter", "detail": str(e)}, status=400)

    out = {"protein": protein, "chain": chain or None}

    # normalizovano samo za ključ keša — handleri i dalje dobijaju sirove
    # vrednosti i sami ih normalizuju, pa se ponašanje ne menja
    key_params = {
        "protein": protein.strip().upper(),
        "chain": chain,
        "type": _norm_stat_type(type_),
        "type_raw": type_ or "",
        "ss": graph_filters["ss"] if graph_filters["ss"] not in (None, "", "prazno", "/") else "",
        "k": graph_filters["k"],
        "max_distance": graph_filters["max_distance"],
        "min_distance": graph_filters["min_distance"],
        "aminonames": ",".join(sorted(graph_filters["aminonames"])),
    }

    for key in requested:
        fn = STAT_HANDLERS.get(key)
        if not fn:
            out[key] = {"error": "unknown stat"}
            continue

        cache_key = _stats_cache_key(key, key_params) if key in STAT_CACHE_PARAMS else None

        if cache_key:
            cached = cache.get(cache_key)
            if cached is not None:
                logger.debug("stats cache hit: %s (%s)", key, cache_key)
                out[key] = cached
                continue

        kwargs = dict(
            protein=protein,
            chain=chain,
            type=type_,
            max_distance=max_distance,
            min_distance=min_distance,
        )
        if key in GRAPH_FILTER_STATS:
            kwargs.update(
                max_distance=graph_filters["max_distance"],
                min_distance=graph_filters["min_distance"],
                ss=graph_filters["ss"],
                k=graph_filters["k"],
                aminonames=graph_filters["aminonames"],
            )

        result = fn(**kwargs)

        if cache_key and _should_cache_stat(result):
            cache.set(cache_key, result, STATS_CACHE_TTL)

        out[key] = result

    return JsonResponse(out)

def _clean_chain(val):
    return val if (val and val != "prazno") else None


def _aa_subtitle(aminonames, prefix):
    """Podnaslov mape sa spiskom aminokiselina, prelomljen po redovima.

    Kod izbora od dvadesetak imena jedan red bi bio širi od same slike, pa bi
    matplotlib odsekao kraj spiska.
    """
    if not aminonames:
        return None
    return prefix + "\n".join(textwrap.wrap(", ".join(aminonames), width=64))


def heatmap_distance(request):
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))
    type_ = request.GET.get("type", "caca")
    if not type_ or type_ == "prazno":
        type_ = "caca"
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    # Kad je označen samo podskup aminokiselina, crta se zasebna mapa samo za
    # njihove međusobne parove. Sortirano i bez duplikata, da isti izbor uvek
    # daje isti ključ keša bez obzira na redosled klikanja.
    aminonames = sorted({a.strip().upper() for a in request.GET.get("aminoname", "").split(",") if a.strip()})

    protein = protein.strip().upper()
    key_raw = f"{protein}|{chain}|{type_}|{','.join(aminonames)}"
    cache_key = "heatmap:" + hashlib.sha256(key_raw.encode("utf-8")).hexdigest()

    cached_png = cache.get(cache_key)
    if cached_png:
        logger.debug("heatmap cache hit: %s", cache_key)
        save_latex_heatmap_png(cached_png, protein, chain, type_, None, None, aminonames)
        return HttpResponse(cached_png, content_type="image/png")

    logger.debug("heatmap cache miss: %s", cache_key)

    df = stat_aa_heatmap_df(
        protein=protein,
        chain=chain,
        type=type_,
        max_distance=None,
        min_distance=None,
        aminonames=aminonames or None,
    )

    if df.empty or "aa1" not in df.columns:
        return JsonResponse({"error": "no_data", "detail": "No distance data found for this protein/chain/type."}, status=404)

    subtitle = _aa_subtitle(aminonames, "selektovane aminokiseline: ")
    imgs = make_heat_maps(df, protein=protein, subtitle=subtitle)
    png = imgs["distance_png"]

    cache.set(cache_key, png, timeout=60 * 60 * 24)
    save_latex_heatmap_png(png, protein, chain, type_, None, None, aminonames)

    return HttpResponse(png, content_type="image/png")


def heatmap_aa_types(request):
    """Toplotna mapa prosečnih rastojanja po VRSTI aminokiseline.

    `heatmap_distance` daje po jedan red i kolonu za svaki reziduum, pa je za
    ceo protein matrica velika koliko i lanac. Ovde je jedinica vrsta: za izbor
    GLU, LYS, THR matrica je 3×3, a ćelija je prosek preko svih parova te dve
    vrste u proteinu.
    """
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))
    type_ = request.GET.get("type", "caca")
    if not type_ or type_ == "prazno":
        type_ = "caca"
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    # sortirano i bez duplikata — isti izbor daje isti ključ keša bez obzira na
    # redosled klikanja, isto kao kod heatmap_distance
    aminonames = sorted({a.strip().upper() for a in request.GET.get("aminoname", "").split(",") if a.strip()})

    protein = protein.strip().upper()
    key_raw = f"aatypes|{protein}|{chain}|{type_}|{','.join(aminonames)}"
    cache_key = "heatmap_aa:" + hashlib.sha256(key_raw.encode("utf-8")).hexdigest()

    cached_png = cache.get(cache_key)
    if cached_png:
        logger.debug("aa-type heatmap cache hit: %s", cache_key)
        save_latex_heatmap_png(cached_png, protein, chain, type_, None, None, aminonames,
                               suffix="aa_type_heatmap.png")
        return HttpResponse(cached_png, content_type="image/png")

    logger.debug("aa-type heatmap cache miss: %s", cache_key)

    df = stat_aa_heatmap_df(
        protein=protein,
        chain=chain,
        type=type_,
        max_distance=None,
        min_distance=None,
        aminonames=aminonames or None,
    )

    if df.empty or "aa1" not in df.columns:
        return JsonResponse({"error": "no_data", "detail": "No distance data found for this protein/chain/type."}, status=404)

    subtitle = _aa_subtitle(aminonames, "izabrane aminokiseline: ")
    png = make_aa_type_heatmap_png(df, protein=protein, subtitle=subtitle)

    if png is None:
        return JsonResponse({"error": "no_data", "detail": "No pairs for the selected amino acids."}, status=404)

    cache.set(cache_key, png, timeout=60 * 60 * 24)
    save_latex_heatmap_png(png, protein, chain, type_, None, None, aminonames,
                           suffix="aa_type_heatmap.png")

    return HttpResponse(png, content_type="image/png")


def segments_histogram(request):
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))

    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    protein = protein.strip().upper()

    key_raw = f"segments|{protein}|{chain}"
    cache_key = "seg:" + hashlib.sha256(key_raw.encode()).hexdigest()

    cached = cache.get(cache_key)
    if cached:
        logger.debug("segments cache hit: %s", cache_key)
        save_latex_segment_hist_png(cached, protein, chain)
        return HttpResponse(cached, content_type="image/png")

    lengths = stat_segment_lengths(protein, chain)
    if not lengths:
        return JsonResponse({"error": "no_data"}, status=404)

    png = make_segment_hist_png(lengths, protein=protein)
    cache.set(cache_key, png, timeout=60 * 60 * 24)
    save_latex_segment_hist_png(png, protein, chain)

    return HttpResponse(png, content_type="image/png")


def segments_histogram_by_ss(request):
    """Isti histogram kao /heatmap/segments/, ali razložen po SS tipovima."""
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))

    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    protein = protein.strip().upper()

    key_raw = f"segments_by_ss|{protein}|{chain}"
    cache_key = "segss:" + hashlib.sha256(key_raw.encode()).hexdigest()

    cached = cache.get(cache_key)
    if cached:
        logger.debug("segments-by-ss cache hit: %s", cache_key)
        save_latex_segment_hist_png(cached, protein, chain, by_ss=True)
        return HttpResponse(cached, content_type="image/png")

    groups = stat_segment_lengths_by_ss(protein, chain)
    if not groups:
        return JsonResponse({"error": "no_data"}, status=404)

    png = make_segment_hist_by_ss_png(groups, protein=protein)
    cache.set(cache_key, png, timeout=60 * 60 * 24)
    save_latex_segment_hist_png(png, protein, chain, by_ss=True)

    return HttpResponse(png, content_type="image/png")


def global_stats(request):
    cached = cache.get("global_stats")
    if cached:
        return JsonResponse(cached)

    try:
        aa_row = run_query("""
            MATCH (a:AminoAcid)
            WITH a.protein AS p, a.chain AS c, count(a) AS len
            RETURN count(DISTINCT p)                   AS proteins,
                   sum(len)                            AS amino_acids,
                   round(avg(len), 1)                  AS avg_len,
                   round(percentileCont(len, 0.5), 1)  AS median_len,
                   min(len)                            AS min_len,
                   max(len)                            AS max_len
        """)

        ss_rows = run_query("""
            MATCH (a:AminoAcid)
            WHERE a.ss IS NOT NULL
            RETURN a.ss AS ss, count(a) AS n
            ORDER BY n DESC
        """)

        top_aa_rows = run_query("""
            MATCH (a:AminoAcid)
            WITH a.name AS name, count(a) AS n
            ORDER BY n DESC
            LIMIT 5
            RETURN name, n
        """)

        dist_rows = run_query("""
            MATCH ()-[r:DISTANCE]->()
            RETURN r.type AS type,
                   count(r)                   AS n,
                   round(avg(r.value), 2)     AS mean,
                   round(stdev(r.value), 2)   AS std,
                   round(min(r.value), 2)     AS min_val,
                   round(max(r.value), 2)     AS max_val
            ORDER BY type
        """)

        aa        = aa_row[0] if aa_row else {}
        total_aa  = aa.get("amino_acids") or 1
        ss_total  = sum(r["n"] for r in ss_rows) or 1

        data = {
            "proteins":        aa.get("proteins", 0),
            "amino_acids":     aa.get("amino_acids", 0),
            "avg_length":      aa.get("avg_len", 0),
            "median_length":   aa.get("median_len", 0),
            "min_length":      aa.get("min_len", 0),
            "max_length":      aa.get("max_len", 0),
            "ss_distribution": [
                {"ss": r["ss"], "n": r["n"],
                 "pct": round(r["n"] * 100 / ss_total, 1)}
                for r in ss_rows
            ],
            "top_aa": [
                {"name": r["name"], "n": r["n"],
                 "pct": round(r["n"] * 100 / total_aa, 1)}
                for r in top_aa_rows
            ],
            "dist_stats": [
                {"type": r["type"], "n": r["n"], "mean": r["mean"],
                 "std": r["std"], "min": r["min_val"], "max": r["max_val"]}
                for r in dist_rows
            ],
        }

        cache.set("global_stats", data, timeout=60 * 60 * 6)
        return JsonResponse(data)

    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)
