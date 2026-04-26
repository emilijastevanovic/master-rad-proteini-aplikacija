import hashlib
import logging

from django.http import HttpResponse, JsonResponse
from django.core.cache import cache
from neo4j.exceptions import Neo4jError

from .neo4j_client import run_query
from .query_builder import build_graph3d_protein
from .stats_queries import stat_aa_composition, stat_aa_seq, stat_aa_heatmap_df, stat_ss_pairs, stat_segment_lengths, stat_distance_summary, stat_outlier_ss
from .protein_analysis import make_heat_maps, make_segment_hist_png

logger = logging.getLogger(__name__)

STAT_HANDLERS = {
    "aa_composition": stat_aa_composition,
    "sequence": stat_aa_seq,
    "ss_pairs": stat_ss_pairs,
    "distance_summary": stat_distance_summary,
    "outlier_ss": stat_outlier_ss,
}

def graph3d_protein(request):
    protein = request.GET.get("protein")
    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    aminoname_raw = request.GET.get("aminoname", "")
    params = {
        "protein": protein,
        "aminonames": [a for a in aminoname_raw.split(",") if a],
        "ss": request.GET.get("ss", "/"),
        "k": request.GET.get("k"),
        "type": request.GET.get("type", "prazno"),
        "max_distance": request.GET.get("max_distance"),
        "min_distance": request.GET.get("min_distance"),
        "chain": request.GET.get("chain"),
    }

    nodes = {}
    edges = []

    # generiše Cypher i parametre
    q, cy_params = build_graph3d_protein(params)

    logger.debug("cypher: %s | params: %s", q, cy_params)

    # izvršavanje Neo4j upita
    try:
        rows = run_query(q, **cy_params)
    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)
    logger.debug("row count: %d", len(rows))

    for r in rows:

        # --- NODES ---
        if "idx" in r:  
            # slučaj bez distance (samo a)
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
            # slučaj DISTANCE (a1, a2)
            n1 = f"{r['protein']}:{r['chain']}:{r['index']}"
            n2 = f"{r['protein2']}:{r['chain2']}:{r['index2']}"

            # NODE 1
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

            # NODE 2
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

            # --- EDGES ---
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
        # default: učitaj samo ovo odmah
        requested = ["aa_composition"]

    out = {"protein": protein, "chain": chain or None}

    for key in requested:
        fn = STAT_HANDLERS.get(key)
        if not fn:
            out[key] = {"error": "unknown stat"}
            continue

        out[key] = fn(
            protein=protein,
            chain=chain,
            type=type_,
            max_distance=max_distance,
            min_distance=min_distance
        )

    return JsonResponse(out)

def _clean_chain(val):
    return val if (val and val != "prazno") else None


def heatmap_distance(request):
    protein = request.GET.get("protein")
    chain = _clean_chain(request.GET.get("chain"))
    type_ = request.GET.get("type", "caca")
    if not type_ or type_ == "prazno":
        type_ = "caca"
    max_distance = request.GET.get("max_distance")
    min_distance = request.GET.get("min_distance")

    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    key_raw = f"{protein}|{chain}|{type_}"
    cache_key = "heatmap:" + hashlib.sha256(key_raw.encode("utf-8")).hexdigest()

    cached_png = cache.get(cache_key)
    if cached_png:
        logger.debug("heatmap cache hit: %s", cache_key)
        return HttpResponse(cached_png, content_type="image/png")

    logger.debug("heatmap cache miss: %s", cache_key)

    df = stat_aa_heatmap_df(
        protein=protein,
        chain=chain,
        type=type_,
        max_distance=max_distance,
        min_distance=min_distance
    )

    if df.empty or "aa1" not in df.columns:
        return JsonResponse({"error": "no_data", "detail": "No distance data found for this protein/chain/type."}, status=404)

    imgs = make_heat_maps(df, protein=protein, threshold=8)
    png = imgs["distance_png"]

    cache.set(cache_key, png, timeout=60 * 60 * 24)

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
        return HttpResponse(cached, content_type="image/png")

    lengths = stat_segment_lengths(protein, chain)
    if not lengths:
        return JsonResponse({"error": "no_data"}, status=404)

    png = make_segment_hist_png(lengths, protein=protein)
    cache.set(cache_key, png, timeout=60 * 60 * 24)

    return HttpResponse(png, content_type="image/png")


def global_stats(request):
    cached = cache.get("global_stats")
    if cached:
        return JsonResponse(cached)

    try:
        # protein/AA agregatne statistike — dužina po lancu
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

        # distribucija sekundarnih struktura
        ss_rows = run_query("""
            MATCH (a:AminoAcid)
            WHERE a.ss IS NOT NULL
            RETURN a.ss AS ss, count(a) AS n
            ORDER BY n DESC
        """)

        # top 5 aminokiselina u celoj bazi
        top_aa_rows = run_query("""
            MATCH (a:AminoAcid)
            WITH a.name AS name, count(a) AS n
            ORDER BY n DESC
            LIMIT 5
            RETURN name, n
        """)

        # statistike rastojanja po tipu — jedan sken svih DISTANCE
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
