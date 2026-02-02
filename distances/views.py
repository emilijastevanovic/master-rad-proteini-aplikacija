from django.shortcuts import render

from django.http import JsonResponse
from django.http import HttpResponse, JsonResponse
from django.core.cache import cache
from .neo4j_client import run_query
from .query_builder import  build_graph3d_protein
from .stats_queries import stat_aa_composition, stat_aa_seq, stat_aa_heatmap_df
from .protein_analysis import make_heat_maps
import hashlib
from neo4j.exceptions import Neo4jError

STAT_HANDLERS = {
    "aa_composition": stat_aa_composition,
    "sequence": stat_aa_seq,
    # "ss_track": stat_ss_track,
    # ...
}

def graph3d_view(request):
    protein = request.GET.get("protein")
    type_ = request.GET.get("type", "caca")
    max_distance = request.GET.get("max_distance")

    try:
        limit = int(request.GET.get("limit", 20000))
    except ValueError:
        limit = 20000
    limit = max(1, min(limit, 50000))  # hard cap

    q = """
    MATCH (a1:AminoAcid)-[r:DISTANCE]->(a2:AminoAcid)
    WHERE a1.protein = $protein AND r.type = $type
    """ 
    if max_distance:
        q += " AND r.value <= $max_distance"

    q += """
    WITH a1, a2, r
    ORDER BY r.value ASC
    LIMIT $limit
    RETURN
        a1.protein AS protein,
        a1.chain AS chain,
        a1.index AS index,
        a1.name AS aa,
        a1.ss AS ss1,
        a1.ca_coordinates.x AS x,
        a1.ca_coordinates.y AS y,
        a1.ca_coordinates.z AS z,
        a2.protein AS protein2,
        a2.chain AS chain2,
        a2.index AS index2,
        a2.name AS aa2,
        a2.ss AS ss2,
        a2.ca_coordinates.x AS x2,
        a2.ca_coordinates.y AS y2,
        a2.ca_coordinates.z AS z2,
        r.type AS type,
        r.value AS value
    """

    cy_params = {"protein": protein, "type": type_, "limit": limit}
    if max_distance:
        cy_params["max_distance"] = float(max_distance)

    try:
        rows = run_query(q, **cy_params)
    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)

    nodes, edges = {}, []

    for r in rows:
        a1_id = f"{r['protein']}:{r['chain']}:{r['index']}"
        a2_id = f"{r['protein2']}:{r['chain2']}:{r['index2']}"
        if a1_id not in nodes:
            nodes[a1_id] = {
                "id": a1_id,
                "protein": r["protein"],
                "chain": r["chain"],
                "index": r["index"],
                "aa": r["aa"],
                "ss": r.get("ss1"),
                "x": r.get("x"),
                "y": r.get("y"),
                "z": r.get("z"),
            }
        if a2_id not in nodes:
            nodes[a2_id] = {
                "id": a2_id,
                "protein": r["protein2"],
                "chain": r["chain2"],
                "index": r["index2"],
                "aa": r["aa2"],
                "ss": r.get("ss2"),
                "x": r.get("x2"),
                "y": r.get("y2"),
                "z": r.get("z2"),
            }

        edges.append({
            "source": a1_id,
            "target": a2_id,
            "type": r["type"],
            "value": r["value"],
        })

    return JsonResponse({
        "nodes": list(nodes.values()),
        "edges": edges
    })

def graph3d_protein(request):

    params = {
        "protein": request.GET.get("protein"),
        "aminoname": request.GET.get("aminoname"),
        "type": request.GET.get("type", "prazno"),
        "max_distance": request.GET.get("max_distance"),
        "min_distance": request.GET.get("min_distance"),
        "chain": request.GET.get("chain"),
    }

    nodes = {}
    edges = []

    # generiše Cypher i parametre
    q, cy_params = build_graph3d_protein(params)

    print("***************")
    print("UPIT")
    print(q)
    print("PARAMS:", cy_params)
    print("***************")

    # izvršavanje Neo4j upita
    try:
        rows = run_query(q, **cy_params)
    except (RuntimeError, Neo4jError) as e:
        return JsonResponse({"error": "neo4j_unavailable", "detail": str(e)}, status=503)
    print("ROW COUNT:", len(rows))

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
    chain = request.GET.get("chain")
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

def heatmap_distance(request):
    protein = request.GET.get("protein")
    chain = request.GET.get("chain")
    type_ = request.GET.get("type", "caca")
    max_distance = request.GET.get("max_distance")
    min_distance = request.GET.get("min_distance")

    if not protein:
        return JsonResponse({"error": "protein parameter is required"}, status=400)

    # CACHE KEY (po svim parametrima koji utiču na rezultat)
    key_raw = f"{protein}|{chain}|{type_}|thr=8"
    cache_key = "heatmap:" + hashlib.sha256(key_raw.encode("utf-8")).hexdigest()

    cached_png = cache.get(cache_key)
    if cached_png:
        print("HEATMAP CACHE HIT", cache_key)
        return HttpResponse(cached_png, content_type="image/png")

    print("HEATMAP CACHE MISS", cache_key)

    df = stat_aa_heatmap_df(
        protein=protein,
        chain=chain,
        type=type_,
        max_distance=max_distance,
        min_distance=min_distance
    )

    imgs = make_heat_maps(df, protein=protein, threshold=8)
    png = imgs["distance_png"]

    cache.set(cache_key, png, timeout=60 * 60 * 24)  # npr 24h

    return HttpResponse(png, content_type="image/png")


