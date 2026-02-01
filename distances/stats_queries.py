import pandas as pd
from django.http import JsonResponse
from .neo4j_client import run_query
from .protein_analysis import make_heat_maps

def stat_aa_composition(protein, chain, type, max_distance, min_distance):

    cypher = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
        RETURN aa.name AS name, count(*) AS cnt
        ORDER BY cnt DESC
    """

    protein = protein.strip().upper()
    print("protein", protein)
    rows = run_query(cypher, protein=protein)

    # (opciono) dodaj procente
    total = sum(r["cnt"] for r in rows) if rows else 0
    for r in rows:
        r["pct"] = (r["cnt"] / total * 100) if total else 0

    return {
        "protein": protein,
        "total": total,
        "aa_composition": rows
    }


def stat_aa_seq(protein, chain, type, max_distance, min_distance):

    cypher1 = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
        RETURN aa.name AS name
        ORDER BY aa.index
    """

    cypher2 = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
        RETURN aa.ss AS ss
        ORDER BY aa.index 
    """

    protein = protein.strip().upper()
    rows1 = run_query(cypher1, protein=protein)
    rows2 = run_query(cypher2, protein=protein)


    return {
        "protein": protein,
        "aa_seq": rows1,
        "aa_ss_seq": rows2
    }


def stat_aa_heatmap_df(protein, chain, type, max_distance, min_distance):

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein=$protein AND a2.protein=$protein
        RETURN a1.name AS aa1, a1.index AS idx1, a2.name AS aa2, a2.index AS idx2, d.value AS distance
        ORDER BY a1.index, a2.index ASC
    """


    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, type=type)
    df = pd.DataFrame(rows)
    print("KOLONE", df.columns)


    return df