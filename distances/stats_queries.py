import logging

import pandas as pd
from .neo4j_client import run_query
from .protein_analysis import make_sspair_stats

logger = logging.getLogger(__name__)


STANDARD_AAS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def stat_aa_composition(protein, chain, type, max_distance, min_distance):

    cypher = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
          AND ($chain IS NULL OR aa.chain = $chain)
        WITH aa.name AS name, count(*) AS cnt
        WITH collect({name: name, cnt: cnt}) AS rows, sum(cnt) AS total
        UNWIND rows AS row
        RETURN row.name AS name, row.cnt AS cnt,
               round(toFloat(row.cnt) * 100 / total, 1) AS pct,
               total
        ORDER BY row.cnt DESC
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, chain=chain or None)

    total = rows[0]["total"] if rows else 0
    present = [{"name": r["name"], "cnt": r["cnt"], "pct": r["pct"]} for r in rows]

    present_names = {r["name"] for r in rows}
    missing = sorted(STANDARD_AAS - present_names)

    return {
        "protein": protein,
        "total": total,
        "aa_composition": present,
        "missing": missing,
    }


def stat_ss_distribution(protein, chain, type, max_distance, min_distance):
    cypher = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
          AND ($chain IS NULL OR aa.chain = $chain)
          AND aa.ss IS NOT NULL
        WITH aa.ss AS ss, count(*) AS cnt
        WITH collect({ss: ss, cnt: cnt}) AS rows, sum(cnt) AS total
        UNWIND rows AS row
        RETURN row.ss AS ss, row.cnt AS cnt,
               round(toFloat(row.cnt) * 100 / total, 1) AS pct,
               total
        ORDER BY row.cnt DESC
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, chain=chain or None)

    total = rows[0]["total"] if rows else 0
    distribution = [{"ss": r["ss"], "cnt": r["cnt"], "pct": r["pct"]} for r in rows]

    return {
        "protein": protein,
        "total": total,
        "ss_distribution": distribution,
    }


def stat_aa_seq(protein, chain, type, max_distance, min_distance):

    # Jedan upit za oba niza — front ih spaja po poziciji (SEQ_CACHE[i]),
    # pa poravnanje mora biti garantovano istim skupom redova. Sa dva
    # odvojena upita i ORDER BY samo po index-u, kod višelančanih proteina
    # izjednačeni indeksi (A:1, B:1, ...) nisu deterministički poređani,
    # što bi tiho pomerilo SS traku u odnosu na sekvencu.
    cypher = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
          AND ($chain IS NULL OR aa.chain = $chain)
        RETURN aa.name AS name, aa.ss AS ss
        ORDER BY aa.chain, aa.index
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, chain=chain or None)

    return {
        "protein": protein,
        "aa_seq": [{"name": r["name"]} for r in rows],
        "aa_ss_seq": [{"ss": r["ss"]} for r in rows],
    }


def stat_aa_heatmap_df(protein, chain, type, max_distance, min_distance):

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND ($chain IS NULL OR (a1.chain = $chain AND a2.chain = $chain))
        RETURN a1.name AS aa1, a1.index AS idx1, a2.name AS aa2, a2.index AS idx2, d.value AS distance
        ORDER BY a1.index, a2.index ASC
    """

    protein = protein.strip().upper()
    rows = run_query(
        cypher,
        protein=protein,
        type=type,
        chain=chain or None,
    )
    df = pd.DataFrame(rows)
    logger.debug("heatmap df columns: %s", df.columns.tolist())

    return df


def stat_ss_pairs(protein, chain, type, max_distance, min_distance):
    type_ = type if (type and type != "prazno") else "caca"

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index < a2.index
          AND a1.ss IS NOT NULL AND a2.ss IS NOT NULL
          AND ($chain IS NULL OR a1.chain = $chain)
        WITH
          CASE WHEN a1.ss <= a2.ss
               THEN a1.ss + ' ' + a2.ss
               ELSE a2.ss + ' ' + a1.ss
          END AS ss_pair,
          d.value AS dist
        RETURN
          ss_pair,
          count(dist)                       AS count,
          round(avg(dist), 3)               AS mean,
          round(min(dist), 3)               AS min,
          round(max(dist), 3)               AS max,
          round(stdev(dist), 3)             AS std,
          round(percentileCont(dist, 0.5), 3) AS median
        ORDER BY count DESC
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, type=type_, chain=chain or None)

    return {"type": type_, "ss_pairs": rows}


def stat_distance_summary(protein, chain, type, max_distance, min_distance):
    type_ = type if (type and type != "prazno") else "caca"

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index < a2.index
          AND ($chain IS NULL OR a1.chain = $chain)
        WITH collect(d.value) AS vals,
             count(d.value)                          AS total,
             avg(d.value)                            AS mean,
             stdev(d.value)                          AS std,
             percentileCont(d.value, 0.10)           AS p10,
             percentileCont(d.value, 0.25)           AS p25,
             percentileCont(d.value, 0.50)           AS median,
             percentileCont(d.value, 0.75)           AS p75,
             percentileCont(d.value, 0.90)           AS p90,
             percentileCont(d.value, 0.95)           AS p95,
             percentileCont(d.value, 0.99)           AS p99,
             min(d.value)                            AS min_val,
             max(d.value)                            AS max_val,
             sum(CASE WHEN d.value <  5.0 THEN 1 ELSE 0 END) AS cnt_lt_5,
             sum(CASE WHEN d.value <  8.0 THEN 1 ELSE 0 END) AS cnt_lt_8,
             sum(CASE WHEN d.value < 10.0 THEN 1 ELSE 0 END) AS cnt_lt_10,
             sum(CASE WHEN d.value < 15.0 THEN 1 ELSE 0 END) AS cnt_lt_15
        RETURN
          total                                              AS count,
          round(mean,   3)                                   AS mean,
          round(std,    3)                                   AS std,
          round(p10,    3)                                   AS p10,
          round(p25,    3)                                   AS p25,
          round(median, 3)                                   AS median,
          round(p75,    3)                                   AS p75,
          round(p90,    3)                                   AS p90,
          round(p95,    3)                                   AS p95,
          round(p99,    3)                                   AS p99,
          round(min_val, 3)                                  AS min,
          round(max_val, 3)                                  AS max,
          round(toFloat(cnt_lt_5)  * 100 / total, 1)        AS pct_lt_5,
          round(toFloat(cnt_lt_8)  * 100 / total, 1)        AS pct_lt_8,
          round(toFloat(cnt_lt_10) * 100 / total, 1)        AS pct_lt_10,
          round(toFloat(cnt_lt_15) * 100 / total, 1)        AS pct_lt_15
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, type=type_, chain=chain or None)

    if not rows or rows[0]["count"] == 0:
        return {"type": type_, "empty": True}

    return {"type": type_, **rows[0]}


def stat_outlier_ss(protein, chain, type, max_distance, min_distance):
    type_ = type if (type and type != "prazno") else "caca"

    cypher_ss = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index < a2.index
          AND ($chain IS NULL OR a1.chain = $chain)
          AND a1.ss IS NOT NULL
        WITH a1.ss AS ss, a1.index AS idx, avg(d.value) AS avg_dist
        WITH ss,
             count(*)                                     AS n_rezidua,
             round(avg(avg_dist), 3)                      AS mean,
             round(percentileCont(avg_dist, 0.5), 3)      AS median,
             round(percentileCont(avg_dist, 0.9), 3)      AS p90,
             round(max(avg_dist), 3)                      AS max
        RETURN ss, n_rezidua, mean, median, p90, max
        ORDER BY mean DESC
    """

    cypher_top = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index < a2.index
          AND ($chain IS NULL OR a1.chain = $chain)
        WITH a1.name AS aa, a1.index AS index, a1.chain AS chain,
             a1.ss AS ss, avg(d.value) AS avg_dist
        ORDER BY avg_dist DESC
        LIMIT 15
        RETURN aa, index, chain, ss, round(avg_dist, 3) AS avg_dist
    """

    protein = protein.strip().upper()
    params = dict(protein=protein, type=type_, chain=chain or None)

    return {
        "type": type_,
        "ss_stats": run_query(cypher_ss, **params),
        "top_residues": run_query(cypher_top, **params),
    }


def stat_segment_lengths(protein, chain):
    cypher = """
        MATCH (aa:AminoAcid)
        WHERE aa.protein = $protein
          AND ($chain IS NULL OR aa.chain = $chain)
          AND aa.ss IS NOT NULL
        RETURN aa.chain AS chain, aa.index AS index, aa.ss AS ss
        ORDER BY aa.chain, aa.index
    """

    protein = protein.strip().upper()
    rows = run_query(cypher, protein=protein, chain=chain or None)
    if not rows:
        return []

    df = pd.DataFrame(rows)
    df["index"] = pd.to_numeric(df["index"], errors="coerce")
    df = df.dropna(subset=["index"]).sort_values(["chain", "index"])

    # nova grupa kad se promeni lanac, ss, ili index nije uzastopan
    df["new_seg"] = (
        (df["chain"] != df["chain"].shift()) |
        (df["ss"] != df["ss"].shift()) |
        (df["index"] != df["index"].shift() + 1)
    )
    df["seg_id"] = df["new_seg"].cumsum()

    return df.groupby(["chain", "seg_id"]).size().tolist()
