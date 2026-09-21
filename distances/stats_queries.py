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

    # oba niza iz jednog upita, front ih spaja po poziciji pa redosled mora
    # biti isti; kod vise lanaca ORDER BY samo po index-u nije dovoljan
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


def stat_aa_heatmap_df(protein, chain, type, max_distance, min_distance, aminonames=None):

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND ($chain IS NULL OR (a1.chain = $chain AND a2.chain = $chain))
          AND ($aminonames IS NULL OR (a1.name IN $aminonames AND a2.name IN $aminonames))
        RETURN a1.chain AS ch1, a1.name AS aa1, a1.index AS idx1,
               a2.chain AS ch2, a2.name AS aa2, a2.index AS idx2,
               d.value AS distance
        ORDER BY a1.chain, a1.index, a2.chain, a2.index ASC
    """

    protein = protein.strip().upper()
    rows = run_query(
        cypher,
        protein=protein,
        type=type,
        chain=chain or None,
        aminonames=aminonames or None,
    )
    df = pd.DataFrame(rows)
    logger.debug("heatmap df columns: %s", df.columns.tolist())

    return df


def stat_ss_pairs(protein, chain, type, max_distance, min_distance):
    type_ = type if (type and type != "prazno") else "caca"

    cypher = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]->(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
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


# zajednicki WHERE za statistiku izabranog opsega, prati filtere iz
# query_builder.build_graph3d_protein da tabele opisu isti skup parova kao 3D prikaz
# granice rastojanja se dodaju posebno, jer se procenat racuna bez njih
_RANGE_BASE_WHERE = """
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND ($chain IS NULL OR (a1.chain = $chain AND a2.chain = $chain))
          AND ($aminonames IS NULL OR (a1.name IN $aminonames AND a2.name IN $aminonames))
          AND ($type IS NULL OR d.type = $type)
          AND ($ss IS NULL OR (a1.ss = $ss AND a2.ss = $ss))
          AND ($k IS NULL OR abs(a1.index - a2.index) >= $k)
"""

_RANGE_LIMITS_WHERE = """
          AND ($max_distance IS NULL OR d.value <= $max_distance)
          AND ($min_distance IS NULL OR d.value >= $min_distance)
"""


def stat_range_summary(protein, chain, type, max_distance, min_distance,
                       ss=None, k=None, aminonames=None):
    """Vraca statistiku rastojanja samo za parove koji prolaze trenutno izabrane filtere."""
    protein = protein.strip().upper()

    params = {
        "protein": protein,
        "chain": chain or None,
        "aminonames": aminonames or None,
        "type": type if (type and type != "prazno") else None,
        "ss": ss if (ss and ss not in ("prazno", "/")) else None,
        "k": int(k) if (k not in (None, "", "prazno")) else None,
        "max_distance": float(max_distance) if (max_distance not in (None, "", "prazno")) else None,
        "min_distance": float(min_distance) if (min_distance not in (None, "", "prazno")) else None,
    }

    # u jednom prolazu i ukupan broj parova bez granica i statistika unutar
    # granica; collect preskace null, pa u vals ostaju samo vrednosti iz opsega
    cypher_summary = """
        MATCH (a1:AminoAcid)-[d:DISTANCE]->(a2:AminoAcid)
    """ + _RANGE_BASE_WHERE + """
        WITH d.value AS v,
             (($max_distance IS NULL OR d.value <= $max_distance)
              AND ($min_distance IS NULL OR d.value >= $min_distance)) AS in_range
        WITH count(v) AS total_base, collect(CASE WHEN in_range THEN v END) AS vals
        UNWIND (CASE WHEN size(vals) = 0 THEN [null] ELSE vals END) AS dist
        RETURN total_base,
               count(dist)                            AS count,
               round(avg(dist), 3)                    AS mean,
               round(stdev(dist), 3)                  AS std,
               round(percentileCont(dist, 0.25), 3)   AS p25,
               round(percentileCont(dist, 0.50), 3)   AS median,
               round(percentileCont(dist, 0.75), 3)   AS p75,
               round(min(dist), 3)                    AS min,
               round(max(dist), 3)                    AS max
    """

    cypher_by_type = """
        MATCH (a1:AminoAcid)-[d:DISTANCE]->(a2:AminoAcid)
    """ + _RANGE_BASE_WHERE + _RANGE_LIMITS_WHERE + """
        WITH d.type AS type, d.value AS v
        RETURN type,
               count(v)                            AS count,
               round(avg(v), 3)                    AS mean,
               round(percentileCont(v, 0.50), 3)   AS median,
               round(min(v), 3)                    AS min,
               round(max(v), 3)                    AS max,
               round(stdev(v), 3)                  AS std
        ORDER BY count DESC
    """

    # par se normalizuje kao u stat_ss_pairs (H E i E H su isti par)
    cypher_by_ss_pair = """
        MATCH (a1:AminoAcid)-[d:DISTANCE]->(a2:AminoAcid)
    """ + _RANGE_BASE_WHERE + _RANGE_LIMITS_WHERE + """
          AND a1.ss IS NOT NULL AND a2.ss IS NOT NULL
        WITH
          CASE WHEN a1.ss <= a2.ss
               THEN a1.ss + ' ' + a2.ss
               ELSE a2.ss + ' ' + a1.ss
          END AS ss_pair,
          d.value AS v
        RETURN ss_pair,
               count(v)                            AS count,
               round(avg(v), 3)                    AS mean,
               round(percentileCont(v, 0.50), 3)   AS median,
               round(min(v), 3)                    AS min,
               round(max(v), 3)                    AS max,
               round(stdev(v), 3)                  AS std
        ORDER BY count DESC
    """

    rows = run_query(cypher_summary, **params)
    summary = rows[0] if rows else {}
    total_base = summary.get("total_base") or 0
    count = summary.get("count") or 0

    filters = {
        "type": params["type"],
        "chain": params["chain"],
        "ss": params["ss"],
        "k": params["k"],
        "max_distance": params["max_distance"],
        "min_distance": params["min_distance"],
        "aminonames": params["aminonames"],
    }

    if count == 0:
        return {
            "empty": True,
            "filters": filters,
            "total_base": total_base,
            "count": 0,
            "by_type": [],
            "by_ss_pair": [],
        }

    ss_pairs = run_query(cypher_by_ss_pair, **params)
    # udeo se racuna u odnosu na parove sa poznatom ss oznakom, a ne na count,
    # jer aminokiseline bez ss ne ulaze u tabelu
    ss_total = sum(r["count"] for r in ss_pairs) or 0
    for r in ss_pairs:
        r["pct"] = round(r["count"] * 100 / ss_total, 1) if ss_total else None

    return {
        "empty": False,
        "filters": filters,
        "total_base": total_base,
        "count": count,
        "pct_of_base": round(count * 100 / total_base, 1) if total_base else None,
        "mean": summary.get("mean"),
        "std": summary.get("std"),
        "p25": summary.get("p25"),
        "median": summary.get("median"),
        "p75": summary.get("p75"),
        "min": summary.get("min"),
        "max": summary.get("max"),
        "by_type": run_query(cypher_by_type, **params),
        "by_ss_pair": ss_pairs,
        "ss_pair_total": ss_total,
    }


def stat_outlier_ss(protein, chain, type, max_distance, min_distance):
    type_ = type if (type and type != "prazno") else "caca"

    # veza se hvata neusmereno, da prosek po reziduumu obuhvati sve partnere,
    # a ne samo one sa visim indeksom
    cypher_ss = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]-(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index <> a2.index
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

    # neusmerena veza iz istog razloga, inace bi najizolovaniji bili reziduumi
    # sa pocetka lanca
    cypher_top = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]-(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index <> a2.index
          AND ($chain IS NULL OR a1.chain = $chain)
        WITH a1.name AS aa, a1.index AS index, a1.chain AS chain,
             a1.ss AS ss, avg(d.value) AS avg_dist
        ORDER BY avg_dist DESC
        LIMIT 15
        RETURN aa, index, chain, ss, round(avg_dist, 3) AS avg_dist
    """

    # obrnuto od prethodnog upita: najmanje prosecno rastojanje imaju
    # aminokiseline iz jezgra proteina
    cypher_closest = """
        MATCH (a1:AminoAcid)-[d:DISTANCE{type:$type}]-(a2:AminoAcid)
        WHERE a1.protein = $protein AND a2.protein = $protein
          AND a1.index <> a2.index
          AND ($chain IS NULL OR a1.chain = $chain)
        WITH a1.name AS aa, a1.index AS index, a1.chain AS chain,
             a1.ss AS ss, avg(d.value) AS avg_dist
        ORDER BY avg_dist ASC
        LIMIT 15
        RETURN aa, index, chain, ss, round(avg_dist, 3) AS avg_dist
    """

    protein = protein.strip().upper()
    params = dict(protein=protein, type=type_, chain=chain or None)

    return {
        "type": type_,
        "ss_stats": run_query(cypher_ss, **params),
        "top_residues": run_query(cypher_top, **params),
        "closest_residues": run_query(cypher_closest, **params),
    }


def _segment_table(protein, chain):
    """Vraca tabelu segmenata, po jedan red za svaki neprekinuti niz iste SS klase."""
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
        return pd.DataFrame(columns=["chain", "seg_id", "ss", "length"])

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

    return (
        df.groupby(["chain", "seg_id"])
          .agg(ss=("ss", "first"), length=("index", "size"))
          .reset_index()
    )


def stat_segment_lengths(protein, chain):
    segments = _segment_table(protein, chain)
    if segments.empty:
        return []

    return segments["length"].tolist()


def stat_segment_lengths_by_ss(protein, chain):
    """Grupise duzine segmenata po SS oznaci, za histograme po SS tipu."""
    segments = _segment_table(protein, chain)
    if segments.empty:
        return {}

    return {
        ss: group["length"].tolist()
        for ss, group in segments.groupby("ss")
    }
