def _present(val):
    return val is not None and val != "" and val not in ("prazno", "/")


def build_graph3d_protein(p):
    params = {"protein": p["protein"]}

    has_ss = _present(p.get("ss"))
    has_k  = _present(p.get("k"))

    # grana samo sa cvorovima, bez edge filtera (k trazi edge granu)
    if (
        not p.get("aminonames") and
        not p.get("max_distance") and
        not p.get("min_distance") and
        not _present(p.get("type")) and
        not has_k
    ):
        q = [
            "MATCH (a:AminoAcid)",
            "WHERE a.protein = $protein AND a.ca_coordinates IS NOT NULL",
        ]

        if has_ss:
            q.append("AND a.ss = $ss")
            params["ss"] = p["ss"]

        if _present(p.get("chain")):
            q.append("AND a.chain = $chain")
            params["chain"] = p["chain"]

        q.append("""
        RETURN
            a.protein AS protein,
            a.chain AS chain,
            a.index AS idx,
            a.name AS aa,
            a.ss AS ss,
            a.ca_coordinates.x AS x,
            a.ca_coordinates.y AS y,
            a.ca_coordinates.z AS z
        """)
        return "\n".join(q), params

    # edge grana
    q = [
        "MATCH (a1:AminoAcid)-[r:DISTANCE]->(a2:AminoAcid)",
        "WHERE a1.protein = $protein",
        # bez Ca koordinata cvor se ne moze iscrtati: x/y/z bi bili null, pa bi
        # jedan NaN u zoomTo() oborio ceo prikaz
        "AND a1.ca_coordinates IS NOT NULL AND a2.ca_coordinates IS NOT NULL",
    ]

    if _present(p.get("chain")):
        q.append("AND a1.chain = $chain AND a2.chain = $chain")
        params["chain"] = p["chain"]

    if p.get("aminonames"):
        q.append("AND a1.name IN $aminonames AND a2.name IN $aminonames")
        params["aminonames"] = p["aminonames"]

    if _present(p.get("type")):
        q.append("AND r.type = $type")
        params["type"] = p["type"]

    if p.get("max_distance"):
        q.append("AND r.value <= $max_distance")
        params["max_distance"] = float(p["max_distance"])

    if p.get("min_distance"):
        q.append("AND r.value >= $min_distance")
        params["min_distance"] = float(p["min_distance"])

    if has_ss:
        q.append("AND a1.ss = $ss AND a2.ss = $ss")
        params["ss"] = p["ss"]

    if has_k:
        q.append("AND abs(a1.index - a2.index) >= $k")
        params["k"] = int(p["k"])

    q.append("""
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
    """)

    return "\n".join(q), params
