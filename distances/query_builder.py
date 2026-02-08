def build_graph3d_protein(p):
    params = {"protein": p["protein"]}

    # node only grana - ako nisu zadati filteri 
    if (
        p["aminoname"] == "" and
        p["max_distance"] == "" and
        p["min_distance"] == "" and
        p["type"] == "prazno" 
    ):
        q = [
            "MATCH (a:AminoAcid)",
            "WHERE a.protein = $protein AND a.ca_coordinates IS NOT NULL",
            """
            RETURN
            a.protein AS protein,
            a.chain AS chain,
            a.index AS idx,
            a.name AS aa,
            a.ss AS ss,
            a.ca_coordinates.x AS x,
            a.ca_coordinates.y AS y,
            a.ca_coordinates.z AS z
            """
        ]
        return "\n".join(q), params

    q = [
        "MATCH (a1:AminoAcid)-[r:DISTANCE]->(a2:AminoAcid)",
        "WHERE a1.protein = $protein",
        "AND a1.index < a2.index"

    ]

    if p.get("chain") and p["chain"] != "prazno":
        q.append("AND a1.chain = $chain AND a2.chain = $chain")
        params["chain"] = p["chain"]

    if p["aminoname"] != "":
        q.append("AND a1.name = $aminoname1 AND a2.name = $aminoname2")
        params["aminoname1"] = p["aminoname"]
        params["aminoname2"] = p["aminoname"]

    if p["type"] != "prazno":
        q.append("AND r.type = $type")
        params["type"] = p["type"]

    if p["max_distance"]:
        q.append("AND r.value <= $max_distance")
        params["max_distance"] = float(p["max_distance"])

    if p["min_distance"]:
        q.append("AND r.value >= $min_distance")
        params["min_distance"] = float(p["min_distance"])

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