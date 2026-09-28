import pandas as pd
from pathlib import Path

NODE_ID_COL   = "node_id"
NODE_NAME_COL = "node_name"
NODE_TYPE_COL = "node_type"

EDGE_SRC_COL  = "x_id"
EDGE_DST_COL  = "y_id"
EDGE_REL_COL  = "relation"

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"


def load_primekg():
    # read the originnal CSV files
    kg_raw = pd.read_csv(DATA_DIR / "kg.csv")
    nodes_raw = pd.read_csv(DATA_DIR / "nodes.csv")
    disease_features = pd.read_csv(DATA_DIR / "disease_features.csv")
    drug_features = pd.read_csv(DATA_DIR / "drug_features.csv")

    # standardize column names to lowercase
    kg_raw.columns = [c.lower() for c in kg_raw.columns]
    nodes_raw.columns = [c.lower() for c in nodes_raw.columns]
    disease_features.columns = [c.lower() for c in disease_features.columns]
    drug_features.columns = [c.lower() for c in drug_features.columns]

    # remove duplicate edges if any
    raw_edge_rows = int(kg_raw.shape[0])
    kg = kg_raw.drop_duplicates()
    duplicate_edges_removed = raw_edge_rows - int(kg.shape[0])

    # remove duplicate nodes if any
    raw_node_rows = int(nodes_raw.shape[0])

    # remove fully identical rows first
    nodes_no_exact_dups = nodes_raw.drop_duplicates()
    exact_duplicate_node_rows_removed = raw_node_rows - int(nodes_no_exact_dups.shape[0])

    # then enforce uniqueness by node_id
    nodes_before_id_dedup = int(nodes_no_exact_dups.shape[0])
    nodes = nodes_no_exact_dups.drop_duplicates(subset=[NODE_ID_COL])
    duplicate_node_ids_removed = nodes_before_id_dedup - int(nodes.shape[0])

    # See which node IDs were duplicated and how many times, by type
    dup_mask = nodes_raw.duplicated(subset=[NODE_ID_COL], keep=False)
    dup_by_type_series = nodes_raw[dup_mask][NODE_TYPE_COL].value_counts(dropna=False)
    duplicate_node_ids_by_type = {str(t): int(v) for t, v in dup_by_type_series.items()}

    # node_info lookup structure
    node_info = {}
    for _, row in nodes.iterrows():
        nid = row[NODE_ID_COL]
        node_info[nid] = {
            "id": nid,
            "name": row.get(NODE_NAME_COL, str(nid)),
            "type": row.get(NODE_TYPE_COL, "unknown"),
        }

    disease_name_to_id = {}
    disease_id_to_name = {}
    for _, row in nodes[nodes[NODE_TYPE_COL] == "disease"].iterrows():
        nid = row[NODE_ID_COL]
        name = str(row.get(NODE_NAME_COL, nid))
        disease_name_to_id[name.lower()] = nid
        disease_id_to_name[nid] = name

    drug_id_to_name = {}
    for _, row in nodes[nodes[NODE_TYPE_COL] == "drug"].iterrows():
        nid = row[NODE_ID_COL]
        name = str(row.get(NODE_NAME_COL, nid))
        drug_id_to_name[nid] = name

    # precompute disease neighbors by type (genes/proteins, effects/phenotypes)
    disease_neighbors_by_type = {
        "gene/protein": {},
        "effect/phenotype": {},
    }

    for _, row in kg.iterrows():
        src = row[EDGE_SRC_COL]
        dst = row[EDGE_DST_COL]

        src_info = node_info.get(src)
        dst_info = node_info.get(dst)
        if not src_info or not dst_info:
            continue

        pairs = []
        if src_info["type"] == "disease":
            pairs.append((src, dst, dst_info["type"]))
        if dst_info["type"] == "disease":
            pairs.append((dst, src, src_info["type"]))

        for disease_id, neighbor_id, n_type in pairs:
            if n_type in disease_neighbors_by_type:
                s = disease_neighbors_by_type[n_type].setdefault(disease_id, set())
                s.add(neighbor_id)


    # summary statistics

    # 1) Node type distribution
    type_counts_series = nodes[NODE_TYPE_COL].value_counts(dropna=False)
    node_type_counts = {str(t): int(cnt) for t, cnt in type_counts_series.items()}
    total_nodes = int(nodes.shape[0])

    # 2) Degree per node
    degree_counts = (
        pd.concat([kg[EDGE_SRC_COL], kg[EDGE_DST_COL]])
        .value_counts()
        .to_dict()
    )

    # Top drugs by degree 
    top_drugs = []
    drug_nodes = nodes[nodes[NODE_TYPE_COL] == "drug"][[NODE_ID_COL, NODE_NAME_COL]]
    for _, row in drug_nodes.iterrows():
        nid = row[NODE_ID_COL]
        deg = int(degree_counts.get(nid, 0))
        top_drugs.append(
            {
                "drug_id": str(nid),
                "drug_name": str(row[NODE_NAME_COL]),
                "degree": deg,
            }
        )
    top_drugs.sort(key=lambda r: r["degree"], reverse=True)

    # 3) Diseases with most drug neighbors
    disease_to_drugs = {}
    for _, row in kg.iterrows():
        src = row[EDGE_SRC_COL]
        dst = row[EDGE_DST_COL]
        src_info = node_info.get(src)
        dst_info = node_info.get(dst)
        if not src_info or not dst_info:
            continue

        # disease-drug pairs 
        if src_info["type"] == "disease" and dst_info["type"] == "drug":
            disease_to_drugs.setdefault(src, set()).add(dst)
        if dst_info["type"] == "disease" and src_info["type"] == "drug":
            disease_to_drugs.setdefault(dst, set()).add(src)

    disease_drug_neighbors = []
    for did, drug_set in disease_to_drugs.items():
        disease_drug_neighbors.append(
            {
                "disease_id": str(did),
                "disease_name": disease_id_to_name.get(did, str(did)),
                "n_drugs": len(drug_set),
            }
        )
    disease_drug_neighbors.sort(key=lambda r: r["n_drugs"], reverse=True)

    # 4) Example relation table
    example_rows = []
    max_examples = 10  

    for _, row in kg.iterrows():
        src = row[EDGE_SRC_COL]
        dst = row[EDGE_DST_COL]
        rel = row.get(EDGE_REL_COL, "")

        src_info = node_info.get(src)
        dst_info = node_info.get(dst)
        if not src_info or not dst_info:
            continue

        # We want disease connected to something that is NOT a disease
        pairs = []
        if src_info["type"] == "disease" and dst_info["type"] != "disease":
            pairs.append((src_info, dst_info))
        if dst_info["type"] == "disease" and src_info["type"] != "disease":
            pairs.append((dst_info, src_info))

        for disease_node, neighbor_node in pairs:
            example_rows.append(
                {
                    "disease_name": disease_node["name"],
                    "disease_type": disease_node["type"],
                    "relation": rel,
                    "neighbor_name": neighbor_node["name"],
                    "neighbor_type": neighbor_node["type"],
                }
            )
            if len(example_rows) >= max_examples:
                break
        if len(example_rows) >= max_examples:
            break

    example_relations = example_rows

    summary = {
        "node_type_counts": node_type_counts,
        "total_nodes": total_nodes,
        "top_drugs": top_drugs,
        "disease_drug_neighbors": disease_drug_neighbors,
        "example_relations": example_relations,
        "cleaning_stats": {
            "raw_node_rows": raw_node_rows,
            "raw_edge_rows": raw_edge_rows,
            "duplicate_node_ids_removed": duplicate_node_ids_removed,
            "duplicate_node_ids_by_type": duplicate_node_ids_by_type,
        },
    }


    return {
        "kg": kg,
        "nodes": nodes,
        "disease_features": disease_features,
        "drug_features": drug_features,
        "node_info": node_info,
        "disease_name_to_id": disease_name_to_id,
        "disease_id_to_name": disease_id_to_name,
        "drug_id_to_name": drug_id_to_name,
        "disease_neighbors_by_type": disease_neighbors_by_type,
        "summary": summary,
    }
