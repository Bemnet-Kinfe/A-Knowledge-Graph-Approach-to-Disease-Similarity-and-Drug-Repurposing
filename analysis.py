
from collections import Counter
from typing import Dict, Any, List, Set   

from Services.data_loader import (
    NODE_ID_COL,
    NODE_NAME_COL,
    NODE_TYPE_COL,
    EDGE_SRC_COL,
    EDGE_DST_COL,
    EDGE_REL_COL,
)

def dataset_summary(kind: str, data: Dict[str, Any]) -> Dict[str, Any]:
    kg = data["kg"]
    nodes = data["nodes"]
    node_info = data["node_info"]

    kind = kind.strip().lower()

    if kind == "node_types":
        counts = Counter(nodes[NODE_TYPE_COL].fillna("unknown"))
        total = int(len(nodes))
        rows = [
            {"type": t, "count": int(c), "pct": float(c) / total}
            for t, c in counts.most_common()
        ]
        return {
            "kind": "node_types",
            "total_nodes": total,
            "rows": rows,
        }

    elif kind == "top_drugs":
        # degree of each drug node 
        deg = Counter()
        for _, row in kg.iterrows():
            src = row[EDGE_SRC_COL]
            dst = row[EDGE_DST_COL]

            src_info = node_info.get(src)
            dst_info = node_info.get(dst)

            if src_info and src_info["type"] == "drug":
                deg[src] += 1
            if dst_info and dst_info["type"] == "drug":
                deg[dst] += 1

        top = deg.most_common(15)
        rows = []
        for drug_id, d in top:
            name = node_info.get(drug_id, {}).get("name", str(drug_id))
            rows.append(
                {
                    "drug_id": drug_id,
                    "drug_name": name,
                    "degree": int(d),
                }
            )
        return {
            "kind": "top_drugs",
            "rows": rows,
        }

    elif kind == "top_diseases_drug_neighbors":
        # number of unique drug neighbors per disease
        from collections import defaultdict

        disease_to_drugs: Dict[Any, Set[Any]] = defaultdict(set)

        for _, row in kg.iterrows():
            src = row[EDGE_SRC_COL]
            dst = row[EDGE_DST_COL]
            src_info = node_info.get(src)
            dst_info = node_info.get(dst)
            if not src_info or not dst_info:
                continue

            if src_info["type"] == "disease" and dst_info["type"] == "drug":
                disease_to_drugs[src].add(dst)
            if dst_info["type"] == "disease" and src_info["type"] == "drug":
                disease_to_drugs[dst].add(src)

        counts = [
            (did, len(drugs)) for did, drugs in disease_to_drugs.items()
        ]
        counts.sort(key=lambda x: x[1], reverse=True)
        top = counts[:15]

        rows = []
        for did, n_drugs in top:
            name = node_info.get(did, {}).get("name", str(did))
            rows.append(
                {
                    "disease_id": did,
                    "disease_name": name,
                    "n_drugs": int(n_drugs),
                }
            )
        return {
            "kind": "top_diseases_drug_neighbors",
            "rows": rows,
        }

    else:
        raise ValueError(f"Unknown summary kind: {kind}")
    
# disease neighborhood explorer   
def disease_profile(disease_name: str, data: Dict[str, Any]) -> Dict[str, Any]:
    kg = data["kg"]
    node_info = data["node_info"]
    disease_name_to_id = data["disease_name_to_id"]
    disease_id_to_name = data["disease_id_to_name"]

    key = disease_name.lower().strip()
    if key not in disease_name_to_id:
        raise KeyError(f"Disease '{disease_name}' not found in nodes table")

    disease_id = disease_name_to_id[key]

    mask = (kg[EDGE_SRC_COL] == disease_id) | (kg[EDGE_DST_COL] == disease_id)
    edges = kg[mask]

    neighbor_counts = Counter()
    neighbors_detail: List[Dict[str, Any]] = []

    for _, row in edges.iterrows():
        src = row[EDGE_SRC_COL]
        dst = row[EDGE_DST_COL]
        relation = row.get(EDGE_REL_COL, "")

        neighbor_id = dst if src == disease_id else src

        info = node_info.get(neighbor_id)
        if info is None:
            continue

        n_type = info["type"]
        n_name = info["name"]

        neighbor_counts[n_type] += 1

        neighbors_detail.append(
            {
                "id": neighbor_id,
                "name": n_name,
                "type": n_type,
                "relation": relation,
            }
        )

    return {
        "disease_id": disease_id,
        "disease_name": disease_id_to_name.get(disease_id, disease_name),
        "neighbor_counts": dict(neighbor_counts),
        "neighbors": neighbors_detail,
    }

# suggest disease names when the quiery is not found in the dataset
def suggest_diseases(query: str, data: Dict[str, Any], max_suggestions: int = 8) -> List[str]:
    nodes = data["nodes"]
    q = query.strip().lower()
    if not q:
        return []

    mask = (
        nodes[NODE_TYPE_COL].str.contains("disease", case=False, na=False)
        & nodes[NODE_NAME_COL].str.lower().str.contains(q, na=False)
    )

    suggestions = (
        nodes.loc[mask, NODE_NAME_COL]
        .dropna()
        .drop_duplicates()
        .sort_values()
        .tolist()
    )

    return suggestions[:max_suggestions]

# disease-disease similarity based on shared genes/phenotypes
def disease_similarity(
    disease_name: str,
    neighbor_type: str,
    k: int,
    data: Dict[str, Any],
) -> Dict[str, Any]:
    nodes = data["nodes"]
    disease_name_to_id = data["disease_name_to_id"]
    disease_id_to_name = data["disease_id_to_name"]
    neighbors_by_type = data["disease_neighbors_by_type"]

    neighbor_type = neighbor_type.strip().lower()
    if neighbor_type not in neighbors_by_type:
        raise ValueError(
            f"Neighbor type '{neighbor_type}' is not supported. "
            "Use 'gene/protein' or 'effect/phenotype'."
        )

    # Look up the query disease
    key = disease_name.lower().strip()
    if key not in disease_name_to_id:
        raise KeyError(f"Disease '{disease_name}' not found in nodes table")

    query_id = disease_name_to_id[key]
    query_name = disease_id_to_name.get(query_id, disease_name)

    # Precomputed neighbor sets for this type
    disease_to_neighbors: Dict[Any, Set[Any]] = neighbors_by_type[neighbor_type]

    # Neighbor set for the query disease
    query_neighbors = disease_to_neighbors.get(query_id, set())
    if not query_neighbors:
        raise ValueError(
            f"No neighbors of type '{neighbor_type}' found for '{query_name}'."
        )

    similarities: List[Dict[str, Any]] = []

    # Iterate over diseases that actually have neighbors of this type
    for other_id, other_neighbors in disease_to_neighbors.items():
        if other_id == query_id:
            continue
        if not other_neighbors:
            continue

        inter = query_neighbors & other_neighbors
        union = query_neighbors | other_neighbors
        if not union:
            continue

        jacc = len(inter) / len(union)
        if jacc <= 0:
            continue

        other_name = disease_id_to_name.get(
            other_id,
            str(other_id),
        )

        similarities.append(
            {
                "disease_id": other_id,
                "disease_name": other_name,
                "jaccard": jacc,
                "n_shared": len(inter),
                "n_query": len(query_neighbors),
                "n_other": len(other_neighbors),
            }
        )

    # Sort by similarity and keep top-k
    similarities.sort(key=lambda x: x["jaccard"], reverse=True)
    topk = similarities[:k]

    return {
        "query_disease_id": query_id,
        "query_disease_name": query_name,
        "neighbor_type": neighbor_type,
        "results": topk,
    }

# Drug repurposing candidates based on shared disease neighbors
def repurposing_candidates(
    base_disease_name: str,
    similar_disease_ids: List[str],
    data: Dict[str, Any],
) -> Dict[str, Any]:

    kg = data["kg"]
    node_info = data["node_info"]
    nodes = data["nodes"]
    disease_name_to_id = data["disease_name_to_id"]
    disease_id_to_name = data["disease_id_to_name"]
    drug_features = data.get("drug_features")  

    # Look up base disease ID 
    key = base_disease_name.lower().strip()
    if key not in disease_name_to_id:
        raise KeyError(f"Disease '{base_disease_name}' not found in nodes table")

    base_id = disease_name_to_id[key]
    base_name = disease_id_to_name.get(base_id, base_disease_name)

    # Drug neighbors of a disease 
    def get_drugs(did: Any) -> Set[Any]:
        mask = (kg[EDGE_SRC_COL] == did) | (kg[EDGE_DST_COL] == did)
        edges = kg[mask]

        neighbors: Set[Any] = set()
        for _, row in edges.iterrows():
            src = row[EDGE_SRC_COL]
            dst = row[EDGE_DST_COL]
            neighbor_id = dst if src == did else src

            info = node_info.get(neighbor_id)
            if info is None:
                continue
            if info["type"].lower() == "drug":
                neighbors.add(neighbor_id)
        return neighbors

    base_drugs = get_drugs(base_id)

    # Validate similar diseases list
    similar_disease_ids = [s for s in similar_disease_ids if s]
    if not similar_disease_ids:
        raise ValueError("No similar diseases provided for repurposing analysis.")

    support_map: Dict[Any, Set[Any]] = {}

    for other_id in similar_disease_ids:
        # Skip if accidentally including the base disease
        if other_id == base_id:
            continue

        other_drugs = get_drugs(other_id)
        if not other_drugs:
            continue

        for drug_id in other_drugs:
            # Exclude drugs already linked to the base disease
            if drug_id in base_drugs:
                continue
            support_map.setdefault(drug_id, set()).add(other_id)

    # Build candidate list 
    candidates: List[Dict[str, Any]] = []

    for drug_id, supporting_diseases in support_map.items():
        info = node_info.get(drug_id, {"name": str(drug_id), "type": "drug"})
        drug_name = info.get("name", str(drug_id))

        # description of the drug from the drug_features table, if available
        description = ""
        if drug_features is not None:
            try:
                rows = drug_features.loc[drug_features[NODE_ID_COL] == drug_id]
                if not rows.empty:
                    if "description" in rows.columns:
                        description = str(rows["description"].iloc[0])
                    elif "drug_description" in rows.columns:
                        description = str(rows["drug_description"].iloc[0])
            except Exception:
                # fail silently if schema is a bit different
                description = ""

        supporting_names = [
            disease_id_to_name.get(did, str(did)) for did in supporting_diseases
        ]

        candidates.append(
            {
                "drug_id": drug_id,
                "drug_name": drug_name,
                "n_support": len(supporting_diseases),
                "supporting_diseases": supporting_names,
                "description": description,
            }
        )

    # Sort in descendinng order for supporting diseases 
    candidates.sort(key=lambda x: (-x["n_support"], x["drug_name"]))

    return {
        "base_disease_id": base_id,
        "base_disease_name": base_name,
        "n_base_drugs": len(base_drugs),
        "n_candidates": len(candidates),
        "candidates": candidates,
    }

