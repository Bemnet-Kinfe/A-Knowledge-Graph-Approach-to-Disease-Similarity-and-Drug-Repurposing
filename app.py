
from flask import Flask, jsonify, render_template, request
from Services.data_loader import load_primekg
from Services.analysis import (
    disease_profile,
    suggest_diseases,
    disease_similarity,
    repurposing_candidates,
    dataset_summary,)

app = Flask(__name__)
DATA = load_primekg()

@app.route("/")
def overview_page():
    # Overview and dataset summary
    summary = DATA["summary"]
    return render_template(
        "home.html",
        cleaning_stats=summary.get("cleaning_stats", {})
    )

@app.route("/neighborhood")
def neighborhood_page():
    # Disease neighborhood exploration
    return render_template("index.html") 

@app.route("/similarity")
def similarity_page():
    # Disease–disease Jaccard similarity
    return render_template("similarity.html")

@app.route("/repurposing")
def repurposing_page():
    # Drug repurposing candidates
    return render_template("repurposing.html")

@app.route("/api/summary")
def api_summary():
    kind = request.args.get("kind", "node_types")
    summary = DATA["summary"]

    if kind == "node_types":
        counts = summary["node_type_counts"]
        total = summary["total_nodes"]
        rows = []
        for t, c in counts.items():
            pct = c / total if total > 0 else 0.0
            rows.append({"type": t, "count": c, "pct": pct})
        rows.sort(key=lambda r: r["count"], reverse=True)
        return jsonify({
            "kind": "node_types",
            "total_nodes": total,
            "rows": rows
        })

    elif kind == "top_drugs":
        rows = summary["top_drugs"][:15]
        return jsonify({
            "kind": "top_drugs",
            "rows": [
                {"drug_name": r["drug_name"], "degree": r["degree"]}
                for r in rows
            ],
        })

    elif kind == "top_diseases_drug_neighbors":
        rows = summary["disease_drug_neighbors"][:15]
        return jsonify({
            "kind": "top_diseases_drug_neighbors",
            "rows": [
                {"disease_name": r["disease_name"], "n_drugs": r["n_drugs"]}
                for r in rows
            ],
        })

    elif kind == "example_relations":
        rows = summary.get("example_relations", [])
        rows = rows[:10]
        return jsonify({
            "kind": "example_relations",
            "rows": rows
        })

    else:
        return jsonify({"error": f"Unknown summary kind: {kind}"}), 400
    
@app.route("/api/disease_profile")
def api_disease_profile():
    disease = request.args.get("disease", "").strip()
    if not disease:
        return jsonify({"error": "Please provide a disease name."}), 400
    try:
        result = disease_profile(disease, DATA)
        return jsonify(result)
    except KeyError:
        candidates = suggest_diseases(disease, DATA)
        return jsonify({
            "error": f"We couldn't find a disease matching '{disease}' in PrimeKG.",
            "suggestions": candidates,
        }), 404
    except Exception as e:
        return jsonify({"error": "Unexpected error", "details": str(e)}), 500
        
@app.route("/api/disease_similarity")
def api_disease_similarity():
    disease = request.args.get("disease", "")
    neighbor_type = request.args.get("neighbor_type", "gene/protein")
    k_raw = request.args.get("k", "10")

    try:
        k = int(k_raw)
    except ValueError:
        k = 10

    try:
        result = disease_similarity(disease, neighbor_type, k, DATA)
        return jsonify(result)
    except KeyError as e:
        suggestions = suggest_diseases(disease, DATA, max_suggestions=8)
        return jsonify(
            {
                "error": str(e),
                "suggestions": suggestions,
            }
        ), 404
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": "Unexpected error", "details": str(e)}), 500


@app.route("/api/repurposing")
def api_repurposing():
    base_disease = request.args.get("disease", "")
    similar_ids_raw = request.args.get("similar_ids", "")

    similar_ids = [s for s in similar_ids_raw.split(",") if s]

    try:
        result = repurposing_candidates(base_disease, similar_ids, DATA)
        return jsonify(result)
    except KeyError as e:
        suggestions = suggest_diseases(base_disease, DATA)
        return jsonify({"error": str(e), "suggestions": suggestions}), 404
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": "Unexpected error", "details": str(e)}), 500


if __name__ == "__main__":
    app.run(debug=True)

