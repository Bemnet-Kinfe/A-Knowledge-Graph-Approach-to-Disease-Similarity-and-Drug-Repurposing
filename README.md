# PrimeKG Explorer
### A Knowledge-Graph Approach to Disease Similarity and Drug Repurposing

**Author:** Bemnet Kinfe
**Course:** BIS 634 – Final Project (December 2025)

PrimeKG Explorer is an interactive web application for exploring the **Precision Medicine Knowledge Graph (PrimeKG)**. It lets users inspect a disease's neighborhood in the graph, find similar diseases based on shared genes or phenotypes, and surface candidate drugs for repurposing.

![PrimeKG Explorer overview page](images/fig16_overview_page.png)

---

## Table of Contents
1. [Dataset](#dataset)
2. [Data Cleaning & Preprocessing](#data-cleaning--preprocessing)
3. [Summary Statistics](#summary-statistics)
4. [Analyses](#analyses)
   - [Disease Neighborhood Explorer](#1-disease-neighborhood-explorer)
   - [Disease–Disease Similarity](#2-diseasedisease-similarity)
   - [Drug Repurposing Candidates](#3-drug-repurposing-candidates)
5. [Validation](#validation)
6. [API Reference](#api-reference)
7. [Web Interface](#web-interface)
8. [Getting Started](#getting-started)
9. [Challenges & Feedback Integration](#challenges--feedback-integration)
10. [Citation & License](#citation--license)

---

## Dataset

**PrimeKG** was developed by Chowdhury and Zitnik (2023) and integrates **20 high-quality biomedical databases and ontologies**. It links **17,000+ diseases** through **4+ million relationships** spanning drugs, genes/proteins, pathways, phenotypes, side effects, and clinical guideline text.

### Why PrimeKG?
PrimeKG represents biomedical knowledge as a connected graph rather than isolated tables. This makes it possible to ask questions that cross disciplinary boundaries, such as:
- Which diseases share molecular or phenotypic mechanisms?
- Where are there treatment gaps or repurposing opportunities?
- How do biological mechanisms relate to what appears in clinical data?

It mirrors how we naturally reason about biology, connecting molecules to mechanisms and mechanisms to patient care, at a scale no individual researcher could manage.

### Acquisition
The dataset was downloaded as CSV files from the [official PrimeKG project page](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) on the Harvard Dataverse. No API access or authentication was required.

### FAIR Assessment

| Principle | Assessment |
|---|---|
| **Findable** | Hosted on Harvard Dataverse with a DOI ([10.7910/DVN/IXA7BM](https://doi.org/10.7910/DVN/IXA7BM)) and rich metadata on scope, version, and provenance. |
| **Accessible** | Freely downloadable over HTTPS without authentication. A dedicated API or standardized query protocol would improve programmatic access. |
| **Interoperable** | Strong semantic interoperability via standard identifiers (e.g., MeSH for diseases, DrugBank IDs for drugs). Technical interoperability is limited by the CSV format; conversion to semantic web formats would ease linking with other knowledge graphs. |
| **Reusable** | Released under CC0 1.0 Public Domain Dedication, allowing unrestricted reuse and redistribution. |

---

## Data Cleaning & Preprocessing

The data consists of four CSV files (`kg.csv`, `nodes.csv`, `disease_features.csv`, `drug_features.csv`), loaded with a custom `load_primekg()` function.

1. **Column name standardization.** All column names were lowercased to avoid case-sensitivity issues across files.
2. **Duplicate removal.**
   - Exact duplicate edges were removed with `kg.drop_duplicates()`.
   - Duplicate nodes were removed based on `node_id`, so each node has a single definition.
   - The raw data contained **129,375 node rows** and **8,100,498 edge rows**. Cleaning removed **39,308 duplicate nodes**, most commonly diseases (12,782) and phenotypes (11,869), leaving **90,067 unique nodes**.
3. **Missing node information.** A `node_info` dictionary stores each node's ID, name, and type. Missing names default to the `node_id` string; missing types default to `"unknown"`. No artificial values are introduced.
4. **Precomputed disease–neighbor sets.** `disease_neighbors_by_type` maps each disease to its `gene/protein` and `effect/phenotype` neighbors. Similarity queries then reduce to fast set operations instead of rescanning the full edge table.
5. **Text normalization for search.** Disease names are stored in a lowercased lookup (`disease_name_to_id[name.lower()]`), and user input is lowercased and stripped before matching.

Overall, preprocessing avoids imputing biomedical relationships, keeping results transparent and free of artificial artifacts.

---

## Summary Statistics

Because PrimeKG is categorical and relational, classical statistics (means, standard deviations) don't apply. Instead, distributional and structural summaries describe how knowledge is represented in the graph.

### Node Type Distribution
The graph is heavily skewed toward molecular entities: **gene/protein (~30%)**, **biological process (~23%)**, and **disease (~10%)**. This makes PrimeKG well suited for mechanism-focused analyses, but it also means diseases and drugs are comparatively underrepresented, which can bias analyses toward gene-centric interpretations.

![Node type distribution](images/fig01_node_type_distribution.png)

### Degree-Based Summaries
Highly connected drugs such as **Quinidine (5,200 edges)**, **Chlorpromazine (5,178)**, and **Clozapine (5,150)** dominate the drug-degree distribution. On the disease side, **gallbladder disease (654 drug neighbors)**, **kidney disease (584)**, and **hypertensive disorder (412)** lead. In both cases, high degree reflects research intensity and clinical ubiquity rather than mechanistic specificity.

![Top drugs by degree](images/fig02_top_drugs_by_degree.png)

![Diseases with most drug neighbors](images/fig03_diseases_most_drug_neighbors.png)

> ⚠️ **Implication:** Hub nodes can disproportionately influence similarity scores and repurposing rankings. Results should be treated as hypothesis-generating.

---

## Analyses

Three linked analyses build on one another, mirroring precision-medicine reasoning from understanding a disease to identifying therapeutic gaps.

### 1. Disease Neighborhood Explorer

**Question:** What does a disease's local neighborhood look like in PrimeKG?

For a selected disease, all directly connected nodes are aggregated and summarized by type (drug, gene/protein, phenotype, etc.), showing how the disease is represented across molecular and clinical dimensions.

**Results:** Well-studied diseases have large, gene-dominated neighborhoods. **Breast cancer** has 1,695 neighbors, including 1,667 genes/proteins but only 4 drugs. In contrast, **central nervous system lupus** has just 13 neighbors, highlighting uneven knowledge coverage. This lets users judge upfront whether downstream similarity or repurposing results are likely to be meaningful.

![Neighborhood explorer: breast cancer](images/fig04_neighborhood_breast_cancer.png)

![Neighborhood explorer: CNS lupus](images/fig05_neighborhood_cns_lupus.png)

### 2. Disease–Disease Similarity

**Question:** Which diseases are most similar to a query disease, based on shared genes or shared phenotypes?

Similarity is computed with **Jaccard similarity** over the precomputed neighbor sets:

```
J(A, B) = |A ∩ B| / |A ∪ B|
```

**Results:** The neighbor type is fully tunable, letting the same disease be viewed from two perspectives:

- **Gene/protein similarity** ranks molecularly related cancers highest (e.g., breast neoplasm, hereditary breast carcinoma, hereditary breast–ovarian cancer syndrome), reflecting shared drivers such as BRCA1/2. Useful for exploring shared molecular etiology and therapeutic targets.
- **Effect/phenotype similarity** shifts toward diseases with overlapping clinical manifestations, even if their molecular bases differ. Useful for exploring diagnostic overlap and comorbidity patterns.

![Similarity by gene/protein: breast cancer](images/fig06_similarity_gene_breast_cancer.png)

![Similarity by phenotype: breast cancer](images/fig07_similarity_phenotype_breast_cancer.png)

**Expectations vs. observations:** Results largely matched expectations. However, Jaccard penalizes differences in set size, so a broadly annotated disease can score lower against a more specific hereditary syndrome even when they share key driver genes. This is a known limitation of set-based similarity in unevenly annotated graphs.

### 3. Drug Repurposing Candidates

**Question:** Given a query disease and its most similar diseases, which drugs treat those neighbors but are not yet connected to the query disease?

Candidates are ranked by the **number of supporting similar diseases**, a simple measure of how consistently a drug appears across related disease contexts.

**Results:**
- **Gene/protein-based:** Oncology drugs such as **Docetaxel**, **Paclitaxel**, and **Capecitabine** rank highest, reflecting strong gene- and pathway-level overlap with molecularly similar cancers.
- **Phenotype-based:** **Carboplatin** and **Doxorubicin** rise to the top, supported by fewer but clinically related diseases.
- Drugs appearing under **both** definitions (e.g., Paclitaxel, Carboplatin, Doxorubicin) are more robust signals, supported by both molecular and clinical neighbors.

![Repurposing by gene/protein: breast cancer](images/fig08_repurposing_gene_breast_cancer.png)

![Repurposing by phenotype: breast cancer](images/fig09_repurposing_phenotype_breast_cancer.png)

**Expectations vs. observations:** Some drugs appear broadly across many diseases, reflecting general use rather than specific mechanisms. Diseases with few neighbors yield few or no candidates, showing that repurposing potential is constrained by existing knowledge density.

> ⚠️ **Disclaimer:** These results do not suggest that candidate drugs are novel treatments. This is a hypothesis-generation tool; literature review and experimental or clinical validation are required before drawing any translational conclusions.

---

## Validation

- **Name matching with fuzzy suggestions.** When a disease name doesn't match the graph, the system returns suggested names instead of failing silently.

  ![Fuzzy suggestions for "lupus"](images/fig10_fuzzy_suggestions_lupus.png)

- **Qualitative similarity checks.** Well-known related subtypes consistently ranked above clearly unrelated diseases (e.g., autoimmune cardiomyopathy returns other cardiomyopathies first).

  ![Similarity check: autoimmune cardiomyopathy](images/fig11_similarity_check_autoimmune_cardiomyopathy.png)

- **Non-circular repurposing.** Candidates are restricted to drugs not already connected to the base disease and supported by at least one similar disease. For hypertension, all 412 known drugs in PrimeKG are correctly excluded from the candidate list.

  ![Repurposing check: hypertension](images/fig12_repurposing_check_hypertension.png)

---

## API Reference

The backend is a **Flask** server exposing a JSON-only API. All computation happens in Python; the front end consumes structured JSON.

### `GET /api/disease_profile`
Returns a disease's local neighborhood: total neighbor counts, counts by node type, and connected nodes with their relationship types.

```
GET /api/disease_profile?disease=breast%20cancer
```

```json
{
  "disease_id": "7254",
  "disease_name": "breast cancer",
  "neighbor_counts": {
    "anatomy": 2,
    "disease": 11,
    "drug": 4,
    "effect/phenotype": 11,
    "gene/protein": 1667
  },
  "neighbors": [
    { "id": "DB14655", "name": "Drostanolone propionate", "relation": "indication", "type": "drug" },
    { "id": "DB00112", "name": "Bevacizumab", "relation": "indication", "type": "drug" },
    { "id": "10", "name": "NAT2", "relation": "disease_protein", "type": "gene/protein" }
  ]
}
```

### `GET /api/disease_similarity`
Computes Jaccard similarity against all other diseases.

| Parameter | Description |
|---|---|
| `disease` | Query disease name |
| `neighbor_type` | `gene/protein` or `effect/phenotype` |
| `k` | Number of top matches to return |

```
GET /api/disease_similarity?disease=breast%20cancer&neighbor_type=gene/protein&k=10
```

```json
{
  "neighbor_type": "gene/protein",
  "query_disease_id": "7254",
  "query_disease_name": "breast cancer",
  "results": [
    {
      "disease_id": "21100",
      "disease_name": "breast neoplasm",
      "jaccard": 0.9925373134328358,
      "n_other": 1069,
      "n_query": 1067,
      "n_shared": 1064
    },
    {
      "disease_id": "16419",
      "disease_name": "hereditary breast carcinoma",
      "jaccard": 0.9860982391102873,
      "n_other": 1076,
      "n_query": 1067,
      "n_shared": 1064
    }
  ]
}
```

### `GET /api/repurposing`
Returns drugs linked to similar diseases but not directly connected to the base disease. Called after the similarity endpoint, using the disease IDs it returns.

```
GET /api/repurposing?disease=breast%20cancer&similar_ids=123,456,789
```

```json
{
  "base_disease_id": "7254",
  "base_disease_name": "breast cancer",
  "candidates": [
    {
      "description": "",
      "drug_id": "DB00437",
      "drug_name": "Allopurinol",
      "n_support": 2,
      "supporting_diseases": ["789", "123"]
    },
    {
      "description": "",
      "drug_id": "DB00188",
      "drug_name": "Bortezomib",
      "n_support": 2,
      "supporting_diseases": ["789", "123"]
    }
  ]
}
```

---

## Web Interface

The front end uses **HTML, CSS, and JavaScript**, with **Plotly.js** for interactive charts. It calls the Flask API with `fetch()` and updates visualizations without page reloads.

| Route | Page |
|---|---|
| `/` | Dataset overview and summary statistics |
| `/neighborhood` | Disease Neighborhood Explorer |
| `/similarity` | Disease–Disease Similarity |
| `/repurposing` | Drug Repurposing Candidates |

![Neighborhood page](images/fig17_neighborhood_page.png)

![Similarity page](images/fig18_similarity_page.png)

![Repurposing page](images/fig19_repurposing_page.png)

### Features
- **Analysis selection** via navigation links and per-page controls.
- **Parameter controls** for disease name, neighbor type (gene/protein vs. phenotype), and Top-K.
- **Robust input handling:** case-insensitive matching, fuzzy suggestions for unmatched names, and safe handling of invalid parameters.
- **Interactive visualizations:** Plotly bar charts with hover tooltips and labeled axes.
- **Expandable result lists:** click a bar to view all connected nodes of that type.

![Expandable results list](images/fig20_results_list_hypertension.png)

---

## Getting Started

```bash
# Clone the repository
git clone https://github.com/Bemnet-Kinfe/compmethods-bk675.git
cd "compmethods-bk675/Final Project"

# Install dependencies
pip install -r requirements.txt

# Download PrimeKG CSVs (kg.csv, nodes.csv, disease_features.csv, drug_features.csv)
# from https://doi.org/10.7910/DVN/IXA7BM and place them in the data folder

# Run the app
python app.py
```

Then open `http://localhost:5000` in your browser.

---

## Challenges & Feedback Integration

### Development challenges
- **Slow similarity queries.** The first version rescanned the full knowledge graph on every request. Precomputing disease-to-neighbor sets at load time reduced similarity calculations to fast set operations.
- **Disease name mismatches.** PrimeKG labels don't always match common usage, causing failed exact-match lookups. A fuzzy suggestion method was added.

### Peer feedback incorporated
- **Data cleaning clarity:** documented how node uniqueness is enforced by removing duplicate `node_id` values.
- **Chart readability:** neighborhood bar charts are now sorted in descending order for easier comparison across node types.

---

## Citation & License

**Dataset:** Chowdhury, S. N., & Zitnik, M. (2023). *PrimeKG: A Precision Medicine Knowledge Graph.* Harvard Dataverse. https://doi.org/10.7910/DVN/IXA7BM

PrimeKG is released under the [Creative Commons CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
