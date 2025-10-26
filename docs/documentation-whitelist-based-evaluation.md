# Documentation: whitelist-based-evaluation.R 

### Use Case: Proceedings Field Study Workflow

## 1. Purpose Note

This repository documents a field study on the **evaluation of conference proceedings** in engineering and computer science, conducted by the AGQ Research Evaluation Team. The goal was to establish **transparent, field-aware formal criteria** for distinguishing between high- and lower-quality conference papers. The workflow combines bibliometric and index-based information to identify systematic, discipline-sensitive differences across proceedings.

The study supports the refinement of **evaluation rules** used by expert panels, demonstrating how combined criteria (e.g., index inclusion, rankings, and citation metrics) affect Level-2 classification outcomes. It serves as a reproducible example for research evaluation processes seeking to balance **factual comparability** and **disciplinary fairness**.

---

## 2. Methodological Background

### 2.1 Motivation
Conference quality signals differ across subdisciplines and publication cultures. Uniform thresholds (e.g., fixed $H5$ values) risk biasing disciplines with smaller communities. This workflow introduces a **two-level quality model** (Level-1 = baseline, Level-2 = high-quality) and empirically tests the effects of extended, index-based rules.

### 2.2 Conceptual Structure
The analytical design integrates three quality dimensions:
1. **Index inclusion** (curated sources such as CPCI, Scopus, SCImago)
2. **Rankings** (CORE conference rankings)
3. **Metrics** (Google Scholar $H5$, SCImago Journal Rank – SJR)

The field study evaluates seven rule variants (R0–R7) built on OR-combinations of these dimensions and examines their effect on Level-2 assignments across disciplines.

---

## 3. Analytical Workflow

### 3.1 Data Preparation Pipeline

The R workflow starts with an **Excel table of conference papers**, each record enriched with metadata from multiple bibliometric sources (Crossref, Scopus, SCImago, and Google Scholar). Each row corresponds to a *Conference Paper* and contains fields such as DOI, proceedings title, and index flags (`is_cpci`, `is_scopus`, `is_scimago`, `is_gs_top20`). This integrated dataset serves as the analytical foundation.

#### (1) Mini-App Pipeline (Python)
A set of modular Python mini-apps was used to extract, identify, and enrich the conference paper data prior to statistical analysis:
1. **Proceedings Extractor** – searches annual research reports for all entries classified as conference proceedings.  
2. **DOI Finder** – retrieves or generates persistent identifiers for each publication.  
3. **Crossref Enricher** – queries Crossref API to retrieve normalized metadata, series titles, ISSNs, and related fields.  
4. **Proceedings Index Matcher** – matches identified proceedings with curated indices (CPCI, Scopus, SCImago, GS Top-20) and records SJR scores where applicable.  

**Output:** An enriched master table consolidating institutional data with verified bibliographic and index-level information.

#### (2) LLM-Aided Harmonization of Field Taxonomies
To ensure consistent disciplinary categorization, the heterogeneous label systems from Scopus (ASJC), SCImago (Categories/Areas), and Google Scholar (Subcategories) were harmonized into **unified disciplinary clusters**.  

A Large Language Model (LLM) assisted in developing 14 semantically coherent “unified clusters.” The model synthesized label proximity based on observed terms across indices and the 951 conference papers. The process incorporated conceptual anchors from international classification systems (OECD FOS, SCImago Areas, DFG Fachsystematik) but remained fully transparent and reproducible.

Human experts validated the model-generated taxonomy, ensuring disciplinary clarity and contextual fit for engineering and computer science. The resulting controlled vocabulary served as the basis for the mapping logic implemented in R.

#### (3) Discipline Consolidation in R
This stage implements the harmonization in code and creates reproducible mappings between index-specific labels and the 14 unified clusters.

**Data Inputs:**  
- `Conference_Papers.xlsx` (enriched conference metadata)  
- `ASJC_Codebook.xlsx` (reference of Scopus ASJC codes)  

**Methodology:**  
- Load and normalize data (`janitor::clean_names()`, lowercase normalization, separator resolution).  
- Identify relevant fields: `Scopus_ASJC`, `Scimago_Category`, `Scimago_Area`, `GS_Category`, and `GS_Subcategory`.  
- Extract and aggregate all field labels into a long-format table (`pivot_longer`, `unnest_longer`).  
- Join ASJC codes with textual descriptions for semantic completeness.  
- Apply the deterministic mapping function `label_to_cluster()`, which uses **regular expressions** to assign each label to one of the 14 unified clusters.

```r
if (str_detect(s, "artificial intelligence|machine learning|pattern recognition"))
  return("AI & Machine Learning")
```

- Generate a **Crosswalk Table** (`discipline_cluster_crosswalk.csv`) linking all observed label names to unified clusters.
- Aggregate paper-level discipline data (`papers_unified_discipline.csv`), allowing fractional assignment where multiple disciplines apply.
- Compute SJR statistics per cluster based on SCImago Conference Series data to establish discipline-specific benchmarks (median, IQR, mean, SD).

**Outputs:**  
1. `discipline_cluster_crosswalk.csv` – mapping of index-specific labels to unified clusters.  
2. `papers_unified_discipline.csv` – enriched paper dataset with harmonized cluster assignments.  
3. `SJR_stats_by_cluster.csv` – computed descriptive statistics for each harmonized discipline.

---

### 3.2 Figures and Statistical Visualizations

The following five R-based plots visualize the relationships between evaluation rules, index coverage, and disciplinary representation.

#### Figure 1 – UpSet Plot: Overlaps between Baseline and Index Rules
**Purpose:** Visualize intersections between the baseline rule and additional index-based rules (CPCI, Scopus, SCImago, GS-Top20).

**Description:** Each vertical bar represents an intersection count, showing how many conference papers are shared by specific rule combinations.

**Technical note:** Generated via `UpSetR` with dichotomous variables (`0/1`), harmonized column naming (`base_line`, `is_cpci`, `is_scopus`, `is_scimago`, `is_gs_top20`).

**Interpretation:** Highlights complementarity between sources—e.g., unique CPCI coverage vs. overlap with Scopus and GS. Large disjoint bars reveal distinctive disciplinary capture.

---

#### Figure 2 – Heatmap: Index Coverage by Harmonized Discipline
**Purpose:** Display disciplinary representation across indices.

**Description:** Rows correspond to harmonized disciplines, columns to indices. The color scale (log-transformed counts) represents the number of proceedings per combination.

**Technical note:** `ggplot2` heatmap using `geom_tile()`, with log10 gradient (`#f7fbff` → `#08306b`). Input derived from long-format index–discipline counts.

**Interpretation:** Reveals which fields dominate each index. Balanced coverage indicates generalizability; gaps suggest potential bias in index curation.

---

#### Figure 3 – Boxplots of SJR Distributions by Harmonized Discipline
**Purpose:** Show variation in citation strength across disciplines.

**Description:** Boxplots summarize SJR values per harmonized cluster. Whiskers denote min/max, the median marks the central tendency.

**Technical note:** Computed from `SJR_stats_by_cluster.csv`. `ggplot2` visualization using log10-scaled y-axis to normalize heavy tails.

**Interpretation:** Fields such as *AI & Machine Learning* or *Networks & Communications* show higher medians, while others (e.g., *Civil Engineering*) cluster lower—reflecting genuine citation intensity differences, not evaluation deficits.

---

#### Figure 4 – ΔL₂ Impact: Effect of Rule Variants on High-Quality Classification
**Purpose:** Quantify relative growth in Level-2 share when new OR-rules are introduced.

**Description:** Bars display the relative increase ($\Delta L_2$) compared to the baseline (R0). Labels indicate both percent gain and absolute new cases.

**Technical note:** Each rule (R1–R7) represents a logical OR expansion (e.g., `base_line == 1 | is_cpci == 1`). Relative increase = $(\text{new L2} / \text{baseline L2})$.

**Interpretation:** Identifies which additional sources yield the largest or most discipline-balanced coverage gains without inflating overall ratings.

---

#### Figure 5 – Pareto Front: Trade-off Between ΔL₂ Growth and Disciplinary Variance
**Purpose:** Assess efficiency of rule variants by jointly minimizing variance across disciplines and maximizing Level-2 growth.

**Description:** Each point represents one rule (R0–R7). The Pareto front connects the non-dominated solutions—those achieving best compromise between fairness and coverage.

**Technical note:** Computed from `sd_results` and `results` tables; standard deviation of discipline counts vs. $\Delta L_2\,(\%)$. Baseline (R0) highlighted by diamond marker.

**Interpretation:** Rules lying on the Pareto front indicate **optimal evaluation configurations**—fairer across disciplines yet still yielding measurable inclusion gains.

---

## 4. Transferability & Outlook

This workflow demonstrates a **scalable evaluation framework** for assessing conference proceedings quality. It can be applied in institutional, regional, or disciplinary contexts where:
- bibliometric indicators differ substantially between subfields,
- expert panels require quantitative pre-structuring before peer review,
- multiple curated and non-curated indices coexist.

The approach supports informed, transparent decisions about weighting schemes (e.g., “5-point vs. 1-point” rules) and offers a reproducible blueprint for comparative evaluation exercises.

---

## 5. Reproducibility Environment

### 5.1 Software Environment
- **R version:** `≥ 4.3.x`
- **Platform:** cross-compatible (Windows, macOS, Linux)
- **Required packages:**
  ```r
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(janitor)
  library(ggplot2)
  library(UpSetR)
  library(grid)
  library(scales)
  ```

### 5.2 Input Data
1. `Conference_Papers.xlsx` – harmonized proceedings with index flags and unified clusters  
2. `ASJC_Codebook.xlsx` – Scopus subject codes reference  
3. `SCIMAGO_IndexData.xlsx` – SCImago Conference Series with SJR values  
4. `discipline_cluster_crosswalk.csv` – label mapping from harmonization step

### 5.3 Execution Order
1. Run discipline harmonization pipeline (creates crosswalk and unified clusters).
2. Compute SJR statistics per cluster.
3. Generate all visualizations sequentially (UpSet, Heatmap, Boxplots, ΔL₂, Pareto Front).
4. Export outputs to CSV and image files for documentation.

### 5.4 Citation
If reused, please cite as:
> Fröhlich, Holger L., (2025). *Proceedings Field Study Workflow: Harmonized Discipline Analysis for Conference Evaluation*. GitHub Repository.

---

**End of Documentation_ProceedingsFieldStudy.md**

