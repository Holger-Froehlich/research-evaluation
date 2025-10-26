# Research Evaluation

### Repository Overview
This repository contains workflows developed in the context of field
studies on research evaluation with a current focus on scholarly publications to enhance transparency and reproducibility.
The intended audience are method-oriented users in the field of research evaluation who seek to further develop their work with expert review panels.

With additional modules for research evaluation to follow, it now includes **two major analytical modules**:

#### (A) Interrater Reliability Analysis
A reproducible pipeline for evaluating interrater agreement in research evaluation contexts, including data parsing, scoring aggregation, and statistical reliability estimation, based on a field study on the evaluation of edited books.

#### (B) Conference Proceedings Field Study *(New in v0.2.0)*
A comprehensive workflow analyzing formal quality indicators for conference proceedings. It harmonizes disciplinary classifications across Scopus, SCImago, and Google Scholar, computes SJR-based statistics, and evaluates rule variants for discipline-aware assessment.

**Key components:**
- `interrater-reliability.R` - analytical workflow and visualization pipeline.
- `documentation-interrater-reliability.md` - full documentation.
- `whitelist-based-evaluation.R` – analytical workflow and visualization pipeline.
- `documentation-whitelist-based-evaluation.md` – full documentation.

### Repository structure

    research-evaluation/
      docs/
        documentation-interrater-reliability.md
        documentation-whitelist-based-evaluation.md   # NEW DOCUMENTATION
      src/
        interrater-reliability.R
        whitelist-based-evaluation.R   # NEW MODULE
      CITATION.cff
      LICENSE
      README.md

### Usage
#### (1) Interrater Reliability Module
1. Clone this repository or download the R scripts.
2. Adapt the script to the requirements of your dataset. For `interrater-reliability.R`, this includes:
   - Adjusting the column names to match your data.  
   - Defining the evaluation criteria relevant for your study (insert criteria names in the script).  
   - Setting the rating scale used by your reviewers and thresholds accordingly.  
   - Specifying the file path of your data source.  
3. Open R or RStudio.  
4. Run the scripts in `src/` to reproduce the analyses and visualizations.  

#### (2) Conference Proceedings Field Study
Run `Proceedings_FieldStudy_Workflow.R` sequentially:
1. Harmonize disciplines via ASJC, SCImago, and GS categories.
2. Compute SJR statistics per unified discipline.
3. Generate the five core plots:
   - UpSet intersections between rules
   - Heatmap of index–discipline coverage
   - SJR boxplots by field
   - ?L2 impact bar chart
   - Pareto front of fairness vs. inclusion

Detailed explanations for each figure and step are in `Documentation_ProceedingsFieldStudy.md`.

### Dependencies

R = 4.3.0
The `interrater-reliability` script requires the following R packages:
`readxl`, `dplyr`, `tidyr`, `ggplot2`, `irr`, `gridExtra`, `cowplot`, `igraph`.

The `whitelist-based-evaluation` script requires the following R packages:
`readxl`, `dplyr`, `stringr`, `tidyr`, `purrr`, `janitor`, `ggplot2`, `UpSetR`, `grid`, `scales`.

## Citation

Please see [`CITATION.cff`](CITATION.cff) for citation information.

## License
Source code is licensed under MIT, documentation under CC-BY 4.0.  
See [LICENSE](LICENSE) for details.
