# Research Evaluation

This repository contains R scripts developed in the context of field
studies on research evaluation with a current focus on scholarly publications.  
The intended audience are method-oriented users in the field of research evaluation who seek to further develop their work with expert review panels.  
The initial contribution is on interrater reliability, based on a field study on the evaluation of edited books, with additional modules for research evaluation to follow.

## Repository structure

    research-evaluation/
      README.md
      docs/
        documentation.md      # Detailed methodological documentation
      src/
        interrater-reliability.R
        heuristic-decisions.R
      CITATION.cff
      LICENSE

## Requirements

The R scripts require the following R packages:
- readxl  
- dplyr  
- tidyr  
- ggplot2  
- irr  
- gridExtra  
- cowplot  
- igraph  

## Usage
1. Clone this repository or download the R scripts.
2. Adapt the script to the requirements of your dataset. For `interrater-reliability.R`, this includes:
   - Adjusting the column names to match your data.  
   - Defining the evaluation criteria relevant for your study (insert criteria names in the script).  
   - Setting the rating scale used by your reviewers and thresholds accordingly.  
   - Specifying the file path of your data source.  
3. Open R or RStudio.  
4. Run the scripts in `src/` to reproduce the analyses and visualizations.  


## Citation

Please see [`CITATION.cff`](CITATION.cff) for citation information.

## License
Source code is licensed under MIT, documentation under CC-BY 4.0.  
See [LICENSE](LICENSE) for details.




