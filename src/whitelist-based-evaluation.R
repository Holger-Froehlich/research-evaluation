# ============================================================
# whitelist-based-evaluation.R
# ------------------------------------------------------------
# Title: Index-Based Evaluation of Publication Venues
# Context: Research Evaluation – HAW BW (AG Qualität in der Forschung)
# Purpose: Harmonization, Index Mapping, and Impact Analysis 
#          across Whitelist Systems (CPCI, Scopus, Scimago, GS Top 20)
#          including rule variants, disciplinary balance, and trade-off visualization.
# Author: Holger L. Froehlich
# Date: 2025-10-25
# ============================================================



# ============================================================
# 0. Libraries
# ============================================================
# Data I/O and cleaning: readxl, janitor
# Data wrangling: dplyr, tidyr, stringr, purrr
# Visualization: ggplot2, UpSetR, grid
# Statistical computation: scales
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(janitor)
library(ggplot2)
library(UpSetR)
library(grid)
library(scales)


# ============================================================
# STEP 1 – Data Preparation & Discipline Harmonization
# ============================================================

# ============================================================
# 1.1 INPUT DATA
# ============================================================
# Technical Note:
# Defines all input sources (conference paper file, ASJC taxonomy) and global parameters controlling the harmonization process.
# Inputs:
#   - PAPERS_XLS: Excel file containing conference paper metadata
#   - PAPERS_SHEET: Sheet name for active dataset
#   - ASJC_XLS: Scopus ASJC taxonomy (code → discipline description)
# Parameters:
#   - n_min, max_clusters: thresholds for cluster inclusion and limit
#   - must_keep: list of mandatory clusters to retain regardless of frequency
#   - fallback_cluster: default label for unmapped items
# These inputs define the boundary conditions for subsequent normalization and clustering.
# IMPORTANT: Replace all file paths below with your local paths

# Input Excel with conference paper data:
PAPERS_XLS  <- "../your-conference-paper-data.xslx" -> #insert your data
PAPERS_SHEET<- "Conference Paper 2024"
ASJC_XLS    <- "../ASJC1.xlsx" #insert the providing des ASJC Codes and scientific field categories (scopus taxonomy)

#define unified  disciplines (for background see documentation)
n_min        <- 5
max_clusters <- 15
must_keep    <- c(
  "AI & Machine Learning",
  "Computer Vision & Graphics",
  "Data Systems (DB/IR/IS)",
  "Software & Systems (SE/OS/Arch.)",
  "Networks & Communications",
  "Signal Processing",
  "Robotics & Control",
  "Electrical & Electronic Engineering",
  "Mechanical & Manufacturing & Aerospace",
  "Civil & Environmental Engineering",
  "Materials & Nanotechnology",
  "Chemical & Process Engineering",
  "Biomedical Engineering / Health Tech",
  "Applied Math & Modeling (eng.-adjacent)"
)
fallback_cluster <- "Engineering/CS – General"

# ============================================================
# 1.2 HELPER FUNCTIONS
# ============================================================
# Technical Note:
# Defines reusable functions for normalization, string parsing, column identification,
# and label-to-cluster mapping:
#   - normalize(): lowercases and trims whitespace in text fields
#   - split_multi(): splits multi-value strings into lists by delimiters
#   - detect_sheet(): reads the first sheet name of an Excel file
#   - pick_col(): locates matching column names using regex patterns with error handling
#   - label_to_cluster(): rule-based mapping from free-text discipline labels to
#                         predefined unified clusters (engineering, CS, applied fields)
# Together, these utilities ensure robust text handling and reproducible mapping logic.

normalize <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[\\s\\u00A0]+", " ") %>%
    str_trim()
}
split_multi <- function(x) {
  ifelse(is.na(x) | x=="", NA_character_, x) %>%
    str_split(pattern = "\\s*[;,/\\|]+\\s*")
}
detect_sheet <- function(xlsx) readxl::excel_sheets(xlsx)[1]

# Column identification and error messages
pick_col <- function(df, patterns, required = TRUE, default = NULL) {
  cand <- names(df)[Reduce(`|`, lapply(patterns, function(p) str_detect(names(df), regex(p, ignore_case = TRUE))))]
  if (length(cand) >= 1) return(cand[1])
  if (!is.null(default) && default %in% names(df)) return(default)
  if (required) stop(sprintf("Spalte nicht gefunden. Gesucht nach Mustern: %s. Vorhandene Spalten: %s",
                             paste(patterns, collapse=", "),
                             paste(names(df), collapse=", ")), call. = FALSE)
  return(NA_character_)
}

# ------------------------------------------------------------
# Regex-based label mapping to unified discipline clusters
# ------------------------------------------------------------
label_to_cluster <- function(lbl) {
  s <- normalize(lbl)
  if (is.na(s) || s=="") return(NA_character_)
 
   # --- Computer Science ---
  if (str_detect(s, "artificial intelligence|machine learning|computational intelligence|pattern recognition"))
    return("AI & Machine Learning")
  if (str_detect(s, "computer vision|image processing|computer graphics|graphics"))
    return("Computer Vision & Graphics")
  if (str_detect(s, "information systems|database|databases|information retrieval|data mining|information science"))
    return("Data Systems (DB/IR/IS)")
  if (str_detect(s, "software|programming|operating systems|distributed systems|hardware and architecture|computer architecture|embedded software"))
    return("Software & Systems (SE/OS/Arch.)")
  if (str_detect(s, "computer networks|networks|telecommunications|communications|networking"))
    return("Networks & Communications")
  if (str_detect(s, "signal processing|speech|audio"))
    return("Signal Processing")
  if (str_detect(s, "robotics|control|automation|mechatronics|control and optimization|control systems"))
    return("Robotics & Control")
  
  # --- Engineering ---
  if (str_detect(s, "electrical and electronic|electrical|electronics|microelectronics|vlsi|semiconductor|power electronics"))
    return("Electrical & Electronic Engineering")
  if (str_detect(s, "mechanical|aerospace|manufacturing|industrial and manufacturing|additive manufacturing"))
    return("Mechanical & Manufacturing & Aerospace")
  if (str_detect(s, "civil|structural|geotechnical|building and construction|construction|transportation|hydraulic|environmental engineering"))
    return("Civil & Environmental Engineering")
  if (str_detect(s, "materials|nanoscience|nanotechnology|ceramics|composites|metallurgy"))
    return("Materials & Nanotechnology")
  if (str_detect(s, "chemical engineering|process"))
    return("Chemical & Process Engineering")
  if (str_detect(s, "biomedical engineering|bioengineering|medical devices|biomaterials"))
    return("Biomedical Engineering / Health Tech")
 
   # --- Adjacent Applied Fields ---
  if (str_detect(s, "applied mathematics|applied math|modeling and simulation|simulation|acoustics|computational mechanics"))
    return("Applied Math & Modeling (eng.-adjacent)")
 
   # --- Broader categories ---
  if (str_detect(s, "biochemistry|biotechnology|molecular|genetics|neuroscience|biology|life science"))
    return("Life Sciences")
  if (str_detect(s, "medicine|medical|public health|clinical"))
    return("Health & Medical")
  if (str_detect(s, "business|management|economics|finance|operations research"))
    return("Business/Economics/Management")
  if (str_detect(s, "earth|geology|geoscience|energy") & !str_detect(s, "engineering"))
    return("Energy & Earth Sciences")
  
  # --- Fallback ---
  if (str_detect(s, "computer science")) return("Software & Systems (SE/OS/Arch.)")
  if (str_detect(s, "engineering"))      return("Engineering/CS – General")
  return(NA_character_)
}

# ============================================================
# 1.3 Label Normalization & Extraction
# ============================================================
# Technical Note:
# Loads and cleans the conference paper dataset and ASJC taxonomy,
# extracts all discipline labels from multiple sources (ASJC, SCImago, GS),
# and harmonizes them into a unified long format.
# Steps:
#   1. Read and clean the input Excel file, add unique .rowid per record.
#   2. Identify available label columns (ASJC, SCImago, GS categories).
#   3. Load ASJC reference table (code → description, area), clean and deduplicate.
#   4. Expand all multi-valued label fields into individual rows.
#   5. Join ASJC codes with descriptions and merge with non-ASJC labels.
# Output: a harmonized label table per paper for subsequent cluster mapping.


## Load and Clean Data
papers_raw <- readxl::read_excel(PAPERS_XLS, sheet = PAPERS_SHEET)
papers <- papers_raw %>%
  janitor::clean_names() %>%
  tibble::rowid_to_column(var = ".rowid")  # schreibt .rowid garantiert

## Harmonize Discipline Labels Across Indices

### identify relevant columns
getcol <- function(pattern) {
  nm <- names(papers)[str_detect(names(papers), regex(pattern, ignore_case = TRUE))]
  if (length(nm)==0) NA_character_ else nm[1]
}
col_asjc     <- getcol("^scopus_asjc$")
col_sci_cat  <- getcol("^scimago_category$|^scimago_categories$")
col_sci_area <- getcol("^scimago_areas?$")
col_gs_cat   <- getcol("^gs_category$")
col_gs_sub   <- getcol("^gs_subcategory$")

# read in of ASJC-reference
asjc_sheet <- detect_sheet(ASJC_XLS)
asjc_ref0  <- readxl::read_excel(ASJC_XLS, sheet = asjc_sheet) %>% janitor::clean_names()

asjc_code_col <- pick_col(asjc_ref0, c("asjc.*code", "^code$"))
asjc_desc_col <- pick_col(asjc_ref0, c("desc", "category", "name"))
# "Area" is optional
asjc_area_col <- tryCatch(pick_col(asjc_ref0, c("area.*desc","area.*name"), required = FALSE), error = function(e) NA_character_)

asjc_ref <- asjc_ref0 %>%
  transmute(
    asjc_code = as.character(.data[[asjc_code_col]]) %>% str_extract("\\d{3,5}"),
    asjc_desc = as.character(.data[[asjc_desc_col]]),
    asjc_area = if (!is.na(asjc_area_col)) as.character(.data[[asjc_area_col]]) else NA_character_
  ) %>%
  filter(!is.na(asjc_code), asjc_code != "") %>%
  distinct(asjc_code, .keep_all = TRUE)


## Long-table of all existing labels per Paper
labels_long <- papers %>%
  transmute(
    .rowid,
    asjc_codes        = if (!is.na(col_asjc)) split_multi(.data[[col_asjc]]) else list(NA_character_),
    scimago_category  = if (!is.na(col_sci_cat)) split_multi(.data[[col_sci_cat]]) else list(NA_character_),
    scimago_area      = if (!is.na(col_sci_area)) split_multi(.data[[col_sci_area]]) else list(NA_character_),
    gs_category       = if (!is.na(col_gs_cat)) split_multi(.data[[col_gs_cat]]) else list(NA_character_),
    gs_subcategory    = if (!is.na(col_gs_sub)) split_multi(.data[[col_gs_sub]]) else list(NA_character_)
  ) %>%
  pivot_longer(cols = -c(.rowid), names_to = "source", values_to = "vals") %>%
  unnest_longer(vals, keep_empty = TRUE) %>%
  mutate(vals = vals %>% as.character() %>% str_trim()) %>%
  filter(!is.na(vals), vals != "") %>%
  distinct(.rowid, source, vals)

## ASJC-Codes → join descriptions
labels_asjc <- labels_long %>%
  filter(source == "asjc_codes") %>%
  mutate(asjc_code = str_extract(vals, "\\d{3,5}")) %>%
  left_join(asjc_ref, by = c("asjc_code" = "asjc_code")) %>%
  transmute(.rowid,
            source = "ASJC",
            original_label = coalesce(asjc_desc, vals),
            area_label     = asjc_area)

## non-ASJC Labels (SCImago/GS)
labels_other <- labels_long %>%
  filter(source != "asjc_codes") %>%
  transmute(.rowid,
            source = toupper(source),
            original_label = vals,
            area_label = NA_character_)

# ============================================================
# 1.4 UNIFIED CLUSTER ASSIGNMENT
# ============================================================
# Technical Note:
# Integrates all label sources into a unified structure and applies rule-based mapping.
# Steps:
#   - Normalize each original label for matching.
#   - Apply label_to_cluster() to both original and area labels.
#   - Merge rule-based results into a preliminary cluster assignment.
# Result: each paper is linked to one or several preliminary unified clusters.


labels_all <- bind_rows(labels_asjc, labels_other) %>%
  mutate(
    original_label_norm = normalize(original_label),
    cluster_rule        = map_chr(original_label, label_to_cluster),
    area_rule           = ifelse(!is.na(area_label), label_to_cluster(area_label), NA_character_),
    prelim_cluster      = coalesce(cluster_rule, area_rule)
  ) %>%
  distinct(.rowid, source, original_label, prelim_cluster)

# ============================================================
# 1.5 Cluster Frequency AND SELECTION
# ============================================================
# Technical Note:
# Quantifies how frequently each preliminary cluster appears and defines the final
# set of clusters to keep under the configured thresholds.
# Steps:
#   1. Count papers per cluster (distinct .rowid).
#   2. Retain clusters meeting minimum frequency (n_min) or listed in must_keep.
#   3. If the total exceeds max_clusters, trim least-frequent clusters dynamically.
#   4. Assign fallback_cluster to residual or unmatched entries.
# Produces the finalized list of clusters (labels_final) for all papers.


cluster_freq <- labels_all %>%
  filter(!is.na(prelim_cluster)) %>%
  distinct(.rowid, prelim_cluster) %>%
  count(prelim_cluster, name = "n_papers") %>%
  arrange(desc(n_papers))

select_clusters <- function(tbl, n_cut, must_keep, max_clusters) {
  keep <- tbl %>%
    mutate(keep = n_papers >= n_cut | prelim_cluster %in% must_keep) %>%
    filter(keep) %>%
    pull(prelim_cluster) %>%
    unique()
  if (length(keep) > max_clusters) {
    base <- tbl %>% arrange(desc(n_papers))
    must <- intersect(must_keep, base$prelim_cluster)
    rest <- setdiff(base$prelim_cluster, must)
    need <- max(0, max_clusters - length(must))
    keep <- c(must, head(rest, need))
  }
  keep
}

n_cur <- n_min
repeat {
  kept <- select_clusters(cluster_freq, n_cur, must_keep, max_clusters)
  if (length(kept) <= max_clusters) break
  n_cur <- n_cur + 1
  if (n_cur > 50) break
}

labels_final <- labels_all %>%
  mutate(
    unified_cluster = case_when(
      prelim_cluster %in% kept ~ prelim_cluster,
      is.na(prelim_cluster)    ~ NA_character_,
      TRUE                     ~ fallback_cluster
    )
  )

# ============================================================
# 1.6 CROSSWALK AND PAPER AGGREGATION
# ============================================================
# Technical Note:
# Builds final output structures linking original labels and papers to unified clusters.
# Steps:
#   1. Create crosswalk table (original_label → unified_cluster) for transparency.
#   2. Aggregate clusters per paper with fractional weighting (1/n_clusters).
#   3. Merge aggregated cluster data back to original paper IDs (.rowid).
#   4. Export both crosswalk and paper-cluster tables as CSV files.
#   5. Generate a console summary reporting coverage and top clusters.
# Outputs:
#   - discipline_cluster_crosswalk.csv
#   - papers_unified_discipline.csv
#   - diagnostic messages for validation.

## Crosswalk Table
discipline_cluster_crosswalk <- labels_final %>%
  filter(!is.na(unified_cluster)) %>%
  distinct(source, original_label, unified_cluster) %>%
  arrange(unified_cluster, source, original_label)

## Paper-Cluster (fractional)
paper_clusters <- labels_final %>%
  filter(!is.na(unified_cluster)) %>%
  distinct(.rowid, unified_cluster) %>%
  group_by(.rowid) %>%
  mutate(n_clusters = n()) %>%
  ungroup() %>%
  mutate(weight = 1/n_clusters)

paper_clusters_aggregated <- paper_clusters %>%
  group_by(.rowid) %>%
  summarise(
    unified_clusters = paste0(unified_cluster, collapse = "; "),
    weights          = paste0(format(round(weight, 4), nsmall = 4), collapse = "; "),
    n_clusters       = n(),
    .groups = "drop"
  )

paper_clusters_aligned <- papers %>%
  select(.rowid) %>%
  left_join(paper_clusters_aggregated, by = ".rowid") %>%
  arrange(.rowid)

## save results
write.csv(discipline_cluster_crosswalk, "discipline_cluster_crosswalk.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(paper_clusters_aligned,       "papers_unified_discipline.csv",   row.names = FALSE, fileEncoding = "UTF-8") # for simplicity the results were merged to the Conference paper Excel sheet and used as input data below 

## Mini-Report
message("# Clustersteuerung")
message(sprintf("- verwendetes n_cut: %s", n_cur))
message(sprintf("- finale Clusterzahl (exkl. NA): %s", length(unique(na.omit(labels_final$unified_cluster)))))

message("# Abdeckung")
covered <- sum(!is.na(paper_clusters_aligned$unified_clusters))
message(sprintf("- Papers mit mind. 1 Cluster: %s / %s", covered, nrow(paper_clusters_aligned)))

message("# Top-Cluster (fractional)")
print(
  paper_clusters %>%
    group_by(unified_cluster) %>%
    summarise(frac_count = sum(weight), .groups="drop") %>%
    arrange(desc(frac_count)) %>%
    head(15)
)

# resulting objects:
# discipline_cluster_crosswalk  (original label -> cluster)
# paper_clusters_aligned        (NA for papers without discipline, due to lack of data in Crossref)


# ============================================================
# STEP 2 – Compute SJR Statistics per Harmonized Discipline
# ============================================================
# Input:
#   - 251016_listen_konsolidiert.xlsx  (Sheet "SCImago_2024")
#   - discipline_cluster_crosswalk.csv  (mapping labels → clusters)
# Output:
#   - SJR_stats_by_cluster.csv          (Median, IQR, mean, etc.)
# ==========================================


# ============================================================
# 2.1 – LOAD SCIMAGO DATA AND CROSSWALK SHEET
# ============================================================
# Technical Note:
# Imports the SCImago index (journal and conference series metadata)
# and the discipline crosswalk generated in Step 1.
# - SCIMAGO_XLS provides SJR scores and subject categories per source.
# - CROSSWALK_CSV maps raw category labels to unified discipline clusters.
# Both datasets are cleaned and normalized to lowercase strings,
# producing a consistent key ('original_label_norm') for joining.

## Paths
SCIMAGO_XLS   <- "..//251016_listen_konsolidiert.xlsx" # load your own index list 
SCIMAGO_SHEET <- "SCImago_2024" #-> insert your own excel sheet name here
CROSSWALK_CSV <- "../discipline_cluster_crosswalk.csv"

## Load data
scimago <- readxl::read_excel(SCIMAGO_XLS, sheet = SCIMAGO_SHEET) %>%
  clean_names()

crosswalk <- read.csv(CROSSWALK_CSV, encoding = "UTF-8", check.names = FALSE) %>%
  clean_names() %>%
  transmute(
    original_label_norm = str_squish(tolower(original_label)),
    unified_cluster
  ) %>%
  distinct()


# =======================================================================================
# 2.2 – FILTER AND JOIN SCIMAGO CONF.SERIES WITH UNIFIED DISCIPLINES FROM CROSSWALK SHEET
# =======================================================================================
# Technical Note:
# Prepares and filters SCImago category data, aligns it with unified discipline clusters.
# Steps:
#  1. Extracts core columns (sourceid, title, sjr, categories, areas).
#  2. Converts SJR values to numeric and filters valid (> 0) entries.
#  3. Merges 'categories' and 'areas' into a unified label field and expands multi-labels.
#  4. Normalizes text (lowercase, trimmed) for consistent matching.
#  5. Retains only categories present in the crosswalk and joins them to obtain unified clusters.
# Result: a harmonized dataset (scimago_filtered) linking each conference series
#         with its SJR score and unified discipline cluster.

## Preprocess SCImago data
scimago_long <- scimago %>%
  select(sourceid, title, sjr, categories, areas) %>%
  mutate(
    sjr = suppressWarnings(as.numeric(sjr))
  ) %>%
  filter(!is.na(sjr), sjr > 0) %>%
  unite("cat_area", categories, areas, sep = "; ", na.rm = TRUE) %>%
  separate_rows(cat_area, sep = "\\s*[;,/\\|]+\\s*") %>%
  mutate(
    cat_area = str_trim(cat_area),
    category_norm = str_squish(tolower(cat_area))
  )

## Keep only categories in crosswalk
allowed_labels <- unique(crosswalk$original_label_norm)

scimago_filtered <- scimago_long %>%
  filter(category_norm %in% allowed_labels) %>%
  inner_join(crosswalk, by = c("category_norm" = "original_label_norm")) %>%
  distinct(sourceid, unified_cluster, sjr)


# ========================================
# 2.3 – Compute SJR statistics per cluster
# ========================================
# Technical Note:
# Aggregates SJR metrics at the level of unified clusters.
# For each cluster, calculates:
#   - n_series: number of matched conference series
#   - sjr_median, sjr_mean, sjr_iqr, sjr_min, sjr_max, sjr_sd
# Produces a ranked summary (by median SJR) representing field-normalized citation quality.
# The resulting table 'sjr_stats_by_cluster' is written to CSV
# and previewed in console for validation.

sjr_stats_by_cluster <- scimago_filtered %>%
  group_by(unified_cluster) %>%
  summarise(
    n_series   = n(),
    sjr_median = median(sjr, na.rm = TRUE),
    sjr_iqr    = IQR(sjr, na.rm = TRUE),
    sjr_mean   = mean(sjr, na.rm = TRUE),
    sjr_min    = min(sjr, na.rm = TRUE),
    sjr_max    = max(sjr, na.rm = TRUE),
    sjr_sd     = sd(sjr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sjr_median))

## Save output
write.csv(sjr_stats_by_cluster,
          "../SJR_stats_by_cluster.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

## Quick summary to console
# ---------------------------
message("JR statistics successfully computed")
print(sjr_stats_by_cluster, n = 20)

# ============================================================
# STEP 3 – Exploratory Visualization & Evaluation
# ============================================================

# ======================
# 3.1 – PLOT 1 Upsetplot
# ======================
# Technical Note:
# Visualizes the intersection of index-based classifications (Baseline, CPCI, Scopus, Scimago, GS Top 20).
# Steps:
#   1. Convert logical or numeric indicators into binary (0/1) membership across five index columns.
#   2. Use UpSetR to plot the size and overlap of these sets.
# Purpose:
#   - Identify how many conference papers are covered by multiple index systems.
#   - Highlight unique vs. shared coverage patterns among indices.
# Output:
#   Interactive-style bar and matrix display of set intersections.


## Load Data
PAPERS_XLS   <- "../251023_Konferenzpaper_JB24_with_doi_with_crossref_with_proceedings_harmonized_disciplines_final.xlsx" # load your data sheet containing the Conference papers with assigend unified disciplines
PAPERS_SHEET <- "Konferenzpaper 2024" # define your excel data sheet

df <- read_excel(PAPERS_XLS, sheet = PAPERS_SHEET) |> janitor::clean_names()

set_cols <- c("base_line","is_cpci","is_scopus","is_scimago","is_gs_top20")

df_plot <- df %>%
  select(all_of(set_cols)) %>%
  mutate(across(everything(), ~ ifelse(.x %in% c(1, TRUE), 1L, 0L))) %>%
  mutate(across(everything(), ~ replace_na(.x, 0L))) %>%
  mutate(across(everything(), as.integer)) %>%
  as.data.frame()

# UpSet-Plot
AGQ_color <- "#006699"
upset(df_plot,
      sets = set_cols,
      order.by = "freq",
      keep.order = TRUE,
      sets.bar.color = AGQ_color,
      main.bar.color = AGQ_color,
      matrix.color = AGQ_color,
      mainbar.y.label = "Intersection Size",
      sets.x.label   = "Set Size",
      point.size = 3.0,   
      line.size = 0.6,
      text.scale = c(1.6, 1.2, 1.1, 1, 1.3, 1.3))

# ======================================================
# 3.2 – PLOT 2 Heatmaps indices x unified disciplines
# ======================================================
# Technical Note:
# Depicts how each index (CPCI, Scopus, Scimago, GS Top 20, Baseline) distributes across harmonized disciplines.
# Steps:
#   1. Expand multi-assigned disciplines into individual rows.
#   2. Pivot index flags into long format (index × discipline).
#   3. Count occurrences and order disciplines by total frequency.
#   4. Plot as logarithmic heatmap (index → x-axis, discipline → y-axis).
# Purpose:
#   - Reveal which indices emphasize which fields.
#   - Assess coverage density and disciplinary bias.
# Output:
#   Heatmap with color intensity scaled to log(count), styled in AGQ corporate blue.

## define theme and Colour
AGQ_color <- "#006699"

theme_AGQ <- function(base_size = 15, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.9),
      axis.text.y = element_text(size = base_size * 0.85),
      axis.title  = element_text(face = "bold", size = base_size * 1.1),
      plot.title  = element_text(face = "bold", size = base_size * 1.2),
      plot.subtitle = element_text(size = base_size * 0.95, color = "grey30"),
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      legend.text  = element_text(size = base_size * 0.85),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

## load data
PAPERS_XLS   <- "../251023_Konferenzpaper_JB24_with_doi_with_crossref_with_proceedings_harmonized_disciplines_final.xlsx" # load your own paper data
PAPERS_SHEET <- "Konferenzpaper 2024"

df <- read_excel(PAPERS_XLS, sheet = PAPERS_SHEET) |> janitor::clean_names()
set_cols <- c("base_line","is_cpci","is_scopus","is_scimago","is_gs_top20")


## prepare column "unified disciplines" 
df_disc <- df %>%
  filter(!is.na(unified_clusters)) %>%
  separate_rows(unified_clusters, sep = ";\\s*") %>%
  mutate(unified_clusters = str_trim(unified_clusters)) %>%
  filter(unified_clusters != "")

## bring Indices into long format
df_long <- df_disc %>%
  select(all_of(set_cols), unified_clusters) %>%
  pivot_longer(all_of(set_cols), names_to = "index", values_to = "in_set") %>%
  filter(in_set %in% c(1, TRUE))

## Counts per index × discipline
heat_data <- df_long %>%
  group_by(index, unified_clusters) %>%
  summarise(count = n(), .groups = "drop")

## sort disciplines according to total frequency
disc_order <- heat_data %>%
  group_by(unified_clusters) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(unified_clusters)

# Heatmap - reduced version
#--------------------------

## exclude NA and disciplines without counts 
heat_data <- heat_data %>%
  filter(!is.na(unified_clusters), unified_clusters != "NA", unified_clusters != "")

## plot
p_heatmap <- ggplot(
  heat_data,
  aes(
    x = index,
    y = factor(unified_clusters, levels = disc_order),
    fill = count
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#08306b"),
    trans = "log10",
    breaks = c(1, 5, 10, 50, 100, 500),
    name = "Count"
  ) +
  labs(
    title = "Disziplinäre Abdeckung der Indizes"
  ) +
  theme_AGQ() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

p_heatmap


# =============================================================
# 3.3 PLOT 3 - boxplots SJR by unified disiplines in scimago (n~650 conf. series x unified disciplines)
# =============================================================
# Technical Note:
# Summarizes citation impact (SJR) distributions across unified disciplines.
# Steps:
#   1. Use aggregated SJR statistics per cluster (median, IQR, min–max).
#   2. Render boxplots using precomputed summary stats with log-scaled y-axis.
#   3. Annotate each box with the number of indexed series (n).
# Purpose:
#   - Compare typical citation strength of conference series by discipline.
#   - Identify fields with high or heterogeneous SJR distributions.
# Output:
#   Horizontal log-scale boxplot ranked by median SJR.

## define theme
AGQ_color <- "#006699"

theme_AGQ <- function(base_size = 15, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x = element_text(size = base_size * 0.9),
      axis.text.y = element_text(size = base_size * 0.85),
      axis.title  = element_text(face = "bold", size = base_size * 1.1),
      plot.title  = element_text(face = "bold", size = base_size * 1.2),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank()
    )
}

## Plot (logarithmic scale on y-axis + counts)
p_sjr_box <- sjr_stats_by_cluster %>%
  mutate(unified_cluster = reorder(unified_cluster, sjr_median)) %>%
  ggplot(aes(
    x = unified_cluster,
    ymin = sjr_min,
    lower = sjr_median - sjr_iqr / 2,
    middle = sjr_median,
    upper = sjr_median + sjr_iqr / 2,
    ymax = sjr_max
  )) +
  geom_boxplot(
    stat = "identity",
    fill = AGQ_color,
    color = "black",
    width = 0.6,
    alpha = 0.7
  ) +
  # log10-Skalierung für SJR-Werte
  scale_y_continuous(trans = "log10", breaks = c(0.1, 0.5, 1, 2, 5, 10, 20)) +
  # Counts als Label
  geom_text(
    aes(
      y = sjr_max * 1.1,   # etwas rechts von der Box
      label = paste0("(n=", n_series, ")")
    ),
    hjust = 0,
    size = 4.2
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "SJR-Verteilungen nach Disziplin",
    x = NULL,
    y = "SJR (log-Skala)"
  ) +
  theme_AGQ() +
  theme(
    plot.margin = margin(10, 60, 10, 10)  # Platz für Labels
  )

p_sjr_box

# ============================
# 3.4 PLOT 4 - ΔL2 Impact Plot
# ============================

# Technical Note:
# Measures the relative growth of high-quality (L₂) classifications under alternative rule sets.
# Steps:
#   1. Define rule variants combining Baseline with CPCI, Scopus, Scimago, and GS Top 20 indicators.
#   2. Evaluate each rule against the Baseline to count newly classified L₂ cases.
#   3. Compute relative increase (ΔL₂ = new / baseline).
#   4. Plot as bar chart with percentage labels and embedded legend of rule definitions.
# Purpose:
#   - Quantify the effect of inclusion rules on overall classification coverage.
#   - Compare alternative heuristics for broadened recognition criteria.
# Output:
#   Bar chart showing ΔL₂ % per rule variant plus textual legend annotation.

## define theme
AGQ_color <- "#006699"

theme_AGQ <- function(base_size = 15, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x = element_text(size = base_size * 0.9),
      axis.text.y = element_text(size = base_size * 0.85),
      axis.title  = element_text(face = "bold", size = base_size * 1.1),
      plot.title  = element_text(face = "bold", size = base_size * 1.2),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

## load data
PAPERS_XLS   <- "../251023_Konferenzpaper_JB24_with_doi_with_crossref_with_proceedings_harmonized_disciplines_final.xlsx" # load your own data
PAPERS_SHEET <- "Konferenzpaper 2024" # define your sheet containing the confernce papers

df <- read_excel(PAPERS_XLS, sheet = PAPERS_SHEET) |> janitor::clean_names()

# ensure, that columns are  0/1
df <- df %>%
  mutate(across(c(base_line, is_cpci, is_scopus, is_scimago, is_gs_top20),
                ~ as.integer(.x %in% c(1, "1", TRUE, "TRUE"))))


## define rule variants
rules <- tribble(
  ~Regel, ~Label, ~Formel,
  "R1", "R0 ODER CPCI", "base_line == 1 | is_cpci == 1",
  "R2", "R0 ODER Scopus", "base_line == 1 | is_scopus == 1",
  "R3", "R0 ODER GS Top 20", "base_line == 1 | is_gs_top20 == 1",
  "R5", "R0 ODER Scimago ODER GS Top 20", "base_line == 1 | is_scimago == 1 | is_gs_top20 == 1",
  "R6", "R0 ODER CPCI ODER Scopus", "base_line == 1 | is_cpci == 1 | is_scopus == 1",
  "R7", "R0 ODER CPCI ODER Scopus ODER GS Top 20", "base_line == 1 | is_cpci == 1 | is_scopus == 1 | is_gs_top20 == 1"
)

## calculation per rule variant
results <- map_dfr(1:nrow(rules), function(i) {
  
  rule_expr <- rules$Formel[i]
  reg_name  <- rules$Regel[i]
  label     <- rules$Label[i]
  
  in_rule <- with(df, eval(parse(text = rule_expr)))
  
  # Baseline-Cases
  base_L2 <- sum(df$base_line == 1, na.rm = TRUE)
  # new cases (now 1, before 0)
  new_L2  <- sum((in_rule == 1) & (df$base_line == 0), na.rm = TRUE)
  
  tibble(
    Regel = reg_name,
    Label = label,
    new_L2 = new_L2,
    base_L2 = base_L2,
    rel_increase = new_L2 / base_L2
  )
})

## Plot
p_delta <- ggplot(results, aes(x = Regel, y = rel_increase)) +
  geom_col(width = 0.7, fill = AGQ_color, color = "black", alpha = 0.85) +
  geom_text(
    aes(label = paste0("+", scales::percent(rel_increase, accuracy = 0.1),
                       "\n(n=", new_L2, ")")),
    vjust = -0.6, size = 4.2, color = "black"
  ) +
  labs(
    title = expression(paste("Impact Regelvarianten [", Delta, " L2 '5fach werten', vs R0]")),
    y = expression(Delta~L[2]~"Zuwachs relativ zu R0"),
    x = NULL,
    caption = expression(Delta~L[2] == frac("Neue Fälle in Regel R[i] (BaseLine==0)",
                                            "Baseline-Fälle (BaseLine==1)")~";"~Delta~L[2] >= 0)
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.2))) +
  theme_AGQ()

## legend box
legend_text <- paste(
  "Regeldefinitionen:",
  "R0 – Baseline (GS H5 ≥ 30 ODER CORE A–B)",
  "R1 – R0 ODER CPCI",
  "R2 – R0 ODER Scopus (Serial with Profile)",
  "R3 – R0 ODER GS Top 20",
  "R5 – R0 ODER Scimago ODER GS Top 20",
  "R6 – R0 ODER CPCI ODER Scopus",
  "R7 – R0 ODER CPCI ODER Scopus ODER GS Top 20",
  sep = "\n"
)

p_final <- p_delta +
  annotation_custom(
    grob = textGrob(
      legend_text,
      x = 0.2,    # weiter links → über R1
      y = 0.98,    # etwas höher im Plotbereich
      just = c("left", "top"),
      gp = gpar(fontsize = 10, fontfamily = "sans")
    )
  )

p_final

# ==============================================================================
# 3.5 PLOT 5 - Pareto Plot - Trade-off Between ΔL₂ Growth and Disciplinary Variance 
# ==============================================================================
# Technical Note:
# Analyzes trade-off between classification growth (ΔL₂ %) and loss of disciplinary balance.
# Steps:
#   1. Compute per-discipline counts and standard deviation (SD_Disziplin) for each rule variant.
#   2. Combine with ΔL₂ growth data from previous plot.
#   3. Identify Pareto-optimal rules minimizing variance for given gain.
#   4. Plot SD (x-axis) vs ΔL₂ % (y-axis) with Pareto-front and Baseline marker.
# Purpose:
#   - Evaluate efficiency of rule variants: maximal L₂ gain with minimal disciplinary distortion.
#   - Visualize balance point between quantitative expansion and field equity.
# Output:
#   Pareto scatterplot with front line, annotated rule points, and summary table.

## prepare discipline distribution
df_disc <- df %>%
  filter(!is.na(unified_clusters), unified_clusters != "NA", unified_clusters != "") %>%
  separate_rows(unified_clusters, sep = ";\\s*") %>%
  mutate(unified_clusters = str_trim(unified_clusters))

disc_list <- sort(unique(df_disc$unified_clusters))

## define rule variants
rules_var <- tribble(
  ~Regel, ~Label, ~Formel,
  "R0", "Baseline (Base_Line == 1)", "base_line == 1",
  "R1", "R0 ODER CPCI", "base_line == 1 | is_cpci == 1",
  "R2", "R0 ODER Scopus", "base_line == 1 | is_scopus == 1",
  "R3", "R0 ODER GS Top 20", "base_line == 1 | is_gs_top20 == 1",
  "R5", "R0 ODER Scimago ODER GS Top 20", "base_line == 1 | is_scimago == 1 | is_gs_top20 == 1",
  "R6", "R0 ODER CPCI ODER Scopus", "base_line == 1 | is_cpci == 1 | is_scopus == 1",
  "R7", "R0 ODER CPCI ODER Scopus ODER GS Top 20", "base_line == 1 | is_cpci == 1 | is_scopus == 1 | is_gs_top20 == 1"
)

## Standardabweichung der discipline distribution per rule variant
sd_results <- map_dfr(1:nrow(rules_var), function(i) {
  rule_expr <- rules_var$Formel[i]
  reg_name  <- rules_var$Regel[i]
  label     <- rules_var$Label[i]
  
  in_rule <- with(df_disc, eval(parse(text = rule_expr)))
  df_rule <- df_disc[in_rule, ]
  
  counts <- df_rule %>%
    count(unified_clusters) %>%
    complete(unified_clusters = disc_list, fill = list(n = 0)) %>%
    arrange(unified_clusters)
  
  tibble(
    Regel = reg_name,
    Label = label,
    SD_Disziplin = sd(counts$n),
    Total = sum(counts$n)
  )
})

## Combine with impact data (Δ L₂)
pareto_df <- sd_results %>%
  select(Regel, SD_Disziplin) %>%
  left_join(results %>% select(Regel, rel_increase), by = "Regel") %>%
  mutate(Delta_L2_pct = rel_increase * 100)

## Pareto-Front
pareto_df <- pareto_df %>%
  arrange(SD_Disziplin, Delta_L2_pct)

pareto_front <- pareto_df %>%
  filter(!duplicated(cummin(Delta_L2_pct))) %>%
  arrange(SD_Disziplin)

## AGQ-Theme
AGQ_color <- "#006699"

theme_AGQ <- function(base_size = 15, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x = element_text(size = base_size * 0.9),
      axis.text.y = element_text(size = base_size * 0.85),
      axis.title  = element_text(face = "plain", size = base_size * 1.1),
      plot.title  = element_text(face = "plain", size = base_size * 1.2),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

baseline_sd <- pareto_df %>%
  filter(Regel == "R0")

## Plot
p_pareto <- ggplot(pareto_df, aes(x = SD_Disziplin, y = Delta_L2_pct, label = Regel)) +
  # alle Regelpunkte
  geom_point(size = 4, color = AGQ_color, alpha = 0.8) +
  geom_text(vjust = -0.7, size = 4) +
  # Pareto-Front
  geom_line(data = pareto_front, aes(x = SD_Disziplin, y = Delta_L2_pct),
            color = AGQ_color, linewidth = 1.2) +
  geom_point(data = pareto_front, color = AGQ_color, size = 5, shape = 21, fill = "white") +
  # Baseline-Marker
  geom_point(data = baseline_sd,
             aes(x = SD_Disziplin, y = Delta_L2_pct),
             shape = 23, size = 6, fill = "white", color = "black", stroke = 1.3) +
  geom_text(data = baseline_sd,
            aes(label = "R0 (Baseline)"),
            vjust = -1, fontface = "bold", size = 4.5) +
  # Referenzlinien (AGQ-Blau)
  geom_vline(xintercept = sd_results$SD_Disziplin [1], color = AGQ_color, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = AGQ_color, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = expression(paste("Zuwachs L2-Bewertung (5-fach) vs.  Disziplinäre Streuung")),
    x = "Standardabweichung der Disziplinverteilung (n = 15 Disziplinen)",
    y = expression(Delta~L[2]~"[Zuwachs 5-fach-Wertungen in %]")
  ) +
  theme_AGQ()

## pareto-optimal rules
pareto_table <- pareto_front %>%
  arrange(SD_Disziplin, Delta_L2_pct) %>%
  select(Regel, SD_Disziplin, Delta_L2_pct)

print(pareto_table)
p_pareto
