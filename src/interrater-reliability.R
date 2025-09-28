# ------------------------------------------------------------
# Interrater Reliability Analysis for Edited Volumes
# ------------------------------------------------------------
# This script covers:
# (A) Balloon plots with ICC and weighted Kappa per criterion
# (B) Heatmaps of pairwise deviations (mean deviation and frequency of strong deviations)
# (C) Network view (Louvain communities) for pairs with strong deviations
#
# Notes:
# - Replace data paths in the two file_path variables with your local files
# - Column names are expected as in the shared data (German field names)
# ------------------------------------------------------------

# Load packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(irr)
library(gridExtra)
library(cowplot)  # for get_legend()
library(igraph)

# ============================================================
# (A) Balloon plots with ICC and Kappa
# ============================================================

# --- Import data (set generic path & sheet) ---
file_path_balloon <- "data/inter_rater.xlsx"
sheet_name        <- "Inter_Rater" # the excel sheet is organized with the following collumns: Item ID; rater.1; criterion.1.1; ...; criterion.n.1; rater.2; criterion.1.2;...; criterion.n.2

data <- read_excel(file_path_balloon, sheet = sheet_name)

# --- Prepare criteria and containers ---
kriterien <- c("Kanal", "Kohaerenz", "Herausgabe", "Auswahlverfahren", "Gesamtrating") # names of the Criteria in our use case
bubble_data_list <- list()
iccs   <- list()
kappas <- list()

# --- Build balloon-plot tables + compute ICC/Kappa per criterion ---
for (kriterium in kriterien) {
  var1 <- data[[paste0(kriterium, "_1")]]
  var2 <- data[[paste0(kriterium, "_2")]]
  
  df <- data.frame(x = var1, y = var2)
  df_table <- as.data.frame(table(df))
  names(df_table) <- c("x", "y", "Freq")
  df_table$Kriterium <- kriterium
  df_table$x <- as.numeric(as.character(df_table$x))
  df_table$y <- as.numeric(as.character(df_table$y))
  bubble_data_list[[kriterium]] <- df_table
  
  icc_res <- icc(data.frame(var1, var2), model = "twoway", type = "consistency", unit = "single")
  iccs[[kriterium]] <- icc_res$value
  
  kappa_res <- kappa2(data.frame(var1, var2), weight = "squared")
  kappas[[kriterium]] <- kappa_res$value
}

# --- Scale bubble sizes globally ---
all_bubbles <- do.call(rbind, bubble_data_list)
max_freq <- max(all_bubbles$Freq)

# --- Create balloon plots (subtitle shows ICC & Kappa) ---
plots <- list()
for (i in seq_along(kriterien)) {
  kriterium <- kriterien[i]
  df_table  <- bubble_data_list[[kriterium]]
  
  subtitle_text <- paste0("ICC = ", round(iccs[[kriterium]], 2),
                          " | Kappa = ", round(kappas[[kriterium]], 2))
  
  p <- ggplot(df_table, aes(x = x, y = y, size = Freq)) +
    geom_point(alpha = 0.7, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
    geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey50") +
    geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey50") +
    scale_size_continuous(range = c(3, 12), limits = c(1, max_freq)) +
    coord_fixed(xlim = c(1, 5), ylim = c(1, 5)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none") +
    labs(
      title = kriterium,
      subtitle = subtitle_text,
      x = "Rater 1",
      y = "Rater 2",
      size = "Frequency"
    )
  plots[[i]] <- p
}

# --- Extract a shared legend (uses 'Gesamtrating' subset) ---
legend_base_plot <- ggplot(bubble_data_list[["Gesamtrating"]], aes(x = x, y = y, size = Freq)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  scale_size_continuous(name = "Frequency", range = c(3, 12), limits = c(1, max_freq)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "left")
legend <- get_legend(legend_base_plot)

# --- Arrange plots + legend (2 rows x 3 cols incl. legend slot) ---
grid.arrange(
  grobs = c(plots, list(legend)),
  nrow = 2,
  ncol = 3,
  top  = "Interrater Analysis: Balloon Plots with ICC and Kappa"
)

# ============================================================
# (B) Heatmaps: mean deviation and frequency of strong deviations
# ============================================================

# --- Import data for deviation analysis (generic path) ---
file_path_dev <- "data/inter_rater_final.xlsx"
df <- read_excel(file_path_dev, sheet = sheet_name)

# --- Compute absolute deviations per criterion (|rater1 - rater2|) ---
df <- df %>%
  mutate(
    Abw_Kanal      = abs(Kanal_1 - Kanal_2),
    Abw_Kohaerenz  = abs(Kohaerenz_1 - Kohaerenz_2),
    Abw_Herausgabe = abs(Herausgabe_1 - Herausgabe_2),
    Abw_Auswahl    = abs(Auswahlverfahren_1 - Auswahlverfahren_2),
    Abw_Gesamt     = abs(Gesamtrating_1 - Gesamtrating_2)
  )

# --- Long format across all criteria ---
abweichungen_long <- df %>%
  select(`Gutachter:in 1`, `Gutachter:in 2`,
         Abw_Kanal, Abw_Kohaerenz, Abw_Herausgabe, Abw_Auswahl, Abw_Gesamt) %>%
  pivot_longer(cols = starts_with("Abw_"),
               names_to = "Kriterium", values_to = "Abweichung") %>%
  rename(Gutachter1 = `Gutachter:in 1`, Gutachter2 = `Gutachter:in 2`) # labels for rater 1 and rater 2

# --- Ensure symmetry (A-B equals B-A) ---
abweichungen_pairs <- bind_rows(
  abweichungen_long,
  abweichungen_long %>% rename(Gutachter1 = Gutachter2, Gutachter2 = Gutachter1)
)

# --- (1) Mean deviation per pair ---
pair_stats_mean <- abweichungen_pairs %>%
  group_by(Kriterium, Gutachter1, Gutachter2) %>%
  summarise(Abw_Mean = mean(Abweichung, na.rm = TRUE), .groups = "drop") %>%
  filter(Gutachter1 != Gutachter2)

# --- (2) Frequency of strong deviations (>= 2 points) per pair ---
pair_stats_freq <- abweichungen_pairs %>%
  group_by(Kriterium, Gutachter1, Gutachter2) %>%
  summarise(Abw_Freq = sum(Abweichung >= 2, na.rm = TRUE), .groups = "drop") %>%
  filter(Gutachter1 != Gutachter2)

# --- Plot: mean deviation ---
p1 <- ggplot(pair_stats_mean, aes(x = Gutachter1, y = Gutachter2, fill = Abw_Mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Abw_Mean, 1)), size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~Kriterium, nrow = 2, ncol = 3) +
  labs(title = "Mean deviation between raters per criterion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Plot: frequency of strong deviations (>= 2) ---
p2 <- ggplot(pair_stats_freq, aes(x = Gutachter1, y = Gutachter2, fill = Abw_Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Abw_Freq), size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~Kriterium, nrow = 2, ncol = 3) +
  labs(title = "Frequency of strong deviations (≥ 2 points) per criterion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Print heatmaps ---
print(p1)
print(p2)

# ============================================================
# (C) Recommendations & network of strong deviations
# ============================================================

# --- Recompute pair stats (kept explicit for clarity) ---
pair_stats_mean <- abweichungen_pairs %>%
  group_by(Kriterium, Gutachter1, Gutachter2) %>%
  summarise(Abw_Mean = mean(Abweichung, na.rm = TRUE), .groups = "drop") %>%
  filter(Gutachter1 != Gutachter2)

pair_stats_freq <- abweichungen_pairs %>%
  group_by(Kriterium, Gutachter1, Gutachter2) %>%
  summarise(Abw_Freq = sum(Abweichung >= 2, na.rm = TRUE), .groups = "drop") %>%
  filter(Gutachter1 != Gutachter2)

# --- Merge stats ---
pair_stats <- pair_stats_mean %>%
  left_join(pair_stats_freq, by = c("Kriterium", "Gutachter1", "Gutachter2"))

# --- Thresholds for recommendations ---
schwelle_mean <- 1.5
schwelle_freq <- 3

# --- Label recommendations (traffic-light) ---
pair_stats <- pair_stats %>%
  mutate(Empfehlung = case_when(
    Abw_Mean >= schwelle_mean & Abw_Freq >= schwelle_freq ~ "Rot",
    Abw_Mean >= schwelle_mean | Abw_Freq >= schwelle_freq ~ "Gelb",
    TRUE ~ "Grün"
  )) %>%
  mutate(Empfehlung = factor(Empfehlung, levels = c("Rot", "Gelb", "Grün")))

cat("\n--- Empfehlungen (pairs of raters) ---\n")
print(pair_stats %>% arrange(Kriterium, desc(Abw_Mean), desc(Abw_Freq)))

# --- Heatmap with traffic-light categories ---
ggplot(pair_stats, aes(x = Gutachter1, y = Gutachter2, fill = Empfehlung)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0("Ø=", round(Abw_Mean, 1), "\nH=", Abw_Freq)), size = 3) +
  scale_fill_manual(values = c("Rot" = "red", "Gelb" = "gold", "Grün" = "lightgreen"), drop = FALSE) +
  facet_wrap(~Kriterium, nrow = 2, ncol = 3) +
  labs(
    title = "Urgency for discussion by criterion",
    fill  = "Recommendation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Build network from pairs with "Rot" recommendation ---
edges <- pair_stats %>%
  filter(Empfehlung == "Rot") %>%
  select(Gutachter1, Gutachter2, Empfehlung, Abw_Mean, Abw_Freq) %>%
  distinct()

# Undirected graph (pairs are symmetric)
g <- graph_from_data_frame(d = edges, directed = FALSE)

# --- Community detection (Louvain) using mean deviation as weight ---
comm <- cluster_louvain(g, weights = E(g)$Abw_Mean)

cat("\n--- Community membership of raters ---\n")
print(membership(comm))

# --- Visualize network ---
set.seed(123)  # reproducible layout
plot(
  comm, g,
  vertex.size = 30,
  vertex.label.cex = 0.9,
  vertex.color = "lightpink",
  edge.color = "red",
  edge.width = 3,
  main = "Rater disagreement network (strong deviations)"
)