# Documentation

## Introduction

This documentation explains the methodological background and analysis
workflow for assessing **interrater reliability** in the evaluation of
edited scholarly volumes.\
Use case: the work was conducted in the context of a field study with
expert reviewers, focusing on the quality criteria defined by the AG
Qualitaet in der Forschung (AGQ).

The analysis follows a logical sequence:\
1. Quantify agreement between reviewers for each criterion using ICC and
weighted Kappa.\
2. Explore deviations between reviewer pairs through heatmaps (mean
deviation and frequency of strong deviations).\
3. Provide actionable recommendations for further discussion on
evaluation criteria using traffic-light coding.\
4. Detect structural disagreement patterns through network analysis,
identifying reviewer subgroups and the criteria requiring further
clarification within these groups.

This sequence was chosen to combine **statistical rigor** with
**practical insights** for reviewer panel management. Reliability
coefficients quantify overall consistency, heatmaps highlight
problematic reviewer pairs, traffic-light coding translates statistical
thresholds into intuitive signals, and network analysis reveals clusters
of systematic disagreement where subgroups of reviewers need to align
their interpretation of evaluation criteria.

------------------------------------------------------------------------

## Script: `interrater-reliability.R`

### Step A: Balloon plots with ICC and Kappa

**Balloon plots of interrater agreement (with ICC and weighted
Kappa).**\
Displayed is the rating coherence for 60 edited scholarly volumes, each
assessed independently by two expert reviewers on a five-point scale (1
= not met, 5 = fully met) across defined quality criteria. For each
criterion (and for the overall assessment), cross-tabulations of rating
pairs are visualized as proportional balloons:\
- **Position** in the grid indicates the combination of the two ratings
(Reviewer 1 vs. Reviewer 2).\
- **Balloon size** reflects the relative frequency of that rating pair.\
- **Solid diagonal line** = perfect agreement; **dashed lines** = ±1
point deviation.

To quantify agreement, two complementary measures were calculated:\
- **Intraclass Correlation Coefficient (ICC, model ICC(C,1), McGraw &
Wong 1996):**\
\[ ICC = `\frac{MS_R - MS_E}{MS_R + (k - 1) \cdot MS_E}`{=tex} \]\
where (MS_R) is the variance between objects (volumes), (MS_E) the error
variance, and (k=2) the number of ratings per volume.\
- **Weighted Cohen's Kappa (quadratic weights):**\
\[ `\kappa`{=tex}*w = 1 -
`\frac{\sum_{i,j} w_{ij} \cdot o_{ij}}{\sum_{i,j} w_{ij} \cdot e_{ij}}`{=tex}
\]\
with (o*{ij}) = observed frequency, (e\_{ij}) = expected frequency based
on marginal distributions, and (w\_{ij} =
`\left`{=tex}(`\frac{i - j}{k - 1}`{=tex}`\right`{=tex})\^2) as the
weight of disagreement.

Importantly, rating pairs were evaluated without fixing raters, meaning
individual reviewer identities were not separated. This avoids bias in
marginal distributions due to systematic differences in rating behavior.
Thus, the Kappa values represent the **average level of agreement in the
process under random allocation of reviewers**.

------------------------------------------------------------------------

### Step B: Heatmaps of reviewer deviations

**Purpose:** Identify where reviewer pairs disagree most strongly.

**Heatmaps of pairwise deviations.**\
These plots examine the extent and frequency of disagreement between
reviewer pairs. For each criterion:\
- **Heatmap 1 (mean deviation):** shows the average difference in
ratings between two reviewers.\
- **Heatmap 2 (frequency of strong deviations):** counts how often
deviations of ≥ 2 points occurred.

Cells are colored from white (low disagreement) to red (high
disagreement), with numeric annotations for clarity.\
- **X-axis / Y-axis:** Reviewer 1 vs. Reviewer 2.\
- **Cell color:** magnitude of deviation (mean or frequency).\
- **Cell label:** rounded mean difference or absolute count.

**Interpretation:**\
Red cells indicate reviewer pairs with consistently high or frequent
disagreements, suggesting a need for calibration or closer discussion.
These visualizations provide an intuitive overview of where reliability
problems are concentrated.

------------------------------------------------------------------------

### Step C: Recommendations and traffic-light visualization

**Purpose:** Translate statistical results into actionable signals for
reviewer management.

**Traffic-light heatmaps for reviewer alignment.**\
To transform statistical thresholds into actionable signals, the
following cutoffs were applied:\
- Mean deviation ≥ 1.5 points\
- Frequency of strong deviations ≥ 3 cases (definition of strong
deviation: ≥ 2 points)

Cells are coded as:\
- **Red:** both thresholds exceeded → urgent need for discussion\
- **Yellow:** one threshold exceeded → potential need for clarification\
- **Green:** no thresholds exceeded → no immediate action required

**Output:** Heatmap with categorical colors and overlaid text values
(mean deviation and frequency).

**Interpretation:**\
The traffic-light visualization provides a concise management tool for
reviewer panels. Red-highlighted reviewer pairs should meet to reconcile
differences, yellow pairs may require clarification, and green pairs
indicate satisfactory alignment.

------------------------------------------------------------------------

### Step D: Network analysis of strong disagreements

**Purpose:** Detect structural patterns of disagreement across reviewer
panels.

**Network graph of reviewer disagreement communities.**\
This analysis examines the structural patterns of disagreement by
focusing on reviewer pairs with "red" recommendations.

**Method:**\
- Build an undirected graph where nodes represent reviewers and edges
represent strong disagreements.\
- Edge weight = mean deviation between the pair.\
- Apply Louvain community detection to identify clusters of reviewers.

**Output:** Network visualization with reviewers as nodes, red edges for
strong disagreements, and community clusters highlighted.

**Interpretation:**\
The network graph shows how disagreements are distributed across the
panel. Communities of reviewers with systematic divergences are
revealed, highlighting where group-level discussions are necessary in
addition to bilateral clarifications. This complements the heatmaps by
showing disagreement as a structural property of the reviewer network.

------------------------------------------------------------------------

## Interpretation of results based on our field study with 60 edited books

-   ICC and Kappa values were near zero or negative, indicating poor
    reliability.\
-   Heatmaps confirmed substantial divergence between several reviewer
    pairs.\
-   Traffic-light visualization highlighted specific pairs needing
    urgent alignment.\
-   Network analysis revealed disagreement communities, providing
    guidance for targeted reviewer calibration.

**Conclusion:** The interrater reliability analysis demonstrated that
current quality criteria require clearer operationalization. The
combination of statistical and visual methods provides a robust
framework for diagnosing and addressing weaknesses observed in the field
study.

------------------------------------------------------------------------
