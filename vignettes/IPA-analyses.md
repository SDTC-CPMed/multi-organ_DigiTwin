## Ingenuity pathway analysis

IPA is a commercial software, but you can request a free trial here
(<https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/>).

To reproduce the analyses for generation of curation file, pathway
analysis, and prediction of upstream regulators (URs) from our project
([CIA analysis](CIA-analyses.html) and [IMIDs
analysis](IMIDs-analyses.html)), follow this pipeline.

1.  Create a project (eg. CIA) in the Project Manager to upload the
    lists of genes. Choose the list of genes of interest, dependent on
    which part of the analysis. For generation of curation-file for
    prediction of inter-organ interactions, include the genes predicted
    as URs by the NicheNet analysis, or all genes from the main data
    input. For pathway analysis or UR prediction, include the DEGs
    (differentially expressed genes) of interest. *Note that IPA has a
    limitation to 5000 genes in the input-file. If >5000 significant
    DEGs has been identified, they need to be prioritized. In this
    study, we never reached >5000 significant DEGs, so no prioritization
    was required.*

<img src="../vignettes/figures/ipa1.png" style="width:45.0%" />

1.  Upload the genes of interest into the project “Dataset Files”. For
    pathway and UR prediction, include the DEGs corresponding LogFCs
    (and q-values if available). Based on these data, choose the ID
    “Mouse gene symbol” or “Human gene symbol”, and the observation
    names “Expr Log Ratio” (LogFC) and “Expr False Discovery Rate”
    (q-val). Keep both LogFC and q-val as the same group, “observation
    1”. ”Save” and name the dataset.

<img src="../vignettes/figures/ipa2.png" style="width:65.0%" />

### Generate curation file

1.  In the top tab tools, the “Location” of the molecules/genes is
    specified. The results can be downloaded by clicking
    <img src="../vignettes/figures/ipa3.png" style="width:4.0%" /> . The
    location of the molecules/genes were in our study confirmed using
    the [Human Protein Atlas](https://www.proteinatlas.org/).

<img src="../vignettes/figures/ipa4.png" style="width:65.0%" />

For downstream codes to run smoothly, save the output file to
“../data/IPA/curation/curation_file.txt”.

### Pathway analysis and prediction of URs

**For each list of DEGs, perform step 3 - 9, for pathway analysis and
prediction of URs in IPA.**

1.  In the lower right corner, click “Analyze/Filter Dataset” and then
    “Core Analysis” to perform IPA analysis of the data.

<img src="../vignettes/figures/ipa5.png" style="width:65.0%" />

1.  Click ”next”, to get to the settings.

<img src="../vignettes/figures/ipa6.png" style="width:45.0%" />

1.  In the settings. based on this dataset, define “General settings -
    Species” = Mouse or Human, “Node Types” = All, “Data Sources” = All,
    “Tissues&Cell Lines” = All, and “Mutation” = All.

2.  Run the analyses by “Run Analysis”.

<img src="../vignettes/figures/ipa7.png" style="width:65.0%" />

1.  All the performed analyses can be found in the “CIA” project under
    “Analyses”. To see and export the results for further analyses,
    choose your current analysis.

<img src="../vignettes/figures/ipa8.png" style="width:40.0%" />

1.  In the top tab tools, go to “Upstream Analysis” to show the UR
    prediction results. The results can be downloaded by clicking
    <img src="../vignettes/figures/ipa3.png" style="width:4.0%" /> .

<img src="../vignettes/figures/ipa9.png" style="width:65.0%" />

For downstream codes to run smoothly, save the output files in
“../data/IPA/pathway_analysis/”, named according to following pattern
**???**.

1.  Next, go to “Upstream Analysis”, in the top tab tools, to show the
    UR prediction results. The results can be downloaded by clicking
    <img src="../vignettes/figures/ipa3.png" style="width:4.0%" /> .

<img src="../vignettes/figures/ipa10.png" style="width:65.0%" />

For downstream codes to run smoothly, save the output files in
“../data/IPA/pathway_analysis/”, named according to following pattern
**???**.

**Add example output table from IPA, and potential code to make readable
by R**
