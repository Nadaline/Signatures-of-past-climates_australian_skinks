# Working from figures, processed inputs or reads

## A. Draw a figure using saved values

From the repository root, run `mkdir -p work && Rscript analysis/08_demography/plot_figure7_archived.R work` (base R with Cairo). It reads the [paired comparison and SFS cells](../results/summary_tables/moments/Figure7_A_model_and_SFS_cells.csv) and the shared palette, validates their dimensions and residuals, and writes `Figure_7_baseline_demography.pdf` and `.png` to `work/`. It does not fit a demographic model. Other renderers listed in the [figure index](figures.md) need saved model objects, spatial inputs or individual summaries that are not all included here.

## B. Refit from processed genotypes

Request the [approved analytical inputs](data_availability.md), including the ordered genotype and sample keys. For structure, start with the reference panel and the [nested sNMF source](../analysis/02_population_structure/hierarchical_snmf.R), which reads the existing global ancestry and aligned genotype input and returns nested ancestry and cross-entropy summaries. For the landscape models use the background panel, environmental table and matched resistance matrix. Run the [contemporary GDM](../analysis/05_contemporary_GDM/gdm_definitive_spatial_pipeline.R) first to produce the held-out fits, folds and conditional comparisons; the [historical comparison](../analysis/06_palaeoclimate/gdm_stage6_historical_landscape.R) reads those model and fold objects and produces historical prediction summaries. The two sources accept `LAMPRO_CHELSA_DIR`. Verify that IDs are unique and that matrix row and column orders match the metadata before calculating distances.

The PAUP block needs the approved 14,766-site NEXUS alignment. TreeMix needs the separate allele-count input and geographic assignment rules. The pairwise moments fits require the original projected, folded spectra, masks, model definitions and fit settings; the Figure 7 cell table is not a replacement for them. R with LEA, dartR and the declared GDM/spatial packages, PAUP*, TreeMix and Python moments serve different parts of the study. R7 identifies PAUP* v4.0a build 169, and the joint 4D check records moments 1.6.0. Neither verifies versions for other runs. The [example configuration](../config/example_paths.yml) describes paths; it is not loaded automatically. Record `sessionInfo()` and external tool versions when repeating an analysis.

## C. Reconstruct from reads

The original README records earlier Stacks and alignment commands, but these do not establish the full route to the reference panel used here. Raw reads, barcodes, reference/index provenance and the full filtering record are not supplied together. A fresh checkout is therefore not sufficient for reads-to-results reconstruction.

