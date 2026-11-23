# Figure sources and captions

The seven main PNGs match the embedded images in the local R7 chapter copy. The README uses their [thumbnails](../assets/thumbnails/) for browsing; use the PDF for the full plate.

| Figure | Exports | Rendering source | Caption |
|---|---|---|---|
| 1 | [PDF](../results/figures/Figure_1_geographic_context.pdf) / [PNG](../results/figures/Figure_1_geographic_context.png) | [source](../analysis/01_genomic_filtering/Fig1_Geographic_context_v3.R) | Sampling points sized by number of individuals; overlapping map panels include 12 common localities. NVIS depicts modern vegetation, not an inferred genetic barrier. |
| 2 | [PDF](../results/figures/Figure_2_nuclear_structure.pdf) / [PNG](../results/figures/Figure_2_nuclear_structure.png) | [source](../analysis/02_population_structure/Fig2_Nuclear_structure.R) | PCoA and broad nuclear sNMF structure at K = 4. Bars summarize individuals; maps display locality means. |
| 3 | [PDF](../results/figures/Figure_3_svdquartets.pdf) / [PNG](../results/figures/Figure_3_svdquartets.png) | [source](../analysis/03_phylogenetics/Fig3_SVDQuartets_geography_v4.R) | Unrooted tree from 170 individuals and 14,766 sites. Bootstrap labels mark ≥70% and support dots ≥95%. Colors are post hoc ancestry annotation. |
| 4 | [PDF](../results/figures/Figure_4_contemporary_gdm.pdf) / [PNG](../results/figures/Figure_4_contemporary_gdm.png) | [source](../analysis/05_contemporary_GDM/Contemporary_GDM_revised.R) | Contemporary GDM and spatial holdouts. BIO4 axes display °C; the fitted curves were unchanged. |
| 5 | [PDF](../results/figures/Figure_5_palaeoclimate.pdf) / [PNG](../results/figures/Figure_5_palaeoclimate.png) | [source](../analysis/06_palaeoclimate/Figure5_palaeoclimate_final.R) | Conditional PaleoMean support across spatial folds and mean squared-error differences by locality. Not a map of divergence dates. |
| 6 | [PDF](../results/figures/Figure_6_treemix.pdf) / [PNG](../results/figures/Figure_6_treemix.png) | [source](../analysis/07_TreeMix/TreeMix_revised.R) | Likelihoods from ten runs per m and the archived m = 4 graph, model lampro_m4_r1; displayed edges are not established events. |
| 7 | [PDF](../results/figures/Figure_7_baseline_demography.pdf) / [PNG](../results/figures/Figure_7_baseline_demography.png) | [source](../analysis/08_demography/plot_figure7_archived.R) | ΔAIC for IM versus SC and observed minus expected mass for the fitted SC folded SFS. Grey denotes masked bins. The colour scale uses asinh spacing with labels in the original residual units. |

The [baseline demographics supplement](../results/figures/Figure_S_demography_baseline.pdf) does not include exploratory expanded-model fits. Consult [input provenance](../data/provenance.md) before running a renderer. For Figure 7, the [source table](../results/summary_tables/moments/Figure7_A_model_and_SFS_cells.csv) retains the original values and cell mask.

## Supplementary figures

I retain the archive filenames; working S numbers are not a new numbering decision. Inputs are identified for provenance and may require permission to access.

| Topic | PDF | Source and inputs |
|---|---|---|
| Three-axis PCoA | [Figure](../results/figures/supplementary/FigS_PCoA_three_axes.pdf) | [Figure 2 renderer](../analysis/02_population_structure/Fig2_Nuclear_structure.R); PCoA scores and eigenvalues |
| Nested sNMF selection | [Figure](../results/figures/supplementary/FigureS1_Nested_sNMF_selection.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); selection and replicate summaries |
| Nested sNMF ancestry | [Figure](../results/figures/supplementary/FigureS1b_Nested_sNMF_ancestry.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); nested ancestry assignments |
| Linkage disequilibrium | [Figure](../results/figures/supplementary/FigureS2_LD_audit.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); LD-by-distance summary |
| Climate selection stability | [Figure](../results/figures/supplementary/FigureS3_Climate_stability.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); climate selection stability summary |
| Historical fold sensitivity | [Figure](../results/figures/supplementary/FigureS4_Historical_fold_sensitivity.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); corrected leave-one-fold-out summary |
| Tasmania sensitivity | [Figure](../results/figures/supplementary/FigureS4b_Tasmania_sensitivity.pdf) | [Archive renderer](../analysis/09_figures/lampro_manuscript_figures.R); Tasmania RMSE sensitivity |
| Separate LIG and LGM tests | [Figure](../results/figures/supplementary/Supplementary_Figure_S7_isolated_LIG_LGM_tests.pdf) | [Archive renderer](../analysis/09_figures/plot_isolated_paleoclimate_tests.R); isolated-period null summaries |
| Random versus spatial validation | [Figure](../results/figures/supplementary/Figure_S_Random_vs_Spatial_CrossValidation.pdf) | [Archive renderer](../analysis/09_figures/plot_heredity.R); [national](../results/summary_tables/validation/national_cv_metrics.csv) and [regional](../results/summary_tables/validation/regional_matched_cv_metrics.csv) summaries |
| TreeMix residual covariance | [Figure](../results/figures/supplementary/treemix_residuals_m4.pdf) | [TreeMix renderer](../analysis/07_TreeMix/TreeMix_revised.R); empirical and model covariance for `lampro_m4_r1` |

I did not recover the global/nested SVDQuartets S1c export from the local or Drive chapter packages. The [joint moments 4D checks](moments_4d_status.md) are not a demographic figure.

