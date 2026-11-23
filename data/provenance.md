# Inputs and provenance

| Role | Inputs used in the chapter |
|---|---|
| Reference panel | `Reference/ref5.rds` and matching `ref5_exact_demographic.vcf.gz`. These are observed genotype inputs, not the imputed variants. |
| Landscape | `data/derived/lampro_landscape_background.rds`; present CHELSA (1981–2010), historical climate inputs and archived effective-resistance definitions and spatial folds. Do not infer a missing version or GCM. |
| Sampling | Archived ABS geometry and NVIS 7.0 vegetation, plus exact coordinates retained outside this package. |
| Structure and phylogeny | Saved ancestry coefficients, PCoA scores, reference-aligned NEXUS, PAUP block, Newick and split-support summaries. Individual inputs and the sample key are not distributed here. |
| TreeMix | Curated 14-group allele-count input, `lampro_m4_r1.*`, likelihoods and TreeMix plotting functions. The aggregate Q and likelihoods are staged. |
| Demography | Baseline projected folded SFS, archived model parameters and fitted spectra. Candidate A needs only the staged A model/cell summaries and the palette. |
| Diversity | Observed reference/background genotype panels, Q, curated individual metadata and TreeMix crosswalk. Group summaries and pairwise matrices are included in the results directory. |

[File checksums and original paths](../docs/artifact_provenance.tsv) identify the copies used here. See [data availability](../docs/data_availability.md) to request inputs.

