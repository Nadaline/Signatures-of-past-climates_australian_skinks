#!/usr/bin/env Rscript
# Nested sNMF analyses within broad genetic clusters
# Estimate fine-scale structure separately within selected components of the global K = 4 solution. Comparing nested solutions helps describe regional variation that is not visible at the range-wide scale.
# Keep individual order identical between the genotype matrix and ancestry table. Evaluate cross-entropy, replicate agreement and component sizes together when interpreting a nested K.
#

#### Packages ####

suppressPackageStartupMessages({library(LEA); library(mclust); library(ggplot2)})

#### Directories ####

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
parent_project <- normalizePath(file.path(project, ".."), mustWork = TRUE)
out <- file.path(project, "results", "structure", "hierarchical_snmf")
projects_dir <- file.path(out, "lea_projects")
figures <- file.path(out, "figures")
dir.create(projects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
set.seed(20260813)
n_replicates <- 10L
n_cpu <- min(n_replicates, max(1L, parallel::detectCores() - 1L))
q_threshold <- 0.70

#### Load global K4 structure ####

q4 <- read.csv(file.path(project, "tables", "snmf_selected_qmatrix_k4.csv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(names(q4), c("id", "K1", "K2", "K3", "K4")), nrow(q4) == 170)
q4$id <- as.character(q4$id)
q4$global_component <- paste0("K", max.col(q4[, c("K1", "K2", "K3", "K4")]))
q4$global_max_q <- apply(q4[, c("K1", "K2", "K3", "K4")], 1, max)
q4$parent_eligible <- q4$global_max_q >= q_threshold
write.csv(q4, file.path(out, "global_k4_ancestry_audit.csv"), row.names = FALSE)

#### Load the matching sNMF genotype input ####

# The global portable LEA project records n = 170 and L = 19,352. LEA's .geno
# format has one locus per line and one individual per character column. The
# columns have no embedded IDs, so their order is validated against the K4 Q
# table before each parent subset is written.
global_geno <- file.path(parent_project, "ref5.geno")
geno_lines <- readLines(global_geno)
stopifnot(length(geno_lines) == 19352, all(nchar(geno_lines) == nrow(q4)))

#### Prepare parent groups ####

parent_counts <- aggregate(id ~ global_component + parent_eligible, q4, length)
names(parent_counts)[3] <- "n_individuals"
write.csv(parent_counts, file.path(out, "hierarchical_snmf_parent_counts.csv"), row.names = FALSE)
parent_sizes <- setNames(vapply(paste0("K", 1:4), function(k) sum(q4$global_component == k & q4$parent_eligible), integer(1)), paste0("K", 1:4))
child_k <- lapply(parent_sizes, function(n) {
  if (n < 15) 1:2 else 1:min(5, floor(n / 4))
})
# I choose the hierarchy only after inspecting entropy, replicate agreement and
# component sizes. K4-child K=3 has the lowest mean entropy but isolates two
# individuals; I retain the stable K=2 representation as the conservative level.
final_child_k <- c(K1 = 1L, K2 = 1L, K3 = 3L, K4 = 2L)

#### Run sNMF within each parent component ####

run_rows <- list()
model_rows <- list()
assignment_rows <- list()

for (parent in names(parent_sizes)) {
  selected <- which(q4$global_component == parent & q4$parent_eligible)
  parent_ids <- q4$id[selected]
  parent_geno <- file.path(projects_dir, paste0(parent, "_Q070.geno"))
  child_project_file <- file.path(projects_dir, paste0(parent, "_Q070.snmfProject"))
  if (file.exists(child_project_file)) {
    child_project <- load.snmfProject(child_project_file)
  } else {
    writeLines(vapply(geno_lines, function(x) paste0(substring(x, selected, selected), collapse = ""), character(1)), parent_geno)
    child_project <- snmf(parent_geno, K = child_k[[parent]], ploidy = 2, entropy = TRUE,
                          repetitions = n_replicates, project = "new", CPU = n_cpu)
  }
  ce_by_k <- list()
  q_by_k <- list()
  for (k in child_k[[parent]]) {
    ce <- cross.entropy(child_project, K = k)
    ce_by_k[[as.character(k)]] <- ce
    q_runs <- lapply(seq_along(ce), function(run) Q(child_project, K = k, run = run))
    q_by_k[[as.character(k)]] <- q_runs
    run_rows[[paste(parent, k)]] <- data.frame(parent = parent, child_K = k, replicate = seq_along(ce), cross_entropy = as.numeric(ce))
    mean_ari <- if (k == 1 || length(q_runs) < 2) 1 else {
      assignments <- lapply(q_runs, max.col)
      pair_ari <- combn(seq_along(assignments), 2, function(z) adjustedRandIndex(assignments[[z[1]]], assignments[[z[2]]]))
      mean(pair_ari, na.rm = TRUE)
    }
    model_rows[[paste(parent, k)]] <- data.frame(parent = parent, n_parent = length(selected), child_K = k,
      mean_cross_entropy = mean(ce), sd_cross_entropy = sd(ce), min_cross_entropy = min(ce),
      mean_replicate_ARI = mean_ari, n_replicates = length(ce), status = "candidate")
  }
  selected_k <- final_child_k[[parent]]
  best_run <- which.min(ce_by_k[[as.character(selected_k)]])
  q_child <- q_by_k[[as.character(selected_k)]][[best_run]]
  colnames(q_child) <- paste0(parent, letters[seq_len(ncol(q_child))])
  assignment_rows[[parent]] <- cbind(q4[selected, c("id", "global_component", "global_max_q")],
    selected_child_K = selected_k, child_component = colnames(q_child)[max.col(q_child)],
    child_max_q = apply(q_child, 1, max), as.data.frame(q_child))
}

#### Inspect model selection ####

runs <- do.call(rbind, run_rows)
models <- do.call(rbind, model_rows)
assignment_columns <- unique(unlist(lapply(assignment_rows, names)))
assignments <- do.call(rbind, lapply(assignment_rows, function(x) {
  missing_columns <- setdiff(assignment_columns, names(x))
  x[missing_columns] <- NA_real_
  x[, assignment_columns, drop = FALSE]
}))
models$cross_entropy_improvement_from_K1 <- ave(models$mean_cross_entropy, models$parent, FUN = function(x) x[1] - x)
models$selected_by_minimum_cross_entropy <- models$mean_cross_entropy == ave(models$mean_cross_entropy, models$parent, FUN = min)
models$selected_final <- models$child_K == final_child_k[models$parent]
models$status <- ifelse(models$selected_final & models$parent %in% c("K1", "K2"), "NO_REPRODUCIBLE_SUBSTRUCTURE",
  ifelse(models$selected_final & models$parent == "K3", "REPRODUCIBLE_SUBSTRUCTURE",
  ifelse(models$selected_final & models$parent == "K4", "REPRODUCIBLE_SUBSTRUCTURE_CONSERVATIVE", "not_selected")))
write.csv(runs, file.path(out, "hierarchical_snmf_runs.csv"), row.names = FALSE)
write.csv(models, file.path(out, "hierarchical_snmf_model_selection.csv"), row.names = FALSE)
write.csv(assignments, file.path(out, "hierarchical_snmf_assignments.csv"), row.names = FALSE)
write.csv(data.frame(parent = names(final_child_k), selected_child_K = as.integer(final_child_k), decision = c("no child split", "no child split", "three reproducible child components", "two-child conservative split; child K=3 isolates two individuals")), file.path(out, "hierarchical_snmf_final_hierarchy.csv"), row.names = FALSE)
write.csv(data.frame(status = "NOT_AVAILABLE", reason = "No sample- or locality-level Chapple mtDNA clade crosswalk was found in the project files; nuclear child K was selected independently."), file.path(out, "hierarchical_snmf_chapple_correspondence.csv"), row.names = FALSE)

#### Figures ####

coords <- read.csv(file.path(stage5 <- file.path(project, "results", "landscape", "gdm_definitive"), "validation", "spatial_folds.csv"), stringsAsFactors = FALSE)
coords$id <- as.character(coords$id)
plot_data <- merge(q4[, c("id", "global_component", "global_max_q")], coords[, c("id", "X", "Y")], by = "id", all.x = TRUE)
plot_data <- merge(plot_data, assignments[, c("id", "child_component", "child_max_q")], by = "id", all.x = TRUE)
plot_data$hierarchical_component <- ifelse(is.na(plot_data$child_component), plot_data$global_component, plot_data$child_component)
map <- ggplot(plot_data, aes(X, Y, colour = hierarchical_component, shape = global_component)) +
  geom_point(size = 2.8, alpha = .9) + coord_equal() + theme_classic() +
  labs(x = "Longitude", y = "Latitude", colour = "Nested nuclear component", shape = "Global K4 component")
ggsave(file.path(figures, "hierarchical_snmf_map.png"), map, width = 7.5, height = 7, dpi = 300)

#### Save summary ####

summary_lines <- c(
  "# Hierarchical sNMF", "",
  "I used the established global K4 Q matrix as the broad nuclear structure and fitted child sNMF models only to individuals with maximum global ancestry Q >= 0.70. The global K4 components remain macrostructure; child components are evaluated independently and are not named after mitochondrial clades.", "",
  "The input genotype file was `../ref5.geno`, matching the global LEA project dimensions of 170 individuals and 19,352 loci. The output tables report cross-entropy across ten replicates and replicate agreement as adjusted Rand index for each parent/child-K combination.", "",
  "I freeze the first nuclear hierarchy as K1: no reproducible child subdivision; K2: no reproducible child subdivision; K3: three reproducible child components; and K4: two-child conservative subdivision. In K4, child K=3 minimized mean cross-entropy but created a component of two individuals; I therefore retain the stable child K=2 solution rather than promote a numerically fragile rare component.", "",
  "I have not used Chapple et al. (2011) to select a child K. A direct sample/locality crosswalk between the historical mtDNA clades and these samples was not present in the project files at the time of this analysis, so no agreement statistic is reported. The correspondence table will remain explicitly unavailable rather than inferred from a map."
)
writeLines(summary_lines, file.path(out, "hierarchical_snmf_summary.md"))
