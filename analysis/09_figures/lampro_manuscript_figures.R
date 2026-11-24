#!/usr/bin/env Rscript
# Main and supplementary figures for the Lampropholis study
# Assemble geographic, nuclear-structure, landscape, TreeMix and demographic panels from the relevant analysis outputs. Keep the scientific meanings of scales, group labels and model summaries consistent across figures.
# Use saved results for figure construction and refer to the analysis scripts for the fitting procedures. Genetic clusters, geographic groups and mitochondrial clades are different classifications.
#

#### Packages ####

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(cowplot); library(png); library(maps); library(ape); library(ggtree)
})

#### Directories ####

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
out <- file.path(project, "figures", "manuscript")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
figure_suffix <- "_v2"

#### Manuscript theme and project colours ####

palette <- read.csv(file.path(project, "project_palette.csv"), stringsAsFactors = FALSE)
ancestry_colours <- setNames(palette$color[match(paste0("K", 1:4), palette$unit)], paste0("K", 1:4))
nested_colours <- c(K1a = "#1B6B5A", K2a = "#C49A00", K3a = "#C8553D", K3b = "#E58B72",
                    K3c = "#8F3B3F", K4a = "#007C91", K4b = "#78C6C8")
history_colours <- c(negative = "#9A4A43", neutral = "#F3F1EA", positive = "#426B5A")
theme_lampro <- theme_classic(base_size = 10) + theme(
  axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 11),
  plot.subtitle = element_text(size = 8.5), legend.title = element_text(face = "bold"),
  panel.border = element_rect(colour = "grey55", fill = NA, linewidth = .35),
  strip.background = element_rect(fill = "grey94", colour = "grey55"), strip.text = element_text(face = "bold")
)
save_figure <- function(plot, name, width, height) {
  pdf_file <- file.path(out, paste0(name, figure_suffix, ".pdf"))
  png_file <- file.path(out, paste0(name, figure_suffix, ".png"))
  if (file.exists(pdf_file) && file.exists(png_file)) return(invisible(NULL))
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = 600)
}
save_svd_figure <- function(plot, name, width, height) {
  ggsave(file.path(out, paste0(name, "_v3.pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out, paste0(name, "_v3.png")), plot, width = width, height = height, dpi = 600)
}
as_plot_image <- function(path) {
  image <- readPNG(path)
  ggplot() + annotation_raster(image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + coord_fixed() + theme_void()
}
australia <- map_data("world", region = c("Australia", "Tasmania"))

#### Read previously estimated genetic and landscape results ####

pca <- read.csv(file.path(project, "population_structure", "pcoa_scores.csv"), stringsAsFactors = FALSE)
# The previously saved object was generated with dartR::gl.pcoa. Its eigenvalues, rather
# than the plotted coordinates, define the reported PCoA variation percentages.
pcoa_archived <- readRDS(file.path(project, "..", "Reference", "pcoa.ref7.rds"))
pcoa_positive_eigenvalues <- pcoa_archived$eig[pcoa_archived$eig > 0]
pcoa_variance <- 100 * pcoa_positive_eigenvalues / sum(pcoa_positive_eigenvalues)
q4 <- read.csv(file.path(project, "tables", "snmf_selected_qmatrix_k4.csv"), stringsAsFactors = FALSE)
q4$global_component <- paste0("K", max.col(q4[, paste0("K", 1:4)]))
pca$id <- as.character(pca$id); q4$id <- as.character(q4$id)
pca <- left_join(pca, q4[, c("id", "global_component")], by = "id")
gcr <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "tables", "gcr_seven_subsets.csv"))
cv <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "validation", "spatial_cv_model_comparison.csv"))
nulls <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "validation", "spatial_group_null_models.csv"))
null_summary <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "validation", "spatial_group_null_summary.csv"))
isplines <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "tables", "final_isplines.csv"))
historical_cv <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "validation", "historical_model_comparison.csv"))
historical_fold <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "null_corrected", "corrected_observed_fold_effects.csv"))
historical_null <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "null_corrected", "corrected_historical_null_statistics.csv"))
historical_null_summary <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "null_corrected", "corrected_historical_null_summary.csv"))
fold_sites <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "validation", "spatial_folds.csv"))
local_gain <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "spatial_effect", "tables", "gdm_historical_gain_by_locality.csv"))
tasmania_sensitivity <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "spatial_effect", "tables", "gdm_tasmania_rmse_sensitivity.csv"))
treemix <- read.csv(file.path(project, "treemix", "likelihood_summary.csv"))
treemix_runs <- read.csv(file.path(project, "treemix", "likelihood_all_runs.csv"))
stage8b <- read.csv(file.path(project, "results", "demography", "moments", "tables", "stage8b_model_comparison.csv"))
hierarchy <- read.csv(file.path(project, "results", "structure", "hierarchical_snmf", "hierarchical_snmf_model_selection.csv"))
hierarchy_runs <- read.csv(file.path(project, "results", "structure", "hierarchical_snmf", "hierarchical_snmf_runs.csv"))
hierarchy_q <- read.csv(file.path(project, "results", "structure", "hierarchical_snmf", "hierarchical_snmf_assignments.csv"), stringsAsFactors = FALSE)
sample_lineages <- read.csv(file.path(project, "population_structure", "snmf_k4_individual_ancestry_long.csv"), stringsAsFactors = FALSE) %>%
  distinct(id, pop) %>% transmute(id = as.character(id), chapple_lineage = pop)
hierarchy_q <- hierarchy_q %>% rename(nested_K = selected_child_K, nested_component = child_component, nested_max_Q = child_max_q) %>%
  mutate(id = as.character(id)) %>% left_join(sample_lineages, by = "id")
ld_summary <- read.csv(file.path(project, "data", "derived", "ld_distance_summary.csv"), stringsAsFactors = FALSE)
climate_stability <- read.csv(file.path(project, "results", "landscape", "gdm_definitive", "validation", "climate_selection_stability.csv"))
leave_one_fold <- read.csv(file.path(project, "results", "landscape", "paleoclimate", "null_corrected", "corrected_leave_one_fold_out.csv"))
convergence <- read.csv(file.path(project, "results", "demography", "moments", "tables", "stage8b_final_convergence.csv"))

#### Figure 1 - Broad ancestry turnover ####

p1a <- ggplot(pca, aes(PCo1, PCo2, colour = global_component)) +
  geom_point(size = 2.1, alpha = .9) + scale_colour_manual(values = ancestry_colours, drop = FALSE) + theme_lampro +
  labs(title = "Nuclear PCoA", x = sprintf("PCoA axis 1 (%.1f%%)", pcoa_variance[1]),
       y = sprintf("PCoA axis 2 (%.1f%%)", pcoa_variance[2]), colour = "Dominant global ancestry component")

# Alpha blends a low-Q assignment into the beige land background. Thus point
# opacity, not point size, represents the strength of each ancestry coefficient.
q4_map <- left_join(q4, fold_sites[, c("id", "X", "Y")], by = "id") %>%
  group_by(X, Y) %>% summarise(K1 = mean(K1), K2 = mean(K2), K3 = mean(K3), K4 = mean(K4), .groups = "drop") %>%
  pivot_longer(cols = K1:K4, names_to = "component", values_to = "ancestry")
p1b <- ggplot() +
  geom_polygon(data = australia, aes(long, lat, group = group), fill = "#EAE4D9", colour = "grey65", linewidth = .25) +
  geom_point(data = q4_map, aes(X, Y, fill = component, alpha = ancestry), shape = 21, colour = "#EAE4D9", stroke = .12, size = 3.7) +
  coord_quickmap(xlim = c(135, 154), ylim = c(-44, -18)) + facet_wrap(~component, ncol = 2) +
  scale_fill_manual(values = ancestry_colours, guide = "none") + scale_alpha(range = c(.08, 1), limits = c(0, 1)) + theme_lampro +
  labs(title = "Spatial K4 ancestry", x = "Longitude", y = "Latitude", alpha = "Mean Q")
figure1 <- (p1a | p1b) + plot_annotation(tag_levels = "A")
save_figure(figure1, "Figure1_population_structure", 11.5, 4.8)

#### Figure 2 - Contemporary landscape turnover ####

gcr$model <- factor(gcr$model, levels = c("G", "C", "R", "CG", "GR", "CR", "CGR"), labels = c("G", "C", "R", "G + C", "G + R", "C + R", "G + C + R"))
p2a <- ggplot(gcr, aes(model, deviance_explained, fill = model)) + geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("G" = "#3F6F58", "C" = "#9B6A3B", "R" = "#577590", "G + C" = "#6F8C6D", "G + R" = "#477F88", "C + R" = "#8A7860", "G + C + R" = "#2F4858")) +
  theme_lampro + labs(title = "Descriptive full-data G/C/R fit", x = NULL, y = "Explained deviance (%)")
cv_summary <- cv %>% filter(model %in% c("G", "GR", "GC", "GRC")) %>% group_by(model) %>% summarise(mean_rmse = mean(rmse), .groups = "drop")
cv_summary$model <- factor(cv_summary$model, levels = c("G", "GR", "GC", "GRC"), labels = c("G", "G + R", "G + C", "G + R + C"))
p2b <- ggplot(cv_summary, aes(model, mean_rmse)) + geom_col(fill = "#4C6A5D") +
  geom_text(aes(label = sprintf("%.4f", mean_rmse)), vjust = -0.35, size = 3) + theme_lampro +
  expand_limits(y = max(cv_summary$mean_rmse) * 1.09) + labs(title = "Spatially held-out prediction", x = NULL, y = "Mean held-out RMSE")
null_long <- bind_rows(data.frame(group = "Climate", gain = nulls$climate_gain_rmse, permutation = nulls$permutation),
                        data.frame(group = "Resistance", gain = nulls$resistance_gain_rmse, permutation = nulls$permutation)) %>%
  group_by(group, permutation) %>% summarise(gain = mean(gain), .groups = "drop")
null_labels <- data.frame(group = c("Climate", "Resistance"), observed = c(.0068804389, .0005532859),
                          label = c("Observed ΔRMSE = 0.00688\np = 0.141", "Observed ΔRMSE = 0.00055\np = 0.019"))
p2c <- ggplot(null_long %>% filter(permutation > 0), aes(gain)) + geom_histogram(bins = 35, fill = "grey75", colour = "white") +
  geom_vline(data = null_labels, aes(xintercept = observed), colour = "#9A4A43", linewidth = .8) +
  geom_label(data = null_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.1, size = 2.6, label.size = .15) +
  facet_wrap(~group, scales = "free") + theme_lampro + labs(title = "Spatial null distributions", x = "Conditional held-out RMSE gain", y = "Permutations")
ispline_long <- bind_rows(lapply(c("Geographic", "tempRange", "tempSeason", "precAnn", "topographicResistance"),
  function(v) data.frame(x = isplines[[paste0("x_", v)]], y = isplines[[paste0("y_", v)]], predictor = v)))
ispline_long$predictor <- factor(ispline_long$predictor, levels = c("Geographic", "topographicResistance", "precAnn", "tempRange", "tempSeason"),
  labels = c("Geographic separation (decimal degrees)", "Effective topographic resistance (Circuitscape units)",
             "Annual precipitation (mm)", "Mean diurnal temperature range (°C)", "Temperature seasonality (BIO4 CHELSA units)"))
p2d <- ggplot(ispline_long, aes(x, y)) + geom_line(colour = "#2F4858", linewidth = .8) +
  facet_wrap(~predictor, scales = "free_x", ncol = 2) + theme_lampro + labs(title = "Fitted GDM I-splines", x = "Predictor distance", y = "Accumulated fitted turnover")
figure2 <- (p2a | p2b) / (p2c | p2d) + plot_annotation(tag_levels = "A")
save_figure(figure2, "Figure2_contemporary_GDM", 12, 9)

#### Figure 3 - Historical landscape signal ####

historical_summary <- historical_cv %>% filter(model %in% c("M2", "M3_PaleoMean")) %>% group_by(model) %>% summarise(mean_rmse = mean(rmse), .groups = "drop")
historical_summary$model <- factor(historical_summary$model, levels = c("M2", "M3_PaleoMean"), labels = c("M2: present landscape", "M3: + PaleoMean"))
p3a <- ggplot(historical_summary, aes(model, mean_rmse, fill = model)) + geom_col(show.legend = FALSE, fill = c("#627D98", "#4C6A5D")) + theme_lampro +
  annotate("text", x = 1.5, y = max(historical_summary$mean_rmse) * 1.04, label = "Equal-fold ΔRMSE = 0.002739 (2.2%)", fontface = "bold", size = 3.3) +
  labs(title = "Palaeoclimate improves held-out prediction", x = NULL, y = "Mean held-out RMSE") + theme(axis.text.x = element_text(size = 7))
p3b <- ggplot(historical_fold, aes(factor(fold), observed_Delta_RMSE, fill = observed_Delta_RMSE > 0)) + geom_hline(yintercept = 0, linetype = 2) + geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = history_colours[["positive"]], `FALSE` = history_colours[["negative"]])) + theme_lampro + labs(title = "Historical gain by spatial fold", x = "Held-out fold", y = "Equal-fold ΔRMSE (M2 − M3)")
p3c <- ggplot(historical_null, aes(Delta_RMSE_null)) + geom_histogram(bins = 35, fill = "grey75", colour = "white") +
  geom_vline(xintercept = historical_null_summary$observed_Delta_RMSE, colour = "#9A4A43", linewidth = .9) + theme_lampro +
  labs(title = paste0("Corrected spatial null: p = ", historical_null_summary$p_value), x = "Mean equal-fold ΔRMSE", y = "999 permutations")
p3d <- ggplot() + geom_polygon(data = australia, aes(long, lat, group = group), fill = "#EAE4D9", colour = "grey65", linewidth = .25) +
  geom_point(data = fold_sites, aes(X, Y, colour = factor(fold)), size = 2.3, alpha = .9) + coord_quickmap(xlim = c(135, 154), ylim = c(-44, -18)) + theme_lampro +
  scale_colour_manual(values = c("1" = "#B66B3E", "2" = "#577590", "3" = "#4C6A5D", "4" = "#9A4A43", "5" = "#7E7E7E")) + labs(title = "Frozen spatial folds", x = "Longitude", y = "Latitude", colour = "Fold")
p3e <- ggplot() + geom_polygon(data = australia, aes(long, lat, group = group), fill = "#EAE4D9", colour = "grey65", linewidth = .25) +
  geom_point(data = local_gain, aes(X, Y, colour = mean_pairwise_gain_squared_error, size = abs(mean_pairwise_gain_squared_error)), alpha = .9) + coord_quickmap(xlim = c(135, 154), ylim = c(-44, -18)) +
  scale_colour_gradient2(low = history_colours[["negative"]], mid = history_colours[["neutral"]], high = history_colours[["positive"]], midpoint = 0) + theme_lampro +
  labs(title = "Locality-level held-out historical gain", subtitle = "Descriptive pooled pair-weighted diagnostic; not the equal-fold test statistic", x = "Longitude", y = "Latitude", colour = "Mean gain", size = "Magnitude")
p3s <- ggplot(tasmania_sensitivity, aes(diagnostic, Delta_RMSE, fill = diagnostic)) + geom_hline(yintercept = 0, linetype = 2) + geom_col(show.legend = FALSE) + theme_lampro +
  scale_fill_manual(values = c("all_heldout_pairs" = "#4C6A5D", "without_Tasmania_pairs" = "#6F8C6D", "without_Tasmania_containing_fold" = "#8AA88A")) + labs(x = NULL, y = "Pooled pair-weighted ΔRMSE", title = "Historical-signal sensitivity") + theme(axis.text.x = element_text(angle = 20, hjust = 1))
figure3 <- (p3a | p3b) / (p3c | p3d) / p3e + plot_annotation(tag_levels = "A")
save_figure(figure3, "Figure3_historical_landscape", 12, 12)
save_figure(p3s, "Supplementary_Figure_historical_Tasmania_sensitivity", 7, 4.5)

#### Figure 4 - Reticulate history ####

p4a <- ggplot(treemix_runs, aes(factor(m), loglik)) + geom_point(alpha = .18, colour = "#2F4858", position = position_jitter(width = .12, height = 0)) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#2F4858", linewidth = .8) + stat_summary(fun = mean, geom = "point", colour = "#2F4858", size = 2.2) +
  stat_summary(data = subset(treemix_runs, m == 4), fun = mean, geom = "point", colour = "#9A4A43", size = 3.2) + theme_lampro +
  labs(title = "TreeMix likelihood trajectory", subtitle = "Points are ten replicate fits; m = 4 is a parsimonious working graph", x = "Migration edges (m)", y = "Log-likelihood")
p4b <- as_plot_image(file.path(project, "population_structure", "05_treemix_best_m4.png"))
figure4 <- (p4a | p4b) + plot_annotation(tag_levels = "A")
save_figure(figure4, "Figure4_treemix", 12, 5)

#### Figure 5 - Demographic model comparison ####

stage8b$model <- factor(stage8b$model, levels = c("SI", "IM", "AM", "SC"))
p5a <- ggplot(stage8b, aes(model, delta_AIC, fill = model)) + geom_col(show.legend = FALSE) + facet_wrap(~pair, scales = "free_y") +
  scale_fill_manual(values = c(SI = "#7E7E7E", IM = "#577590", AM = "#B66B3E", SC = "#4C6A5D")) + theme_lampro +
  labs(title = "Constant-size demographic comparison", subtitle = "ΔAIC is comparable within, not among, population pairs", x = NULL, y = "ΔAIC")
sfs_scaling <- lapply(c("K1_K3", "K1_K4", "K2_K4"), function(pair) {
  observed <- as.matrix(read.csv(file.path(project, "results", "demography", "moments", "figures", "final_model_fits", paste0(pair, "_SC_observed.csv")), header = FALSE))
  expected <- as.matrix(read.csv(file.path(project, "results", "demography", "moments", "figures", "final_model_fits", paste0(pair, "_SC_expected.csv")), header = FALSE))
  factor <- sum(observed) / sum(expected); expected_scaled <- expected * factor
  data.frame(pair = pair, sum_observed = sum(observed), sum_raw_expected = sum(expected), scaling_factor = factor, sum_scaled_expected = sum(expected_scaled))
}) %>% bind_rows()
write.csv(sfs_scaling, file.path(project, "results", "demography", "moments", "tables", "stage8b_sfs_plot_scaling_audit.csv"), row.names = FALSE)
matrix_plot <- function(pair, type, title) {
  observed <- as.matrix(read.csv(file.path(project, "results", "demography", "moments", "figures", "final_model_fits", paste0(pair, "_SC_observed.csv")), header = FALSE))
  expected <- as.matrix(read.csv(file.path(project, "results", "demography", "moments", "figures", "final_model_fits", paste0(pair, "_SC_expected.csv")), header = FALSE))
  expected <- expected * (sum(observed) / sum(expected))
  mat <- if (type == "observed") observed else if (type == "expected") expected else observed - expected
  long <- expand.grid(i = seq_len(nrow(mat)), j = seq_len(ncol(mat))); long$value <- as.vector(mat)
  ggplot(long, aes(j, i, fill = value)) + geom_raster() + coord_equal() + theme_void() + ggtitle(title) +
    scale_fill_gradient2(low = "#577590", mid = "white", high = "#9A4A43", midpoint = if (type == "residual") 0 else median(long$value, na.rm = TRUE))
}
sfs_panels <- lapply(c("K1_K3", "K1_K4", "K2_K4"), function(pair) {
  pair_lab <- gsub("_", " vs ", pair)
  matrix_plot(pair, "observed", paste(pair_lab, "observed")) | matrix_plot(pair, "expected", paste(pair_lab, "expected, scaled")) | matrix_plot(pair, "residual", paste(pair_lab, "observed − scaled expected"))
})
figure5 <- p5a / wrap_plots(sfs_panels, ncol = 1) + plot_annotation(title = "Stage 8B folded-SFS model fit (SC shown; Stage 8C robustness pending)", tag_levels = "A")
save_figure(figure5, "Figure5_demographic_model_comparison", 12, 15)

#### Supplementary Figure S1 - Nested sNMF selection ####

hierarchy <- hierarchy %>% rename(global_component = parent, nested_K = child_K, n_individuals = n_parent)
hierarchy_runs <- hierarchy_runs %>% rename(global_component = parent, nested_K = child_K)
hierarchy$global_component <- factor(hierarchy$global_component, levels = paste0("K", 1:4))
hierarchy_runs$global_component <- factor(hierarchy_runs$global_component, levels = paste0("K", 1:4))
selection_points <- hierarchy %>% filter(selected_final) %>% transmute(global_component, nested_K, mean_cross_entropy)
pS1a <- ggplot(hierarchy_runs, aes(nested_K, cross_entropy)) + geom_point(colour = "grey55", alpha = .32, position = position_jitter(width = .06, height = 0)) +
  geom_errorbar(data = hierarchy, aes(nested_K, ymin = mean_cross_entropy - sd_cross_entropy, ymax = mean_cross_entropy + sd_cross_entropy), inherit.aes = FALSE, width = .12, colour = "#2F4858") +
  geom_line(data = hierarchy, aes(nested_K, mean_cross_entropy, group = global_component), inherit.aes = FALSE, colour = "#2F4858") + geom_point(data = hierarchy, aes(nested_K, mean_cross_entropy), inherit.aes = FALSE, colour = "#2F4858", size = 1.8) +
  geom_point(data = selection_points, aes(nested_K, mean_cross_entropy), inherit.aes = FALSE, shape = 21, fill = "#9A4A43", colour = "black", size = 3) +
  facet_wrap(~global_component, scales = "free_y", nrow = 1, labeller = labeller(global_component = function(x) paste0("Global component ", x))) + theme_lampro +
  labs(title = "Nested sNMF cross-entropy selection", subtitle = "Light points: replicate fits; bars: ±1 SD; red circle: retained nested K", x = "Candidate nested K", y = "Cross-entropy")
pS1b <- ggplot(hierarchy, aes(nested_K, mean_replicate_ARI)) + geom_line(colour = "#2F4858") + geom_point(colour = "#2F4858", size = 1.8) +
  geom_point(data = selection_points %>% left_join(hierarchy[, c("global_component", "nested_K", "mean_replicate_ARI")], by = c("global_component", "nested_K")), aes(nested_K, mean_replicate_ARI), inherit.aes = FALSE, shape = 21, fill = "#9A4A43", colour = "black", size = 3) +
  facet_wrap(~global_component, scales = "free_y", nrow = 1, labeller = labeller(global_component = function(x) paste0("Global component ", x))) + theme_lampro +
  labs(title = "Nested sNMF replicate stability", subtitle = "Red circle: retained nested K", x = "Candidate nested K", y = "Mean adjusted Rand index")
save_figure(pS1a / pS1b, "Supplementary_Figure_S1_nested_snmf_K_selection", 14, 7)

nested_decisions <- hierarchy %>% transmute(global_component, n_individuals, tested_nested_K = nested_K, mean_cross_entropy,
  cross_entropy_change = cross_entropy_improvement_from_K1, replicate_stability = mean_replicate_ARI,
  selected_nested_K = ifelse(selected_final, nested_K, NA_integer_), decision = ifelse(status == "NO_REPRODUCIBLE_SUBSTRUCTURE", "NO_REPRODUCIBLE_NESTED_STRUCTURE", ifelse(selected_final, "REPRODUCIBLE_NESTED_STRUCTURE", "NOT_SELECTED")))
write.csv(nested_decisions, file.path(out, "nested_snmf_selection_decisions.csv"), row.names = FALSE)

#### Supplementary Figure S1b - Nested nuclear ancestry ####

hierarchy_plot <- left_join(hierarchy_q, fold_sites[, c("id", "X", "Y")], by = "id")
hierarchy_plot$order <- ave(hierarchy_plot$Y, hierarchy_plot$global_component, FUN = function(x) rank(x, ties.method = "first"))
hierarchy_long <- bind_rows(
  hierarchy_plot %>% filter(global_component == "K1") %>% transmute(id, global_component = "K1", order, nested_component = "K1a", Q = K1a),
  hierarchy_plot %>% filter(global_component == "K2") %>% transmute(id, global_component = "K2", order, nested_component = "K2a", Q = K2a),
  hierarchy_plot %>% filter(global_component == "K3") %>% select(id, order, K3a, K3b, K3c) %>% pivot_longer(K3a:K3c, names_to = "nested_component", values_to = "Q") %>% mutate(global_component = "K3"),
  hierarchy_plot %>% filter(global_component == "K4") %>% select(id, order, K4a, K4b) %>% pivot_longer(K4a:K4b, names_to = "nested_component", values_to = "Q") %>% mutate(global_component = "K4")
)
pS1c <- ggplot(hierarchy_long, aes(order, Q, fill = nested_component)) + geom_col(width = 1) + facet_grid(~global_component, scales = "free_x", space = "free_x") + scale_fill_manual(values = nested_colours) + theme_lampro +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = grid::unit(.04, "lines")) +
  labs(title = "Nested nuclear structure within global ancestry components", x = "Individuals ordered south to north", y = "Nested Q", fill = "Nested component")
save_figure(pS1c, "Supplementary_Figure_S1b_nested_snmf_ancestry", 12, 4.8)

#### Supplementary Figure S1c - SVDQuartets, Chapple labels and aligned Q ####

svd_tree <- read.tree(file.path(project, "treemix", "svdquartets_ref5", "svdquartets_ref5_170_random_final.tre"))
# PAUP writes bootstrap percentages in the raw Newick branch-value field. I
# redraw the final consensus as a cladogram because those values are support,
# not evolutionary branch lengths or divergence times.
svd_tree$edge.length <- rep(1, nrow(svd_tree$edge))
chapple_crosswalk <- sample_lineages %>% mutate(
  chapple_clade = ifelse(grepl("^deli[1-9]", chapple_lineage), sub("^deli([1-9]).*$", "deli\\1", chapple_lineage), chapple_lineage),
  chapple_clade = ifelse(chapple_clade == "unknown", "unknown", chapple_clade)
) %>% arrange(chapple_clade, chapple_lineage, id) %>%
  mutate(provenance = "Project population label; user-confirmed reference to Chapple et al. (2011) Figure 6 MCC tree")
write.csv(chapple_crosswalk, file.path(out, "chapple_figure6_sample_crosswalk.csv"), row.names = FALSE)

# Every terminal keeps its sample identifier and the project population label
# that refers to the Chapple Figure 6 MCC tree. The two horizontal Q panels use
# the tree's y-order, so each bar is directly aligned with its corresponding tip.
svd_metadata <- q4 %>% transmute(id = as.character(id), K1, K2, K3, K4) %>%
  left_join(chapple_crosswalk[, c("id", "chapple_lineage", "chapple_clade")], by = "id") %>%
  left_join(hierarchy_q %>% select(id, nested_component, nested_max_Q, K1a, K2a, K3a, K3b, K3c, K4a, K4b), by = "id") %>%
  mutate(label = id, tip_label = ifelse(is.na(chapple_lineage), paste0(id, " | no project population label"), paste0(id, " | ", chapple_lineage)))
global_q_long <- svd_metadata %>% select(label, K1:K4) %>% pivot_longer(K1:K4, names_to = "global_component", values_to = "Q")
nested_q_long <- svd_metadata %>% select(label, K1a, K2a, K3a, K3b, K3c, K4a, K4b) %>%
  pivot_longer(K1a:K4b, names_to = "nested_component", values_to = "Q") %>% filter(!is.na(Q))
svd_tree_data <- data.frame(label = svd_tree$tip.label, stringsAsFactors = FALSE) %>% left_join(svd_metadata, by = "label")
chapple_clade_colours <- setNames(grDevices::hcl.colors(length(unique(na.omit(svd_tree_data$chapple_clade))), "Dynamic"), sort(unique(na.omit(svd_tree_data$chapple_clade))))
pS1d <- ggtree(svd_tree, size = .22) %<+% svd_tree_data +
  geom_tippoint(aes(colour = chapple_clade), size = .45, na.rm = TRUE) +
  geom_tiplab(aes(label = tip_label), size = 1.18, align = TRUE, linesize = .14, na.rm = TRUE) +
  scale_colour_manual(values = chapple_clade_colours, na.value = "grey75", name = "Chapple Figure 6 clade") + theme_tree2()
tip_y <- pS1d$data %>% filter(isTip) %>% select(label, y)
global_q_long <- global_q_long %>% left_join(tip_y, by = "label")
nested_q_long <- nested_q_long %>% left_join(tip_y, by = "label")
tip_limits <- range(tip_y$y) + c(-.5, .5)
pS1d_global <- ggplot(global_q_long, aes(Q, y, fill = global_component)) +
  geom_col(orientation = "y", position = "stack", width = .82) + scale_fill_manual(values = ancestry_colours, name = "Global component") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = tip_limits, expand = c(0, 0)) +
  theme_void() + theme(plot.title = element_text(face = "bold", hjust = .5), legend.position = "bottom") + labs(title = "Global sNMF Q")
pS1d_nested <- ggplot(nested_q_long, aes(Q, y, fill = nested_component)) +
  geom_col(orientation = "y", position = "stack", width = .82) + scale_fill_manual(values = nested_colours, name = "Nested component") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = tip_limits, expand = c(0, 0)) +
  theme_void() + theme(plot.title = element_text(face = "bold", hjust = .5), legend.position = "bottom") + labs(title = "Nested sNMF Q")
pS1d <- (pS1d | pS1d_global | pS1d_nested) + plot_layout(widths = c(3.6, 1, 1)) +
  plot_annotation(title = "Final SVDQuartets consensus (100,000 random quartets; 1,000 bootstraps) aligned with Chapple Figure 6 labels and nuclear ancestry")
save_svd_figure(pS1d, "Supplementary_Figure_S1c_SVDquartets_global_nested_Q", 21, 22)

#### Supplementary diagnostics ####

ld_summary$distance_bin <- factor(ld_summary$distance_bin, levels = c("0-10 kb", "10-25 kb", "25-50 kb", "50-100 kb", "100-250 kb", ">250 kb"))
ld_long <- ld_summary %>% filter(!is.na(median_r2)) %>% select(distance_bin, median_r2, q95, prop_r2_gt_0_1) %>% pivot_longer(-distance_bin, names_to = "metric", values_to = "value")
ld_long$metric <- factor(ld_long$metric, levels = c("median_r2", "q95", "prop_r2_gt_0_1"), labels = c("Median r²", "95th percentile r²", "Proportion r² > 0.1"))
pS2 <- ggplot(ld_long, aes(distance_bin, value, group = 1)) + geom_line(colour = "#577590") + geom_point(colour = "#577590", size = 2) + facet_wrap(~metric, scales = "free_y", ncol = 1) + theme_lampro +
  labs(title = "Same-scaffold genotype LD by physical distance", x = "Distance bin", y = NULL)
save_figure(pS2, "Supplementary_Figure_S2_LD_audit", 7, 7)

climate_stability$label <- paste0(climate_stability$n_folds_selected, "/", climate_stability$n_folds_available)
pS3 <- ggplot(climate_stability, aes(reorder(variable, n_folds_selected), n_folds_selected)) + geom_col(fill = "#3F6F58") + geom_text(aes(label = label), hjust = -0.2, size = 3) + coord_flip() +
  scale_y_continuous(breaks = 0:5, limits = c(0, 5.5)) + theme_lampro + labs(title = "Climate selection stability", x = NULL, y = "Spatial training folds selecting predictor")
save_figure(pS3, "Supplementary_Figure_S3_climate_stability", 7, 5)

pS4 <- ggplot(leave_one_fold, aes(factor(omitted_fold), mean_Delta_RMSE)) + geom_hline(yintercept = 0, linetype = 2) + geom_hline(yintercept = .002739, colour = "#9A4A43", linetype = 3, linewidth = .8) + geom_col(fill = "#4C6A5D") + theme_lampro +
  labs(title = "Historical effect after omitting each fold", subtitle = "Red reference: complete-data equal-fold ΔRMSE = 0.002739", x = "Omitted fold", y = "Mean ΔRMSE")
save_figure(pS4, "Supplementary_Figure_S4_historical_leave_one_fold_out", 7, 4.5)

convergence$model <- factor(convergence$model, levels = c("SI", "IM", "AM", "SC"))
pS5 <- ggplot(convergence, aes(model, number_within_0.1, fill = model)) + geom_col(show.legend = FALSE) + facet_wrap(~pair) +
  scale_fill_manual(values = c(SI = "#7E7E7E", IM = "#577590", AM = "#B66B3E", SC = "#4C6A5D")) + theme_lampro +
  labs(title = "Stage 8B optimization confirmations", subtitle = "Independent optimisation solutions within 0.1 log-likelihood units of the frozen maximum", x = NULL, y = "Number of independent solutions")
save_figure(pS5, "Supplementary_Figure_S5_moments_convergence", 9, 4.5)

#### Graphical abstract ####

# This uses only empirical spatial coordinates, ancestry coefficients, selected
# nested assignments and locality-level historical gain. It represents neither
# inferred refugia nor literal TreeMix migration routes.
global_map <- ggplot() + geom_polygon(data = australia, aes(long, lat, group = group), fill = "#EAE4D9", colour = "grey65", linewidth = .25) +
  geom_point(data = q4_map, aes(X, Y, fill = component, alpha = ancestry), shape = 21, colour = "#EAE4D9", stroke = .1, size = 2.2) + coord_quickmap(xlim = c(135, 154), ylim = c(-44, -18)) +
  scale_fill_manual(values = ancestry_colours, guide = "none") + scale_alpha(range = c(.08, 1), guide = "none") + theme_void() + ggtitle("BROAD ANCESTRY TURNOVER")
past_map <- ggplot() + geom_polygon(data = australia, aes(long, lat, group = group), fill = "#EAE4D9", colour = "grey65", linewidth = .25) +
  geom_point(data = local_gain, aes(X, Y, colour = mean_pairwise_gain_squared_error), size = 2.4) + coord_quickmap(xlim = c(135, 154), ylim = c(-44, -18)) +
  scale_colour_gradient2(low = history_colours[["negative"]], mid = history_colours[["neutral"]], high = history_colours[["positive"]], midpoint = 0, guide = "none") + theme_void() + ggtitle("PAST LANDSCAPE")
nested_plot <- pS1c + theme(legend.position = "none", plot.title = element_text(size = 9), axis.title = element_text(size = 8)) + ggtitle("NESTED NUCLEAR STRUCTURE")
reticulate_plot <- p4a + theme(legend.position = "none", axis.text = element_text(size = 6), axis.title = element_text(size = 7), plot.title = element_text(size = 9), plot.subtitle = element_blank()) + ggtitle("RETICULATION / PRESENT GENOMIC PATTERN")
graphical_abstract <- (past_map | global_map | nested_plot | reticulate_plot) + plot_annotation(title = "Historical landscape signal and hierarchical genomic structure in Lampropholis delicata")
save_figure(graphical_abstract, "Graphical_Abstract", 16, 4.5)

#### Record figure outputs and input checks ####

manifest <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Supplementary Figure S1", "Supplementary Figure S1b", "Supplementary Figure S1c", "Supplementary Figure S2", "Supplementary Figure S3", "Supplementary Figure S4", "Supplementary Figure S5", "Graphical abstract"),
  panel = c("A-B", "A-D", "A-E", "A-B", "A-D", "A-B", "A", "A", "A", "A", "A", "A", "A"),
  source_file = c("pcoa_scores.csv; pcoa.ref7.rds; snmf_selected_qmatrix_k4.csv", "Stage 5 GDM tables and validation", "Stage 6 validation, corrected null and spatial-effect tables", "TreeMix likelihood replicates and m4 graph", "Stage 8B model table and folded-SFS arrays", "hierarchical_snmf runs and model selection", "hierarchical_snmf assignments", "SVDQuartets tree and hierarchical assignments", "LD distance summary", "climate selection stability", "leave-one-fold diagnostic", "Stage 8B convergence table", "empirical spatial, ancestry, nested-Q and historical-gain outputs"),
  status = c("MAIN", "MAIN", "MAIN", "MAIN", "MAIN_PROVISIONAL", rep("SUPPLEMENTARY", 7), "GRAPHICAL_ABSTRACT"),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out, "FIGURE_MANIFEST.csv"), row.names = FALSE)
audit_lines <- c("# Figure audit", "", "I generated the versioned figures from frozen results only. I did not rerun sNMF, GDM, TreeMix, moments, null models or demographic fitting.", "", "- Figure 1 now uses the archived dartR PCoA eigenvalues (28.2% and 10.8%) and K4 map opacity represents mean Q.", "- Figure 2 distinguishes descriptive full-data fit from spatially held-out results and labels observed spatial-null statistics.", "- Figure 3 retains equal-fold inference, labels the pooled pair-weighted locality diagnostic, and places Australian coastline beneath the spatial panels.", "- Figure 4 displays all ten TreeMix replicate fits for each migration-edge number; m = 4 remains a parsimonious working graph.", "- Figure 5 uses expected folded spectra rescaled only for visualization. The likelihood calculations were not changed; Stage 8C remains pending.", "- The SVDQuartets supplementary tree now uses the project population labels documented by the user as references to Chapple et al. (2011) Figure 6 MCC lineages, alongside the nested nuclear Q column. The crosswalk is saved with the figure outputs.")
writeLines(audit_lines, file.path(out, "FIGURE_AUDIT.md"))
