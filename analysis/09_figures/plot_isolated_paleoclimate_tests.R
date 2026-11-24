# Separate effects of Last Interglacial and Last Glacial Maximum conditions
# Plot the historical climate comparisons by period and visualize their performance under spatial validation. The period-specific summaries complement the main combined historical analysis.
# Check the comparison baseline and the relevant conditional permutation result before interpreting an apparent improvement.
#
#### Packages ####

# ggplot2 produces the compact manuscript figure from the previously estimated Stage 6B outputs.
require(ggplot2)
# dplyr keeps the period-specific summaries transparent and readable.
require(dplyr)
# gridExtra places the three scientific comparisons on one publication plate.
require(gridExtra)

#### Directories ####

setwd(Sys.getenv("LAMPROPHOLIS_PROJECT", unset = "."))

dir.create("figures/manuscript", recursive = TRUE, showWarnings = FALSE)

#### Load isolated palaeoclimate tests ####

observed <- read.csv("results/landscape/paleoclimate/null_period_specific/isolated_period_observed_fold_effects.csv")
null <- read.csv("results/landscape/paleoclimate/null_period_specific/isolated_period_null_statistics.csv")
summary <- read.csv("results/landscape/paleoclimate/null_period_specific/isolated_period_null_summary.csv")

period.labels <- c(LGM = "LGM projected climate", LIG = "LIG projected climate")
period.colours <- c(LGM = "#4C78A8", LIG = "#C65D7B")

observed$period <- factor(observed$period, levels = c("LGM", "LIG"))
null$period <- factor(null$period, levels = c("LGM", "LIG"))
summary$period <- factor(summary$period, levels = c("LGM", "LIG"))

#### Build figure ####

fold.figure <- ggplot(observed, aes(factor(fold), Delta_RMSE, fill = period)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = period.colours, labels = period.labels) +
  labs(x = "Spatial fold", y = "Delta RMSE (M2 - M2 + period)", fill = NULL,
       title = "Fold-specific held-out gain") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        legend.position = "top", legend.text = element_text(size = 7),
        axis.title = element_text(size = 8))

null.figure <- ggplot(null, aes(Delta_RMSE, fill = period)) +
  geom_histogram(bins = 34, alpha = 0.72, position = "identity", colour = "white", linewidth = 0.18) +
  geom_vline(data = summary, aes(xintercept = observed_mean_Delta_RMSE, colour = period), linewidth = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = period.colours, labels = period.labels) +
  scale_colour_manual(values = period.colours) +
  facet_wrap(~period, ncol = 1, scales = "free_y", labeller = as_labeller(period.labels)) +
  labs(x = "Mean Delta RMSE under locality/block null", y = "Permutations",
       title = "Observed gain against spatial null") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5), legend.position = "none",
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        axis.title = element_text(size = 8))

summary$label <- paste0("p = ", format.pval(summary$p_value, digits = 2, eps = 0.001),
                        "\n", summary$folds_improved, "/5 folds improved")

summary.figure <- ggplot(summary, aes(period, observed_mean_Delta_RMSE, fill = period)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = label), vjust = -0.25, size = 2.5, lineheight = 0.92) +
  scale_fill_manual(values = period.colours, guide = "none") +
  scale_x_discrete(labels = period.labels) +
  expand_limits(y = max(summary$observed_mean_Delta_RMSE) * 1.6) +
  labs(x = NULL, y = "Mean held-out Delta RMSE", title = "Isolated period effects") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 7), axis.title = element_text(size = 8))

#### Save figure ####

figure.plate <- arrangeGrob(summary.figure, fold.figure, null.figure, ncol = 3, widths = c(0.95, 1, 1.25))

png("figures/manuscript/Supplementary_Figure_S7_isolated_LIG_LGM_tests.png", width = 3300, height = 1200, res = 300)
grid::grid.draw(figure.plate)
dev.off()

pdf("figures/manuscript/Supplementary_Figure_S7_isolated_LIG_LGM_tests.pdf", width = 11.5, height = 4.4)
grid::grid.draw(figure.plate)
dev.off()
