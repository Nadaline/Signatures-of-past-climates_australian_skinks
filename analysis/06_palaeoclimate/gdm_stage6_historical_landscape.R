#!/usr/bin/env Rscript
# Historical landscape predictors in generalized dissimilarity models
# Project the climate transformations fitted using present-day observations onto Last Interglacial and Last Glacial Maximum environmental layers. Evaluate whether these transformed historical differences improve spatially withheld predictions.
# Check climate units before calculating environmental distances. Historical effects are tested beyond geography, contemporary climate and topographic resistance on comparable withheld pairs.
#

#### Packages ####

suppressPackageStartupMessages({library(dartR.base); library(gdm); library(terra); library(raster); library(ggplot2)})

#### Directories ####

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
seed <- if (length(args) >= 2) as.integer(args[2]) else 20260812L
n_perm_final <- if (length(args) >= 3) as.integer(args[3]) else 999L
set.seed(seed)
stage5 <- file.path(project, "results", "landscape", "gdm_definitive")
base <- file.path(project, "results", "landscape", "paleoclimate")
dirs <- setNames(file.path(base, c("tables", "validation", "figures", "models", "logs")), c("tables", "validation", "figures", "models", "logs"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

#### Read contemporary GDM estimates used as the reference model ####

climate_variables <- c("tempRange", "tempSeason", "precAnn")
gl <- readRDS(file.path(project, "data", "derived", "lampro_landscape_background.rds"))
ids <- indNames(gl)
resistance <- as.matrix(read.csv(file.path(project, "results", "connectivity", "source_matrices", "Lampro_resistances_pairwise.csv"), row.names = 1, check.names = FALSE))
full_model <- readRDS(file.path(stage5, "models", "gcr_seven_models.rds"))[["CGR"]]

#### Load Present climate ####

present <- rast(file.path(Sys.getenv("LAMPRO_CHELSA_DIR","data/local/chelsa"), paste0("CHELSA_", c("bio2", "bio4", "bio12"), "_1981-2010_australia.tif")))
names(present) <- climate_variables

#### Load LIG ####

lig_all <- readRDS(file.path(project, "data", "paleoclimate", "paleoclim_lig_5m_australia.rds"))
lig <- rast(lig_all[[c("bio_2", "bio_4", "bio_12")]])
names(lig) <- climate_variables

#### Load LGM ####

lgm <- rast(readRDS(file.path(project, "data", "paleoclimate", "paleoclim_lgm_5m_australia.rds")))
names(lgm) <- climate_variables

#### Harmonize climate layers ####

# PaleoClim BIO2 and BIO4 are expressed in tenths of the current CHELSA units;
# BIO12 is already in millimetres. I apply only this unit conversion.
lig[[c("tempRange", "tempSeason")]] <- lig[[c("tempRange", "tempSeason")]] / 10
lgm[[c("tempRange", "tempSeason")]] <- lgm[[c("tempRange", "tempSeason")]] / 10
writeRaster(lig, file.path(project, "data", "paleoclimate", "lig_bio2_bio4_bio12_harmonized.tif"), overwrite = TRUE)
writeRaster(lgm, file.path(project, "data", "paleoclimate", "lgm_bio2_bio4_bio12_harmonized.tif"), overwrite = TRUE)

#### Check sampled localities ####

coords <- as.data.frame(gl@other$latlon)
names(coords)[match(c("lon", "lat"), names(coords))] <- c("X", "Y")
coords$id <- ids
extract_climate <- function(x) { z <- terra::extract(x, coords[, c("X", "Y")])[, -1, drop = FALSE]; colnames(z) <- climate_variables; z }
present_values <- extract_climate(present)
lig_values <- extract_climate(lig)
lgm_values <- extract_climate(lgm)
common <- ids[ids %in% rownames(resistance) & ids %in% colnames(resistance) & complete.cases(present_values, lig_values, lgm_values)]
gl <- gl[match(common, indNames(gl)), ]
coords <- coords[match(common, coords$id), ]
present_values <- present_values[match(common, ids), , drop = FALSE]
lig_values <- lig_values[match(common, ids), , drop = FALSE]
lgm_values <- lgm_values[match(common, ids), , drop = FALSE]
resistance <- resistance[common, common]
genomic <- as.matrix(gl.dist.ind(gl, scale = TRUE))[common, common]
env <- cbind(coords, present_values)
env$site <- seq_len(nrow(env))
folds <- read.csv(file.path(stage5, "validation", "spatial_folds.csv"))
env$fold <- folds$fold[match(env$id, folds$id)]
stopifnot(nrow(env) == 170, !anyNA(env$fold), isTRUE(all.equal(resistance, t(resistance))))

#### GDM helpers ####

matrix_data <- function(x, index) { z <- x[index, index, drop = FALSE]; out <- data.frame(site = seq_along(index), z, check.names = FALSE); names(out)[-1] <- seq_along(index); out }
sitepair <- function(index, historical = NULL) {
  e <- env[index, c("site", "X", "Y", climate_variables)]
  input <- list(bioData = matrix_data(genomic, index), bioFormat = 3, XColumn = "X", YColumn = "Y", siteColumn = "site", predData = e, distPreds = list(matrix_data(resistance, index)))
  if (!is.null(historical)) input$distPreds <- c(input$distPreds, list(matrix_data(historical, index)))
  z <- do.call(formatsitepair, input)
  names(z) <- sub("matrix_1", "topographicResistance", names(z), fixed = TRUE)
  names(z) <- sub("matrix_2", "historicalClimate", names(z), fixed = TRUE)
  z
}
fit_model <- function(index, historical = NULL) gdm(sitepair(index, historical), geo = TRUE)
performance <- function(model, table) { predicted <- as.numeric(predict(model, table)); observed <- table$distance; data.frame(n_pairs = nrow(table), rmse = sqrt(mean((observed - predicted)^2)), mae = mean(abs(observed - predicted)), correlation = if (sd(predicted) > 0) cor(observed, predicted) else NA_real_) }

#### Transform climate with GDM splines ####

# I interpolate the climate I-splines fitted by M2. This uses the contemporary
# genomic transformation but does not refit a genomic response to palaeoclimate.
transformed_distance <- function(model, values) {
  sp <- isplineExtract(model)
  out <- matrix(0, nrow(values), nrow(values))
  for (v in climate_variables) {
    d <- abs(outer(values[, v], values[, v], "-"))
    out <- out + matrix(approx(sp$x[, v], sp$y[, v], xout = as.vector(d), rule = 2)$y, nrow(d), ncol(d))
  }
  diag(out) <- 0
  out
}

#### Pairwise palaeoclimate distances ####

c_present <- transformed_distance(full_model, present_values)
c_lig <- transformed_distance(full_model, lig_values)
c_lgm <- transformed_distance(full_model, lgm_values)
paleo_mean <- (c_lig + c_lgm) / 2
paleo_change <- abs(c_lig - c_lgm)
pair <- which(lower.tri(c_present), arr.ind = TRUE)
pair_table <- data.frame(id1 = common[pair[, 1]], id2 = common[pair[, 2]], C_present = c_present[pair], C_LIG = c_lig[pair], C_LGM = c_lgm[pair], PaleoMean = paleo_mean[pair], PaleoChange = paleo_change[pair])
write.csv(pair_table, file.path(dirs["tables"], "palaeo_pairwise_climate_distances.csv"), row.names = FALSE)
write.csv(cor(pair_table[, 3:7], method = "spearman"), file.path(dirs["tables"], "historical_predictor_correlations.csv"), row.names = TRUE)
write.csv(data.frame(metric = c("PaleoMean", "PaleoChange"), formula = c("(C_LIG + C_LGM) / 2", "abs(C_LIG - C_LGM)")), file.path(dirs["tables"], "historical_predictor_definition.csv"), row.names = FALSE)

#### Spatial cross-validation ####

cv_rows <- list(); gain_rows <- list()
for (f in sort(unique(env$fold))) {
  train <- which(env$fold != f); test <- which(env$fold == f)
  m2 <- fit_model(train)
  lig_distance <- transformed_distance(m2, lig_values)
  lgm_distance <- transformed_distance(m2, lgm_values)
  mean_distance <- (lig_distance + lgm_distance) / 2
  models <- list(M2 = m2, M2_LGM = fit_model(train, lgm_distance), M2_LIG = fit_model(train, lig_distance), M3_PaleoMean = fit_model(train, mean_distance))
  predictors <- list(M2 = NULL, M2_LGM = lgm_distance, M2_LIG = lig_distance, M3_PaleoMean = mean_distance)
  for (nm in names(models)) cv_rows[[paste(f, nm)]] <- cbind(fold = f, model = nm, n_train = length(train), n_test = length(test), train_deviance = models[[nm]]$explained, performance(models[[nm]], sitepair(test, predictors[[nm]])))
  base_rmse <- cv_rows[[paste(f, "M2")]]$rmse
  gain_rows[[as.character(f)]] <- data.frame(fold = f, delta_LGM_rmse = base_rmse - cv_rows[[paste(f, "M2_LGM")]]$rmse, delta_LIG_rmse = base_rmse - cv_rows[[paste(f, "M2_LIG")]]$rmse, delta_PaleoMean_rmse = base_rmse - cv_rows[[paste(f, "M3_PaleoMean")]]$rmse)
}
cv <- do.call(rbind, cv_rows); gains <- do.call(rbind, gain_rows)
write.csv(cv, file.path(dirs["validation"], "historical_model_comparison.csv"), row.names = FALSE)
write.csv(gains, file.path(dirs["validation"], "historical_cv_incremental_effects.csv"), row.names = FALSE)

#### Historical null models ####

permute_localities <- function(index, values) {
  out <- values
  locality <- paste(round(env$X, 5), round(env$Y, 5), sep = "_")
  for (b in unique(env$fold[index])) {
    sites <- unique(locality[index][env$fold[index] == b])
    if (length(sites) > 1) { shuffled <- setNames(sample(sites), sites); destination <- which(env$fold == b); source <- match(unname(shuffled[locality[destination]]), locality); out[destination, ] <- values[source, ] }
  }
  out
}
null_rows <- list()
for (f in sort(unique(env$fold))) {
  train <- which(env$fold != f); test <- which(env$fold == f)
  m2 <- fit_model(train)
  lig_distance <- transformed_distance(m2, lig_values); lgm_distance <- transformed_distance(m2, lgm_values); mean_distance <- (lig_distance + lgm_distance) / 2
  m3 <- fit_model(train, mean_distance)
  observed <- performance(m2, sitepair(test))$rmse - performance(m3, sitepair(test, mean_distance))$rmse
  null_rows[[paste(f, 0)]] <- data.frame(fold = f, permutation = 0, paleo_gain_rmse = observed)
  for (k in seq_len(n_perm_final)) {
    lig_null_values <- permute_localities(train, lig_values)
    lgm_null_values <- permute_localities(train, lgm_values)
    mean_null <- (transformed_distance(m2, lig_null_values) + transformed_distance(m2, lgm_null_values)) / 2
    null_model <- fit_model(train, mean_null)
    null_rows[[paste(f, k)]] <- data.frame(fold = f, permutation = k, paleo_gain_rmse = performance(m2, sitepair(test))$rmse - performance(null_model, sitepair(test, mean_distance))$rmse)
  }
}
nulls <- do.call(rbind, null_rows)
write.csv(nulls, file.path(dirs["validation"], "historical_null_models.csv"), row.names = FALSE)
observed <- mean(nulls$paleo_gain_rmse[nulls$permutation == 0])
null_mean <- aggregate(paleo_gain_rmse ~ permutation, nulls[nulls$permutation > 0, ], mean)
null_summary <- data.frame(metric = "PaleoMean", observed_mean_rmse_gain = observed, null_mean = mean(null_mean$paleo_gain_rmse), null_median = median(null_mean$paleo_gain_rmse), null_q025 = quantile(null_mean$paleo_gain_rmse, .025), null_q975 = quantile(null_mean$paleo_gain_rmse, .975), p_value = (1 + sum(null_mean$paleo_gain_rmse >= observed)) / (n_perm_final + 1), nperm = n_perm_final)
write.csv(null_summary, file.path(dirs["validation"], "historical_null_summary.csv"), row.names = FALSE)

#### Figures ####

cv_plot <- ggplot(cv, aes(model, rmse, group = fold)) + geom_hline(yintercept = 0, colour = "white") + geom_line(colour = "grey70") + geom_point(aes(colour = model), size = 2) + theme_classic() + guides(colour = "none") + labs(x = NULL, y = "Held-out RMSE")
ggsave(file.path(dirs["figures"], "historical_spatial_cv.png"), cv_plot, width = 7, height = 5, dpi = 240)
gain_long <- reshape(gains, varying = names(gains)[-1], v.names = "gain", timevar = "metric", times = c("LGM", "LIG", "PaleoMean"), direction = "long")
gain_plot <- ggplot(gain_long, aes(metric, gain, group = fold)) + geom_hline(yintercept = 0, linetype = 2) + geom_line(colour = "grey70") + geom_point(size = 2) + theme_classic() + labs(x = NULL, y = "Held-out RMSE gain over M2")
ggsave(file.path(dirs["figures"], "historical_fold_gains.png"), gain_plot, width = 7, height = 5, dpi = 240)
null_plot <- ggplot(null_mean, aes(paleo_gain_rmse)) + geom_histogram(bins = 35, fill = "grey75", colour = "white") + geom_vline(xintercept = observed, colour = "#8F3B3F", linewidth = .9) + theme_classic() + labs(x = "Mean held-out RMSE gain", y = "Spatial-null permutations")
ggsave(file.path(dirs["figures"], "historical_null.png"), null_plot, width = 7, height = 5, dpi = 240)

#### Save results ####

write.csv(data.frame(seed = seed, n_perm_final = n_perm_final, primary_metric = "PaleoMean", final_climate_variables = paste(climate_variables, collapse = ";"), genomic_input = file.path(project, "data", "derived", "lampro_landscape_background.rds")), file.path(dirs["logs"], "parameters.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(dirs["logs"], "sessionInfo.txt"))
