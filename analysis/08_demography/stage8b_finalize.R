#### Directories ####

project_dir <- normalizePath(".")
table_dir <- file.path(project_dir, "results", "demography", "moments", "tables")

#### Load validated optimization results ####

pilot <- read.csv(file.path(table_dir, "stage8_clean_pilot_runs.csv"), stringsAsFactors = FALSE)
pilot <- pilot[pilot$optimization_status == "success", ]
pilot$optimizer <- "PILOT_FMIN"
pilot$logLik <- as.numeric(pilot$log_likelihood)
pilot$source_stage <- "Stage8_clean_pilot"
pilot <- pilot[, c("pair", "model", "optimizer", "optimized_parameters", "logLik", "source_stage")]

stage8a <- read.csv(file.path(table_dir, "stage8a_cross_optimizer_runs.csv"), stringsAsFactors = FALSE)
stage8a <- stage8a[stage8a$optimizer_status == "success", ]
stage8a$logLik <- as.numeric(stage8a$final_logLik)
stage8a$source_stage <- "Stage8A"
stage8a <- stage8a[, c("pair", "model", "optimizer", "optimized_parameters", "logLik", "source_stage")]

stage8b <- read.csv(file.path(table_dir, "stage8b_confirmatory_runs.csv"), stringsAsFactors = FALSE)
stage8b <- stage8b[stage8b$optimizer_status == "success", ]
stage8b$key <- paste(stage8b$pair, stage8b$model, stage8b$replicate, stage8b$optimizer, sep = "::")
stage8b <- stage8b[!duplicated(stage8b$key, fromLast = TRUE), ]
stage8b$logLik <- as.numeric(stage8b$logLik)
stage8b$source_stage <- "Stage8B"
stage8b <- stage8b[, c("pair", "model", "optimizer", "optimized_parameters", "logLik", "source_stage")]

all_runs <- rbind(pilot, stage8a, stage8b)

#### Freeze reproducible maxima ####

groups <- split(all_runs, interaction(all_runs$pair, all_runs$model, drop = TRUE))
maxima <- list(); convergence <- list(); limiting <- list(); parameters <- list(); comparison <- list()
for (key in names(groups)) {
  x <- groups[[key]]
  x <- x[order(-x$logLik), ]
  best <- x[1, ]
  within01 <- x[x$logLik >= best$logLik - 0.1, ]
  within1 <- x[x$logLik >= best$logLik - 1, ]
  within2 <- x[x$logLik >= best$logLik - 2, ]
  params <- as.numeric(strsplit(best$optimized_parameters, ";", fixed = TRUE)[[1]])
  lower <- switch(best$model, SI=c(.01,.01,1e-4), IM=c(.01,.01,1e-4,1e-5), AM=c(.01,.01,1e-4,1e-4,1e-5), SC=c(.01,.01,1e-4,1e-4,1e-5))
  upper <- switch(best$model, SI=c(100,100,10), IM=c(100,100,10,50), AM=c(100,100,10,10,50), SC=c(100,100,10,10,50))
  position <- (log(params)-log(lower))/(log(upper)-log(lower))
  bound_status <- if(any(position <= .02 | position >= .98)) "near or at log-scale bound" else "no log-scale boundary support"
  k <- switch(best$model, SI = 3, IM = 4, AM = 5, SC = 5)
  maxima[[length(maxima) + 1]] <- data.frame(pair=best$pair, model=best$model, final_best_logLik=best$logLik, final_best_parameters=best$optimized_parameters, source_stage=best$source_stage, source_optimizer=best$optimizer, reproducible_n=nrow(within01))
  convergence[[length(convergence) + 1]] <- data.frame(pair=best$pair, model=best$model, final_best_logLik=best$logLik, number_independent_confirmations=nrow(within01), number_within_0.1=nrow(within01), number_within_1=nrow(within1), likelihood_convergence=if(nrow(within01)>=2) "CONVERGED" else "SINGLE_SOLUTION", parameter_convergence="see parameter table", boundary_status=bound_status, limiting_model_status="see limiting audit")
  if (best$model == "SI") {
    limiting_status <- "not applicable"
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter=c("nu1","nu2","T"), value=params)
  } else if (best$model == "IM") {
    limiting_status <- if(params[4] <= 1e-4) "approaches SI" else "no limiting collapse"
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter=c("nu1","nu2","T","M"), value=params)
  } else if (best$model == "AM") {
    limiting_status <- if(params[5] <= 1e-4 || params[3] <= 1e-4) "approaches simpler history" else "no limiting collapse"
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter=c("nu1","nu2","T_mig","T_iso","M"), value=params)
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter="T_total", value=params[3]+params[4])
  } else {
    limiting_status <- if(params[5] <= 1e-4 || params[4] <= 1e-4) "approaches simpler history" else "no limiting collapse"
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter=c("nu1","nu2","T_iso","T_sc","M"), value=params)
    parameters[[length(parameters)+1]] <- data.frame(pair=best$pair, model=best$model, parameter="T_total", value=params[3]+params[4])
  }
  limiting[[length(limiting)+1]] <- data.frame(pair=best$pair, model=best$model, limiting_model_status=limiting_status, boundary_status=bound_status)
  comparison[[length(comparison)+1]] <- data.frame(pair=best$pair, model=best$model, k=k, logLik=best$logLik, AIC=2*k-2*best$logLik)
}

maxima <- do.call(rbind, maxima); convergence <- do.call(rbind, convergence); limiting <- do.call(rbind, limiting); parameters <- do.call(rbind, parameters); comparison <- do.call(rbind, comparison)
comparison <- do.call(rbind, lapply(split(comparison, comparison$pair), function(x) { x$delta_AIC <- x$AIC-min(x$AIC); w <- exp(-.5*x$delta_AIC); x$Akaike_weight <- w/sum(w); x$rank <- rank(x$AIC, ties.method="first"); x }))

write.csv(maxima, file.path(table_dir, "stage8b_final_model_maxima.csv"), row.names=FALSE)
write.csv(convergence, file.path(table_dir, "stage8b_final_convergence.csv"), row.names=FALSE)
write.csv(comparison, file.path(table_dir, "stage8b_model_comparison.csv"), row.names=FALSE)
write.csv(limiting, file.path(table_dir, "stage8b_limiting_model_audit.csv"), row.names=FALSE)
write.csv(parameters, file.path(table_dir, "stage8b_relative_parameter_summary.csv"), row.names=FALSE)

#### K1-K3 SC identifiability ####

sc <- all_runs[all_runs$pair == "K1_K3" & all_runs$model == "SC", ]
sc <- sc[sc$logLik >= max(sc$logLik)-0.1, ]
mat <- do.call(rbind, lapply(sc$optimized_parameters, function(z) as.numeric(strsplit(z, ";", fixed=TRUE)[[1]])))
sc_summary <- data.frame(parameter=c("nu1","nu2","T_iso","T_sc","M"), minimum=apply(mat,2,min), median=apply(mat,2,median), maximum=apply(mat,2,max), log_range=apply(log(mat),2,function(z) diff(range(z))))
sc_total <- data.frame(parameter="T_total", minimum=min(mat[,3]+mat[,4]), median=median(mat[,3]+mat[,4]), maximum=max(mat[,3]+mat[,4]), log_range=diff(range(log(mat[,3]+mat[,4]))))
write.csv(rbind(sc_summary, sc_total), file.path(table_dir, "stage8b_k1_k3_sc_identifiability.csv"), row.names=FALSE)
