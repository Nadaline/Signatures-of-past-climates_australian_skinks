#### Recolour existing plots; never fit models ####
arg <- commandArgs(TRUE)
kind <- arg[1]
stage <- normalizePath(arg[2])
local <- normalizePath(file.path(stage,"../../../.."))
project <- arg[3]
base <- file.path(local,"results")
boundary <- file.path(base,"Lampropholis_Fig1_Geographic_context/data/abs_states_generalized.geojson")
figures <- file.path(stage,"figures",kind)
dir.create(figures,showWarnings=FALSE)
p <- read.csv(file.path(stage,"palette/chapter2_palette.csv"))
c <- setNames(p$hex,p$key)
macro <- setNames(unname(c[c("violet","ink","lavender","amber")]),paste0("K",1:4))
nested <- c(macro[1:2],c[c("K3a","K3b","K3c","K4a","K4b")])
env <- new.env()
if (kind=="structure") {
  script <- file.path(base,"Fig2_Nuclear_structure_revision/Fig2_Nuclear_structure.R")
  env$commandArgs <- function(...) c(project,boundary,figures)
  z <- read.csv(file.path(dirname(script),"P4_palette.csv"))
  colors <- c(macro,nested[-(1:2)],K1a=macro[[1]],K2a=macro[[2]])
  i <- z$role %in% names(colors)
  z$colour[i] <- unname(colors[z$role[i]])
  write.csv(z,file.path(figures,"P4_palette.csv"),row.names=FALSE)
} else if(kind %in% c("svd","svd_bootstrap","svd_discussion")) {
  script <- file.path(base,"Lampropholis_Fig3_SVDQuartets_geography_v4/Fig3_SVDQuartets_geography_v4.R")
  env$commandArgs <- function(...) c(project,boundary,figures,
    if(kind=="svd_bootstrap") "--bootstrap" else if(kind=="svd_discussion") "--discussion")
  # The bootstrap edition uses the original Figure 3 reference palette,
  # including its nested assignments, mixed symbols and K3 link shades.
  if (kind=="svd") {
    env$macro <- macro
    env$nested <- nested
    env$k3_link <- nested[c("K3a","K3b","K3c")]
  }
} else if(kind=="treemix") {
  script <- file.path(base,"Lampropholis_TreeMix_GDM_final_revision/TreeMix_revised.R")
  helpers <- Sys.glob("~/miniconda3/pkgs/treemix-1.13*/bin/plotting_funcs.R")
  stopifnot(length(helpers)>=1)
  tips <- file.path(base,"Lampropholis_Fig3_SVDQuartets_geography_v4/Fig3_sample_tip_key.csv")
  env$commandArgs <- function(...) c(project,figures,helpers[1],tips)
  env$pal <- macro
  env$line_col <- colorRampPalette(unname(c[c("lavender","violet","ink")]))(256)
} else if(kind=="historical") {
  script <- file.path(base,"Lampropholis_Figure5_palaeoclimate/Figure5_palaeoclimate_final.R")
  env$commandArgs <- function(...) c(project,boundary,figures)
  env$fold_col <- unname(c[c("violet","ink","K4a","lavender","amber")])
  env$pal <- colorRampPalette(unname(c[c("violet","paper","amber")]))(257)
}
for(expr in parse(script)) {
  lhs <- if(is.call(expr) && identical(expr[[1]],as.name("<-"))) deparse(expr[[2]]) else ""
  if (!lhs %in% c(if(kind=="svd") c("macro","nested","k3_link") else
      if(kind=="treemix") c("pal","line_col") else if(kind=="historical") c("fold_col","pal"))) {
    eval(expr,env)
  }
}
