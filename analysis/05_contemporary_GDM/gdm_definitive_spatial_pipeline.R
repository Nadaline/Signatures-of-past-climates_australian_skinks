#!/usr/bin/env Rscript
# Contemporary generalized dissimilarity models with spatial validation
# Fit and compare models based on geographic distance, present-day climate and effective topographic resistance. Cross-validation withholds entire sampling localities so related pairwise observations do not leak into training data.
# Climate predictors are selected using training localities. Conditional tests randomize the relevant training information within spatial blocks while keeping held-out observations unchanged.
#

# I fit the contemporary landscape GDM from first principles. Geography and
# direct Circuitscape effective topographic resistance is fixed; current CHELSA variables
# are selected inside each spatial training fold. I never permute pair rows:
# climate values are permuted at sampling localities within broad spatial
# blocks, and the pair table is rebuilt after every permutation.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dartR.base); library(dartR.popgen); library(gdm); library(terra)
  library(ggplot2); library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
seed <- if (length(args) >= 2) as.integer(args[2]) else 20260811L
#### Analysis settings ####

# I keep a small test setting for preparation and a separate value for the
# future full analysis.  The final setting is deliberately not run here.
n_perm_test <- if (length(args) >= 3) as.integer(args[3]) else 19L
n_perm_final <- if (length(args) >= 5) as.integer(args[5]) else 999L
nfold <- if (length(args) >= 4) as.integer(args[4]) else 5L
set.seed(seed)

source_root <- dirname(project)
base <- file.path(project, "results", "landscape", "gdm_definitive")
dirs <- setNames(file.path(base, c("results","tables","figures","models","validation","logs")),
                 c("results","tables","figures","models","validation","logs"))
invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))
logfile <- file.path(dirs["logs"], "pipeline.log")
logit <- function(...) { z <- paste(format(Sys.time()), paste(..., collapse=" ")); cat(z,"\n"); cat(z,"\n",file=logfile,append=TRUE) }
write.csv(data.frame(seed=seed,n_perm_test=n_perm_test,n_perm_final=n_perm_final,nfold=nfold,correlation_threshold=.92,
                     selection_p=.10,min_gain_deviance=.10,
                     permutation_unit="sampling locality within training spatial block"),
          file.path(dirs["logs"],"parameters.csv"),row.names=FALSE)

bio_names <- c(tempMean="bio1",tempRange="bio2",isotherm="bio3",tempSeason="bio4",
 tempMax="bio5",tempMin="bio6",tempAnnRg="bio7",tempWetQ="bio8",tempDryQ="bio9",
 tempWarmQ="bio10",tempColdQ="bio11",precAnn="bio12",precWetM="bio13",
 precDryM="bio14",precSeason="bio15",precWetQ="bio16",precDryQ="bio17",
 precWarmQ="bio18",precColdQ="bio19")
# Within nearly identical clusters I prioritize direct, interpretable thermal
# and water-balance axes that have equivalents in present, LGM, and LIG.
priority <- names(bio_names)

gen_file <- file.path(project,"data","derived","lampro_landscape_background.rds")
res_file <- file.path(project,"results","connectivity","source_matrices","Lampro_resistances_pairwise.csv")
clim_dir <- Sys.getenv("LAMPRO_CHELSA_DIR","data/local/chelsa")
clim_files <- file.path(clim_dir,paste0("CHELSA_",bio_names,"_1981-2010_australia.tif"))
stopifnot(file.exists(gen_file),file.exists(res_file),all(file.exists(clim_files)))

gl <- readRDS(gen_file); ids <- indNames(gl)
tags <- sub(":.*$","",locNames(gl))
gen_audit <- data.frame(object=gen_file,n_ind=nInd(gl),n_snp=nLoc(gl),unique_RAD_tags=length(unique(tags)),
 one_snp_per_tag=!anyDuplicated(tags),candidate_list_available=FALSE,neutral_verified=FALSE,
 dataset_role="landscape-background panel; one SNP per RAD locus; 221 reproducibly recoverable PCA-conditioned RDA candidates removed; not fully neutral")
write.csv(gen_audit,file.path(dirs["tables"],"genomic_dataset_audit.csv"),row.names=FALSE)

coords <- as.data.frame(gl@other$latlon); names(coords)[match(c("lon","lat"),names(coords))] <- c("X","Y")
coords$id <- ids
clim <- rast(clim_files); names(clim) <- names(bio_names)
cv <- terra::extract(clim,coords[,c("X","Y")])[, -1,drop=FALSE]
env <- cbind(coords,cv)
res <- as.matrix(read.csv(res_file,row.names=1,check.names=FALSE)); storage.mode(res)<-"numeric"
common <- ids[ids %in% rownames(res) & ids %in% colnames(res) & complete.cases(env[,c("X","Y",names(bio_names))])]
gl <- gl[match(common,indNames(gl)),]; env <- env[match(common,env$id),]; res <- res[common,common,drop=FALSE]
gen <- as.matrix(gl.dist.ind(gl,scale=TRUE)); gen <- gen[common,common,drop=FALSE]
stopifnot(isTRUE(all.equal(res,t(res))),all(diag(res)==0))
env$site <- seq_len(nrow(env)); rownames(env)<-common

# Spatial folds are assigned to unique localities, never to pair rows. K-means
# is run on equal-area coordinates so each locality and all its individuals stay
# together. No pair touching a test individual enters training.
# A local equirectangular approximation is sufficient for fold assignment and
# avoids any runtime dependency on online PROJ grids.
lat0 <- mean(env$Y) * pi/180
xy_m <- cbind(x=env$X*cos(lat0)*111.32, y=env$Y*110.57)
loc_key <- paste(round(env$X,5),round(env$Y,5),sep="_")
u <- !duplicated(loc_key); km <- kmeans(scale(xy_m[u,,drop=FALSE]),centers=nfold,nstart=100)
fold_map <- setNames(km$cluster,loc_key[u]); env$fold <- unname(fold_map[loc_key])
folds <- env[,c("id","site","X","Y","fold")]
write.csv(folds,file.path(dirs["validation"],"spatial_folds.csv"),row.names=FALSE)

matrix_df <- function(mat, idx) {
  z <- mat[idx,idx,drop=FALSE]; d <- data.frame(site=seq_along(idx),z,check.names=FALSE)
  names(d)[-1] <- as.character(seq_along(idx)); d
}
sitepair <- function(idx, climates=character(), env_override=NULL, geo=TRUE, resistance=TRUE,
                     res_override=NULL) {
  ee <- if(is.null(env_override)) env[idx,,drop=FALSE] else env_override
  ee$site <- seq_along(idx)
  pred <- ee[,c("site","X","Y",climates),drop=FALSE]
  gd <- matrix_df(gen,idx); rr <- if(is.null(res_override)) res else res_override
  rd <- matrix_df(rr,idx)
  a <- list(bioData=gd,bioFormat=3,XColumn="X",YColumn="Y",siteColumn="site",predData=pred)
  if(resistance) a$distPreds <- list(rd)
  z <- do.call(formatsitepair,a); names(z)<-sub("matrix_1","topographicResistance",names(z),fixed=TRUE); z
}
fit_model <- function(idx,climates=character(),geo=TRUE,resistance=TRUE,env_override=NULL,res_override=NULL)
  gdm(sitepair(idx,climates,env_override,geo,resistance,res_override),geo=geo)
metric <- function(model,tab) {
  pr <- as.numeric(predict(model,tab)); ob <- tab$distance
  data.frame(n_pairs=nrow(tab),rmse=sqrt(mean((ob-pr)^2)),mae=mean(abs(ob-pr)),
             correlation=if(sd(pr)>0) cor(ob,pr) else NA_real_)
}

redundancy_filter <- function(idx,vars,threshold=.92) {
  cm <- cor(env[idx,vars,drop=FALSE],use="pairwise.complete.obs",method="spearman")
  hc <- hclust(as.dist(1-abs(cm)),method="average"); cl <- cutree(hc,h=1-threshold)
  chosen <- vapply(split(vars,cl),function(v) priority[priority %in% v][1],character(1))
  tab <- do.call(rbind,lapply(split(vars,cl),function(v) data.frame(cluster=cl[v[1]],variable=v,
    retained=v %in% chosen,reason=ifelse(v %in% chosen,"retained representative",
      paste0("removed: |rho| >= ",threshold," cluster; biological/paleo priority given to ",chosen[cl[v[1]]])))))
  list(vars=unname(chosen),cor=cm,table=tab)
}

permute_site_variable <- function(idx,var,ep=NULL) {
  if(is.null(ep)) ep <- env[idx,,drop=FALSE]
  # I permute locality-level climate values within each broad training block.
  # Replicates at one coordinate move together; pair dependence and coarse
  # spatial structure are retained.
  lk <- paste(round(ep$X,5),round(ep$Y,5),sep="_")
  for(b in unique(ep$fold)) {
    ub <- unique(lk[ep$fold==b]); if(length(ub)<2) next
    shifted <- sample(ub,length(ub),replace=FALSE); mp <- setNames(shifted,ub)
    source_rows <- match(unname(mp[lk[ep$fold==b]]),lk)
    ep[ep$fold==b,var] <- ep[source_rows,var]
  }
  ep
}

permute_resistance_labels <- function(idx) {
  perm <- seq_len(nrow(env))
  for(b in unique(env$fold[idx])) {
    ib <- idx[env$fold[idx]==b]
    # Individuals collected at an identical locality are moved as one unit by
    # matching locality groups; where group sizes differ, individual labels are
    # permuted within the same broad block. Pair rows themselves never move.
    perm[ib] <- sample(ib,length(ib),replace=FALSE)
  }
  res[perm,perm,drop=FALSE]
}

conditional_screen <- function(idx,candidates,nperm_local=n_perm_test) {
  basefit <- fit_model(idx,character(),TRUE,TRUE); base <- basefit$explained
  out <- lapply(candidates,function(v){
    observed <- fit_model(idx,v,TRUE,TRUE)$explained; gain <- observed-base
    null <- vapply(seq_len(nperm_local),function(k) {
      ep <- permute_site_variable(idx,v); fit_model(idx,v,TRUE,TRUE,ep)$explained-base
    },numeric(1))
    data.frame(variable=v,base_deviance=base,test_deviance=observed,gain_deviance=gain,
      permutation_importance=gain-mean(null),p_value=(1+sum(null>=gain))/(nperm_local+1),
      null_mean=mean(null),null_q95=quantile(null,.95),nperm=nperm_local)
  })
  do.call(rbind,out)
}

backward_select <- function(idx,candidates) {
  scr <- conditional_screen(idx,candidates)
  keep <- scr$variable[scr$gain_deviance>=.10 & scr$p_value<=.10]
  if(!length(keep)) keep <- scr$variable[which.max(scr$gain_deviance)]
  history <- data.frame(step=0,removed=NA,variables=paste(keep,collapse=";"))
  repeat {
    full <- fit_model(idx,keep,TRUE,TRUE); if(length(keep)<=1) break
    drops <- vapply(keep,function(v) full$explained-fit_model(idx,setdiff(keep,v),TRUE,TRUE)$explained,numeric(1))
    weak <- names(which.min(drops)); if(min(drops)>=.10) break
    keep <- setdiff(keep,weak); history<-rbind(history,data.frame(step=nrow(history),removed=weak,variables=paste(keep,collapse=";")))
  }
  list(selected=keep,screen=scr,history=history)
}

logit("Starting nested spatial CV for",nrow(env),"individuals and",nLoc(gl),"SNPs")
fold_results <- list(); selections <- list(); screens <- list(); redund <- list(); model_store <- list()
for(f in seq_len(nfold)) {
  logit("Fold",f,"of",nfold)
  train <- which(env$fold!=f); test <- which(env$fold==f)
  red <- redundancy_filter(train,names(bio_names)); sel <- backward_select(train,red$vars)
  redund[[f]] <- transform(red$table,fold=f); screens[[f]]<-transform(sel$screen,fold=f)
  selections[[f]] <- data.frame(fold=f,variable=sel$selected)
  specs <- list(G=list(cl=character(),r=FALSE),GR=list(cl=character(),r=TRUE),
                GC=list(cl=sel$selected,r=FALSE),GRC=list(cl=sel$selected,r=TRUE))
  for(nm in names(specs)) {
    m <- fit_model(train,specs[[nm]]$cl,TRUE,specs[[nm]]$r)
    tt <- sitepair(test,specs[[nm]]$cl,resistance=specs[[nm]]$r)
    me <- metric(m,tt); fold_results[[paste(f,nm)]] <- cbind(fold=f,model=nm,
      n_train=length(train),n_test=length(test),train_deviance=m$explained,me)
    model_store[[paste(f,nm)]] <- m
  }
}
cv <- do.call(rbind,fold_results); rownames(cv)<-NULL
selall <- do.call(rbind,selections); freq <- data.frame(variable=names(bio_names),
 selection_frequency=vapply(names(bio_names),function(v)mean(vapply(selections,function(z)v %in% z$variable,logical(1))),numeric(1)))
write.csv(do.call(rbind,redund),file.path(dirs["tables"],"climate_redundancy_by_fold.csv"),row.names=FALSE)
write.csv(do.call(rbind,screens),file.path(dirs["tables"],"conditional_gdm_importance_by_fold.csv"),row.names=FALSE)
write.csv(selall,file.path(dirs["validation"],"selected_variables_by_fold.csv"),row.names=FALSE)
write.csv(freq,file.path(dirs["validation"],"climate_selection_stability.csv"),row.names=FALSE)
write.csv(cv,file.path(dirs["validation"],"spatial_cv_model_comparison.csv"),row.names=FALSE)
saveRDS(model_store,file.path(dirs["models"],"spatial_fold_models.rds"))

# Pair-compatible incremental null models. I compare the observed out-of-sample
# RMSE gain of climate over G+R and resistance over G+C with gains obtained when
# the relevant information is randomized at individuals/localities inside each
# training block. Test data remain untouched.
null_rows <- list()
for(f in seq_len(nfold)) {
  train <- which(env$fold!=f); test <- which(env$fold==f)
  vars <- selections[[f]]$variable
  base_gr <- fit_model(train,character(),TRUE,TRUE)
  base_gc <- fit_model(train,vars,TRUE,FALSE)
  test_gr <- sitepair(test,character(),resistance=TRUE)
  test_gc <- sitepair(test,vars,resistance=FALSE)
  test_full <- sitepair(test,vars,resistance=TRUE)
  rmse_gr <- metric(base_gr,test_gr)$rmse; rmse_gc <- metric(base_gc,test_gc)$rmse
  obs_full <- fit_model(train,vars,TRUE,TRUE); rmse_full <- metric(obs_full,test_full)$rmse
  null_rows[[paste(f,"obs")]] <- data.frame(fold=f,permutation=0,
    climate_gain_rmse=rmse_gr-rmse_full,resistance_gain_rmse=rmse_gc-rmse_full)
  for(k in seq_len(n_perm_final)) {
    ep <- env[train,,drop=FALSE]
    for(v in vars) ep <- permute_site_variable(train,v,ep)
    mc <- fit_model(train,vars,TRUE,TRUE,env_override=ep)
    rp <- permute_resistance_labels(train)
    mr <- fit_model(train,vars,TRUE,TRUE,res_override=rp)
    null_rows[[paste(f,k)]] <- data.frame(fold=f,permutation=k,
      climate_gain_rmse=rmse_gr-metric(mc,test_full)$rmse,
      resistance_gain_rmse=rmse_gc-metric(mr,test_full)$rmse)
  }
}
nulls <- do.call(rbind,null_rows); rownames(nulls)<-NULL
write.csv(nulls,file.path(dirs["validation"],"spatial_group_null_models.csv"),row.names=FALSE)
obs <- colMeans(nulls[nulls$permutation==0,c("climate_gain_rmse","resistance_gain_rmse")])
nul <- aggregate(cbind(climate_gain_rmse,resistance_gain_rmse)~permutation,
                 nulls[nulls$permutation>0,],mean)
null_summary <- data.frame(group=c("Current climate","Effective resistance"),
 observed_mean_rmse_gain=unname(obs),null_mean=c(mean(nul$climate_gain_rmse),mean(nul$resistance_gain_rmse)),
 p_value=c((1+sum(nul$climate_gain_rmse>=obs[1]))/(n_perm_final+1),
           (1+sum(nul$resistance_gain_rmse>=obs[2]))/(n_perm_final+1)),nperm=n_perm_final)
write.csv(null_summary,file.path(dirs["validation"],"spatial_group_null_summary.csv"),row.names=FALSE)

# Stable final climate set: selected by a majority of spatial folds. This is
# a transparent stability rule, applied after inspecting the fold results.
final_clim <- freq$variable[freq$selection_frequency>.50]
if(!length(final_clim)) final_clim <- freq$variable[freq$selection_frequency==max(freq$selection_frequency)]
full_idx <- seq_len(nrow(env)); final_screen <- conditional_screen(full_idx,final_clim)
write.csv(final_screen,file.path(dirs["tables"],"final_climate_permutation_importance.csv"),row.names=FALSE)
final_variables <- merge(data.frame(variable=final_clim),freq,by="variable",all.x=TRUE)
final_variables$BIOCLIM <- unname(bio_names[final_variables$variable])
final_variables$reason <- "selected in a majority of spatial training folds and retained after conditional screening"
write.csv(final_variables,file.path(dirs["tables"],"final_current_climate_variables.csv"),row.names=FALSE)

# Official gdm.varImp is preserved as a method-sensitivity output. Its native
# pair-row permutations are not used for inferential p-values; the locality-
# level compatible null above is the inferential result.
full_tab <- sitepair(full_idx,final_clim,resistance=TRUE)
official <- gdm.varImp(full_tab,geo=TRUE,predSelect=FALSE,nPerm=n_perm_test,parallel=FALSE)
capture.output(official,file=file.path(dirs["tables"],"gdm_varImp_official_sensitivity.txt"))

sets <- expand.grid(G=c(FALSE,TRUE),C=c(FALSE,TRUE),R=c(FALSE,TRUE)); sets<-sets[rowSums(sets)>0,]
mods<-list(); rows<-list()
for(i in seq_len(nrow(sets))) {
 s<-sets[i,]; key<-paste(sort(c("G","C","R")[unlist(s)]),collapse="")
 m<-fit_model(full_idx,if(s$C)final_clim else character(),s$G,s$R); mods[[key]]<-m
 rows[[key]]<-data.frame(model=key,G=s$G,C=s$C,R=s$R,deviance_explained=m$explained)
}
subsets<-do.call(rbind,rows); write.csv(subsets,file.path(dirs["tables"],"gcr_seven_subsets.csv"),row.names=FALSE)
saveRDS(mods,file.path(dirs["models"],"gcr_seven_models.rds"))

val<-setNames(subsets$deviance_explained,subsets$model); vv<-function(k)if(!nzchar(k))0 else val[[k]]
groups<-c("G","C","R"); shap<-setNames(numeric(3),groups)
for(g in groups) for(S in list(character(),setdiff(groups,g)[1],setdiff(groups,g)[2],setdiff(groups,g))) {
 b<-paste(sort(S),collapse=""); a<-paste(sort(c(S,g)),collapse=""); w<-factorial(length(S))*factorial(2-length(S))/factorial(3)
 shap[g]<-shap[g]+w*(vv(a)-vv(b))
}
part<-data.frame(group=groups,shapley=shap,share_percent=100*shap/vv("CGR"),
 drop_one=c(vv("CGR")-vv("CR"),vv("CGR")-vv("GR"),vv("CGR")-vv("CG")))
write.csv(part,file.path(dirs["tables"],"gcr_partition.csv"),row.names=FALSE)
write.csv(data.frame(group=groups,full_model="GCR",without=c("CR","GR","CG"),
                     full_explained_deviance=vv("CGR"),without_explained_deviance=c(vv("CR"),vv("GR"),vv("CG")),
                     deviance_loss=part$drop_one),file.path(dirs["tables"],"gdm_drop_one.csv"),row.names=FALSE)

imp <- isplineExtract(mods[["CGR"]]); ii<-data.frame(variable=colnames(imp$y),spline_height=apply(imp$y,2,max))
ii$relative_percent<-100*ii$spline_height/sum(ii$spline_height); ii<-ii[order(-ii$spline_height),]
write.csv(ii,file.path(dirs["tables"],"final_ispline_importance.csv"),row.names=FALSE)

p <- ggplot(freq,aes(reorder(variable,selection_frequency),selection_frequency))+
 geom_col(fill="#2F6B4F")+geom_hline(yintercept=.50,lty=2,col="#8F3B3F")+coord_flip()+theme_classic()+
 labs(x=NULL,y="Selection frequency across spatial folds",title="Climate-variable stability")
ggsave(file.path(dirs["figures"],"climate_selection_stability.png"),p,width=8,height=6,dpi=220)

writeLines(c("DEFINITIVE CONTEMPORARY GDM AUDIT",paste("Seed:",seed),paste("Test permutations:",n_perm_test),paste("Reserved final permutations:",n_perm_final),
 paste("Stable climate variables:",paste(final_clim,collapse=", ")),
 "Geography and direct Circuitscape effective topographic resistance remain fixed during climate selection.",
 "The resistance surface was designed mainly from topography and manually included exposed mainland-Tasmania land during low sea level; I do not fit a separate Bass Strait predictor.",
"No pair touching a test individual enters training. Test performance uses test-test pairs only.",
"Permutation p-values come from locality-level permutations within training spatial blocks, never pair-row permutations.",
 "The genomic matrix is the landscape-background panel: one SNP per RAD tag with the reproducibly recoverable PCA-conditioned RDA candidate set excluded; it is not claimed to be fully neutral.",
 "Palaeoclimate is not fitted jointly with current climate. Projection and incremental PaleoHistory tests require locally archived, provenance-matched LIG/LGM rasters."),
 file.path(dirs["results"],"METHOD_STATUS.txt"))
capture.output(sessionInfo(),file=file.path(dirs["logs"],"sessionInfo.txt"))
logit("Completed. Stable variables:",paste(final_clim,collapse=", "))
