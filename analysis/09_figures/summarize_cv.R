# 1. Read predictions from random and spatial folds ----
args <- commandArgs(TRUE); out <- args[1]
dir <- file.path(out,"tables")
r <- do.call(rbind,lapply(sort(list.files(file.path(out,"checkpoints"),"^prediction_..\\.csv$",full.names=TRUE)),
                         read.csv,colClasses=c(id1="character",id2="character")))
s <- read.csv(file.path(dir,"spatial_cv_test_pair_predictions.csv"),colClasses=c(id1="character",id2="character"))
a <- read.csv(file.path(dir,"random_fold_assignments.csv"))
done <- sort(unique(r$repetition))
stopifnot(length(done)>=5,!anyDuplicated(r[,c("repetition","pair_id")]),!anyDuplicated(s$pair_id))
# 2. Match observations by pair and sampling region ----
stable <- c("id1","id2","locality1","locality2","region1","region2")
mid <- match(r$pair_id,s$pair_id); within <- r$region1==r$region2
stopifnot(all(!is.na(mid[within])),all(is.na(mid[!within])),
          all(r[within,stable]==s[mid[within],stable]),
          max(abs(r$observed[within]-s$observed[mid[within]]))<1e-12)
for(n in c("m2","m3")) r[[paste0(n,"_spatial")]] <- s[[n]][mid]
# 3. Calculate prediction error and geographic distances ----
rmse <- function(p,y) sqrt(mean((p-y)^2))
calc <- function(x,sp=FALSE) {
  n <- if(sp) c("m2_spatial","m3_spatial") else c("m2","m3")
  c(m2=rmse(x[[n[1]]],x$observed),m3=rmse(x[[n[2]]],x$observed))
}
region <- list(); nation <- list(); diagnostics <- list(); distances <- list()
add <- function(r,nm,aggregation,design,val,n_pairs,n_units) {
  data.frame(repetition=r,target=nm,aggregation=aggregation,design=design,
             rmse_m2=val[["m2"]],rmse_m3=val[["m3"]],gain=val[["m2"]]-val[["m3"]],n_pairs=n_pairs,n_units=n_units)
}
d_km <- function(x,y,xs,ys) {
  rad <- pi/180
  z <- sin((ys-y)*rad/2)^2+cos(y*rad)*cos(ys*rad)*sin((xs-x)*rad/2)^2
  6371.0088*2*asin(sqrt(pmin(1,z)))
}
write_csv <- function(x,name) write.csv(x,file.path(dir,name),row.names=FALSE)
num_ind <- table(c(setNames(s$locality1,s$id1),setNames(s$locality2,s$id2))[
  !duplicated(names(c(setNames(s$locality1,s$id1),setNames(s$locality2,s$id2))))])
local_dist <- function(loc,folds,label,b) {
  data.frame(repetition=b,design=label,locality=loc$locality,fold=folds,region=loc$fold,
    nearest_km=vapply(seq_len(nrow(loc)),function(i) {
      t <- folds!=folds[i]; min(d_km(loc$X[i],loc$Y[i],loc$X[t],loc$Y[t]))
    },numeric(1)))
}
u0 <- subset(a,repetition==1)
distances[[1]] <- local_dist(u0,u0$fold,"spatial",0L)
for(b in done) {
  y <- subset(r,repetition==b)
  u <- subset(a,repetition==b)
  distances[[length(distances)+1L]] <- local_dist(u,u$random_fold,"random",b)
  fold_values <- do.call(rbind,lapply(split(y,y$fold),calc))
  nation[[length(nation)+1L]] <- add(b,"all_test_pairs","equal_random_fold","random",colMeans(fold_values),nrow(y),5)
  nation[[length(nation)+1L]] <- add(b,"all_test_pairs","pair_weighted","random",calc(y),nrow(y),5)
  for(f in 1:5) {
    ff <- subset(u,random_fold==f); pred <- subset(y,fold==f)
    diagnostics[[length(diagnostics)+1L]] <- data.frame(repetition=b,fold=f,localities=nrow(ff),
      individuals=sum(num_ind[ff$locality]),regions_represented=length(unique(ff$fold)),
      lon_min=min(ff$X),lon_max=max(ff$X),lat_min=min(ff$Y),lat_max=max(ff$Y),n_test_pairs=nrow(pred),
      n_within_region=sum(pred$region1==pred$region2))
  }
  keep <- list()
  for(f in 1:5) {
    z <- subset(y,region1==f & region2==f)
    loci <- unique(c(z$locality1,z$locality2))
    ll <- with(z,paste(pmin(locality1,locality2),pmax(locality1,locality2)))
    cross <- length(unique(ll[z$locality1!=z$locality2]))
    good <- length(loci)>=3 && cross>=2
    vals <- if(good) c(calc(z),calc(z,TRUE)) else rep(NA_real_,4)
    if(good) keep[[as.character(f)]] <- z
    region[[length(region)+1L]] <- data.frame(repetition=b,region=f,usable=good,pairs=nrow(z),
      ids=length(unique(c(z$id1,z$id2))),localities=length(loci),distinct_interlocality_pairs=cross,
      matched_coverage=nrow(z)/sum(s$fold==f),random_m2=vals[1],random_m3=vals[2],spatial_m2=vals[3],spatial_m3=vals[4],
      gain_random=vals[1]-vals[2],gain_spatial=vals[3]-vals[4],delta_gain=vals[1]-vals[2]-vals[3]+vals[4],row.names=NULL)
  }
  if(length(keep)) for(sp in c(FALSE,TRUE)) {
    full <- do.call(rbind,keep)
    for(aggregation in c("equal_original_region","pair_weighted")) {
      x <- if(aggregation=="equal_original_region") colMeans(do.call(rbind,lapply(keep,calc,sp=sp))) else calc(full,sp)
      nation[[length(nation)+1L]] <- add(b,"matched_test_pairs",aggregation,if(sp) "spatial" else "random",x,nrow(full),length(keep))
    }
  }
}
nation[[length(nation)+1L]] <- add(0,"all_original_spatial_pairs","equal_original_region","spatial",
                                 colMeans(do.call(rbind,lapply(split(s,s$fold),calc))),nrow(s),5)
national <- do.call(rbind,nation)
regional <- do.call(rbind,region)
write_csv(r,"random_cv_test_pair_predictions.csv")
write_csv(national,"national_cv_metrics.csv")
write_csv(regional,"regional_matched_cv_metrics.csv")
write_csv(do.call(rbind,diagnostics),"random_cv_fold_diagnostics.csv")
write_csv(do.call(rbind,distances),"test_train_geographic_distances.csv")
# 4. Summarise errors and model gains across replicates ----
summary <- function(x) {
  a <- x[is.finite(x)]
  if(!length(a)) return(rep(NA_real_,5))
  unname(c(min(a),quantile(a,c(.025,.50,.975)),max(a)))
}
summarize <- function(x,values,id) {
  ans <- do.call(rbind,lapply(split(x,interaction(x[id],drop=TRUE)),function(d) {
    do.call(rbind,lapply(values,function(v) setNames(
      cbind(d[1,id,drop=FALSE],metric=v,n=sum(is.finite(d[[v]])),t(summary(d[[v]]))),
      c(id,"metric","n","minimum","q025","median","q975","maximum"))))
  }))
  rownames(ans)<-NULL; ans
}
write_csv(summarize(national,c("rmse_m2","rmse_m3","gain"),c("target","aggregation","design")),"random_cv_summary.csv")
write_csv(summarize(regional,c("delta_gain","gain_spatial","gain_random","matched_coverage"),"region"),"regional_cv_summary.csv")
