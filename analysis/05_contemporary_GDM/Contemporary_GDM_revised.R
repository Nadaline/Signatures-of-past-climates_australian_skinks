#### 1. Frozen contemporary fits, folds and conditional nulls ####
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L || (length(args) == 3L && args[3] == "--figure-only"))
project <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = TRUE)
d <- file.path(project,"results/landscape/gdm_definitive")
subsets <- read.csv(file.path(d,"tables/gcr_seven_subsets.csv"))
splines <- read.csv(file.path(d,"tables/final_isplines.csv"))
clim <- read.csv(file.path(d,"tables/final_current_climate_variables.csv"))
cv <- read.csv(file.path(d,"validation/spatial_cv_model_comparison.csv"))
folds <- read.csv(file.path(d,"validation/spatial_folds.csv"))
null <- read.csv(file.path(d,"validation/spatial_group_null_models.csv"))
tests <- read.csv(file.path(d,"validation/spatial_group_null_summary.csv"))
audit <- read.csv(file.path(d,"tables/genomic_dataset_audit.csv"))
stopifnot(audit$n_ind == 170, audit$n_snp == 19131, choose(audit$n_ind,2) == 14365,
  nrow(folds) == 170, setequal(folds$fold,1:5),
  all(table(cv$fold,cv$model) == 1L), setequal(cv$fold,1:5),
  setequal(clim$BIOCLIM,c("bio2","bio4","bio12")), setequal(clim$variable,c("tempRange","tempSeason","precAnn")),
  abs(max(dist(folds[,c("X","Y")]))-max(splines$x_Geographic)) < 1e-9,
  setequal(null$permutation, 0:999), nrow(null)==5000L, all(table(null$fold,null$permutation)==1L))
null_mean <- aggregate(cbind(climate_gain_rmse,resistance_gain_rmse)~permutation,null,mean)
gains <- as.matrix(null_mean[null_mean$permutation > 0,-1])
observed <- as.numeric(null_mean[null_mean$permutation == 0,-1])
tests <- tests[match(c("Current climate","Effective resistance"),tests$group),]
p <- (1+colSums(sweep(gains,2,observed,">=")))/(1+nrow(gains))
cv_mean <- tapply(cv$rmse,cv$model,mean)
stopifnot(max(abs(observed-tests$observed_mean_rmse_gain)) < 1e-13,
  max(abs(observed-c(cv_mean[["GR"]],cv_mean[["GC"]])+cv_mean[["GRC"]])) < 1e-13,
  identical(unname(p),tests$p_value), identical(unname(p),c(.141,.019)),
  all(tests$nperm==999))

#### 2. Comparison table; retain missing validation metrics ####
idx <- match(c("G","C","R","CG","GR","CR","CGR"),subsets$model)
stopifnot(!anyNA(idx), !anyDuplicated(idx))
tab <- data.frame(Model = c("G","C","R","G + C","G + R","C + R","G + C + R"),
  Components = c("Geography","Climate","Resistance","Geography, climate",
    "Geography, resistance","Climate, resistance","Geography, climate, resistance"),
  Training_deviance_explained_pct = subsets$deviance_explained[idx],
  Mean_spatial_holdout_RMSE = unname(cv_mean[c("G","C","R","GC","GR","CR","GRC")]))
if (length(args) == 2L) {
write.csv(tab,file.path(out,"Supplementary_GDM_model_comparison.csv"),row.names=FALSE,na="")
word_src <- tempfile(fileext=".md")
word <- file.path(out,"Supplementary_GDM_model_comparison.docx")
md_rows <- sprintf("| %s | %s | %.2f | %s |",tab$Model,tab$Components,
  tab$Training_deviance_explained_pct,ifelse(is.na(tab$Mean_spatial_holdout_RMSE),"—",
  sprintf("%.4f",tab$Mean_spatial_holdout_RMSE)))
writeLines(c("**Supplementary table. Contemporary GDM model comparison**","",
  "| Model | Components included | Training deviance explained (%) | Mean spatial hold-out RMSE |",
  "| :--------------- | :-------------------------------- | ----------------------: | ----------------------: |",
  md_rows,"",
  "Note: Deviance explained refers to the models fitted to the full contemporary dataset. RMSE is the equally weighted mean across five spatial hold-out folds defined by sampling locality. A dash indicates that no spatial hold-out RMSE was archived for that model. G, geographic separation; C, present-day climate; R, effective topographic resistance."),word_src)
pandoc <- Sys.which("pandoc")
stopifnot(nzchar(pandoc), system2(pandoc,c("-f","markdown","-t","docx",shQuote(word_src),"-o",shQuote(word))) == 0L)
}

#### 3. Five saved I-splines and two conditional permutation panels ####
ink <- "#251462"
electric_violet <- "#6811DB"
amber <- "#EBB201"
names_x <- c("Geographic distance","Topographic resistance","Annual precipitation",
  "Diurnal temperature range","Temperature seasonality")
units_x <- c("Decimal degrees","Circuitscape resistance","BIO12 (mm yr⁻¹)","BIO2 (°C)","BIO4 (°C)")
vars <- c("Geographic","topographicResistance","precAnn","tempRange","tempSeason")
stopifnot(all(paste0("x_",vars) %in% names(splines)),all(paste0("y_",vars) %in% names(splines)))
yt <- c(0,.1,.2,.3,.4)
draw <- function() {
  par(mar=rep(0,4),family="sans",xaxs="i",yaxs="i")
  plot.new(); plot.window(c(0,450),c(0,450))
  text(12,433,"A  Fitted GDM I-splines",adj=0,cex=1.16,font=2,col=ink)
  text(15,283,"Cumulative fitted turnover",srt=90,cex=.96,col=ink)
  text(12,149,"B  Conditional permutation tests",adj=0,cex=1.16,font=2,col=ink)
  text(15,77,"Permutations",srt=90,cex=.96,col=ink)
  text(238,12,"Held-out RMSE improvement",cex=.96,col=ink)
  pos <- rbind(c(27,163,287,417),c(169,305,287,417),c(311,447,287,417),
    c(98,234,160,290),c(240,376,160,290))
  par(fg=ink,col.axis=ink,col.lab=ink,col.main=ink,las=1,mgp=c(1.55,.43,0),tcl=-.23)
  for (j in seq_along(vars)) {
    par(new=TRUE,fig=pos[j,]/450,mar=c(3.0,3.03,2.6,.61),xaxs="r",yaxs="i")
    xx <- splines[[paste0("x_",vars[j])]]
    yy <- splines[[paste0("y_",vars[j])]]
    stopifnot(all(diff(xx)>=0), all(diff(yy)>=-1e-12))
    plot(xx,yy,type="l",ylim=c(0,.44),axes=FALSE,xlab=units_x[j],ylab="",
      col=amber,lwd=1.48,cex.lab=.87)
    tick <- pretty(range(xx),3)
    tick <- tick[tick>=min(xx) & tick<=max(xx)]
    # The local BIO4 input is in hundredths of °C; leave the spline and tick positions unchanged.
    tick_label <- if (vars[j] == "tempSeason") tick/100 else tick
    axis(1,at=tick,labels=format(tick_label,scientific=FALSE,trim=TRUE),cex.axis=.84,lwd=.59)
    axis(2,at=yt,labels=formatC(yt,digits=1,format="f"),cex.axis=.83,lwd=.59)
    mtext(names_x[j],side=3,adj=0,line=.92,cex=.93,font=2)
  }
  hp <- rbind(c(27,229,21,140),c(239,444,21,140))
  htitle <- c("Climate: G + R  →  G + C + R","Resistance: G + C  →  G + C + R")
  for (j in seq_len(2L)) {
    par(new=TRUE,fig=hp[j,]/450,mar=c(2.45,3.05,3.0,.78),xaxs="r",yaxs="i")
    hist_j <- hist(gains[,j],breaks=32,plot=FALSE)
    hmax <- max(hist_j$counts)
    xr <- range(hist_j$breaks,observed[j])
    plot(hist_j,freq=TRUE,axes=FALSE,ann=FALSE,ylim=c(0,hmax*1.14),xlim=xr,
      col=electric_violet,border="white",lwd=.43)
    abline(v=tests$observed_mean_rmse_gain[j],col=amber,lwd=1.45)
    uu <- par("usr")
    text(uu[2]-.028*diff(uu[1:2]),hmax*1.08,sprintf("P = %.3f",tests$p_value[j]),
      adj=c(1,.5),col=ink,cex=.90)
    at <- pretty(xr,4)
    axis(1,at=at[at>=xr[1] & at<=xr[2]],cex.axis=.88,lwd=.59)
    axis(2,cex.axis=.84,lwd=.59)
    mtext(htitle[j],side=3,adj=0,line=.85,cex=.91,font=2)
  }
}
cairo_pdf(file.path(out,"Contemporary_GDM_revised.pdf"),width=450/72,height=450/72,
  pointsize=9,bg="white")
draw(); invisible(dev.off())
png(file.path(out,"Contemporary_GDM_revised.png"),width=450/72,height=450/72,
  units="in",type="cairo",res=600,pointsize=9,bg="white")
draw(); invisible(dev.off())
print(tests,digits=16,row.names=FALSE)
