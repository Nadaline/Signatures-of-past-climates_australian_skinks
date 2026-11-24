# Figure 5: additional predictive contribution of historical climates
# Display the contribution of Last Interglacial and Last Glacial Maximum environments to genetic dissimilarity after accounting for contemporary landscape predictors. This script visualizes model comparisons using saved predictions.
# Compare the same held-out localities across models; palaeoclimate predictors describe environmental contrasts and do not date population splits.
#
#### 1. Historical-climate model results ####
library(sf)
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
project <- normalizePath(args[1], mustWork = TRUE)
basemap <- normalizePath(args[2], mustWork = TRUE)
out <- normalizePath(args[3], mustWork = TRUE)
d <- file.path(project, "results/landscape/paleoclimate")
cv <- read.csv(file.path(d, "validation/historical_model_comparison.csv"))
fold <- read.csv(file.path(d, "null_corrected/corrected_observed_fold_effects.csv"))
null <- read.csv(file.path(d, "null_corrected/corrected_historical_null_statistics.csv"))
obs <- read.csv(file.path(d, "null_corrected/corrected_historical_null_summary.csv"))
sites <- read.csv(file.path(project, "results/landscape/gdm_definitive/validation/spatial_folds.csv"))
gain <- read.csv(file.path(d, "spatial_effect/tables/gdm_historical_gain_by_locality.csv"))
tas <- read.csv(file.path(d, "spatial_effect/tables/gdm_tasmania_rmse_sensitivity.csv"))
pairs <- read.csv(file.path(d, "spatial_effect/tables/gdm_pairwise_historical_gain.csv"))
fold <- fold[order(fold$fold), ]
ma <- cv[cv$model == "M2", ]; ma <- ma[order(ma$fold), ]
mb <- cv[cv$model == "M3_PaleoMean", ]; mb <- mb[order(mb$fold), ]
rmse <- c(mean(ma$rmse), mean(mb$rmse))
n_exceed <- sum(null$Delta_RMSE_null >= obs$observed_Delta_RMSE)
stopifnot(identical(fold$fold, 1:5), identical(ma$fold, fold$fold), identical(mb$fold, fold$fold),
  max(abs(ma$rmse - mb$rmse - fold$observed_Delta_RMSE)) < 1e-12,
  abs(rmse[1] - rmse[2] - obs$observed_Delta_RMSE) < 1e-12,
  obs$n_perm == 999, setequal(null$permutation_id, seq_len(obs$n_perm)),
  nrow(null) == obs$n_perm, (1+n_exceed)/(1+nrow(null)) == obs$p_value, obs$p_value == .035,
  nrow(sites) == 170L, !anyDuplicated(sites$id))
loc <- unique(sites[, c("X", "Y", "fold")])
k <- function(x) paste(x$X, x$Y, sep = "/")
stopifnot(nrow(loc) == 101L, nrow(gain) == nrow(loc), !anyDuplicated(k(loc)), !anyDuplicated(k(gain)),
  setequal(k(loc), k(gain)), all(loc$fold[match(k(gain), k(loc))] == gain$fold))
tas <- tas[match(c("all_heldout_pairs", "without_Tasmania_pairs", "without_Tasmania_containing_fold"),
  tas$diagnostic), ]
stopifnot(!anyNA(tas$Delta_RMSE), all(tas$Delta_RMSE > 0), all(is.finite(gain$mean_pairwise_gain_squared_error)))
# Verify the plotted summaries from saved prediction errors, without fitting a model.
rmse_pair <- function(p) sqrt(mean(p$squared_error_M2)) - sqrt(mean(p$squared_error_M3))
tas_fold <- unique(pairs$fold[pairs$lat1 < -39 | pairs$lat2 < -39])
stopifnot(identical(tas_fold, 1L))
subsets <- list(pairs, pairs[pairs$lat1 >= -39 & pairs$lat2 >= -39, ], pairs[pairs$fold != tas_fold, ])
stopifnot(all(vapply(subsets, nrow, integer(1)) == tas$n_pairs),
  max(abs(vapply(subsets, rmse_pair, numeric(1)) - tas$Delta_RMSE)) < 1e-12,
  max(abs(vapply(split(pairs, pairs$fold), rmse_pair, numeric(1)) - fold$observed_Delta_RMSE)) < 1e-12)
se_diff <- pairs$squared_error_M2 - pairs$squared_error_M3
endpoint <- rbind(data.frame(X=pairs$lon1, Y=pairs$lat1, value=se_diff),
  data.frame(X=pairs$lon2, Y=pairs$lat2, value=se_diff))
local_means <- tapply(endpoint$value, k(endpoint), mean)
stopifnot(setequal(names(local_means), k(gain)),
  max(abs(unname(local_means[k(gain)]) - gain$mean_pairwise_gain_squared_error)) < 1e-12)

#### 2. Shared projected coast and display scales ####
crs_map <- 9473
states <- st_read(basemap, quiet = TRUE)
states <- states[states$state_code_2016 %in% as.character(1:8), ]
states <- st_cast(suppressWarnings(st_collection_extract(
  st_make_valid(st_transform(states, crs_map)), "POLYGON")), "MULTIPOLYGON")
xy <- function(x) unname(st_coordinates(st_transform(st_as_sf(x, coords = c("X","Y"),
  crs = 4326), crs_map)))
loc_xy <- xy(loc)
gain_xy <- loc_xy[match(k(gain), k(loc)), , drop = FALSE]
bb <- c(xmin = min(loc_xy[,1]) - 125000, ymin = min(loc_xy[,2]) - 115000,
  xmax = max(loc_xy[,1]) + 100000, ymax = max(loc_xy[,2]) + 115000)
clip <- st_as_sfc(st_bbox(bb, crs = st_crs(crs_map)))
land <- suppressWarnings(st_intersection(st_geometry(states), clip))
coast <- suppressWarnings(st_intersection(st_boundary(st_geometry(states)), clip))
st_lab <- xy(data.frame(X=c(146.5, 145.8, 143.4, 137.3, 149.35),
  Y=c(-22.1, -32.7, -36.75, -32.05, -42.95)))

ink <- "#251462"
violet <- "#6811DB"
amber <- "#EBB201"
fold_col <- c(violet, "#277B9E", "#C78724", "#A4526C", "#B8D6E2")
fold_pch <- c(21, 22, 24, 23, 25)
z <- gain$mean_pairwise_gain_squared_error
lmt <- ceiling(max(abs(z))*1000)/1000
pal <- colorRampPalette(c("#A45244","#F7F4F2","#6811DB"))(257)
gain_col <- function(x) pal[1L + round(256*(x/lmt+1)/2)]
gain_cex <- function(x) sqrt(.52^2 + (1.45^2 - .52^2)*abs(x)/lmt)
page <- c(450, 567)
fig <- function(x1,x2,y1,y2) c(x1,x2,y1,y2)/rep(page, each = 2)
axis_col <- "#716687"

#### 3. Single six-panel plate ####
draw <- function() {
  par(mar=rep(0,4),family="sans",xaxs="i",yaxs="i",fg=ink,col.axis=ink,col.lab=ink)
  plot.new(); plot.window(c(0,page[1]),c(0,page[2]))
  text(c(10,155,304), 550, c("A  Predictive performance","B  Prediction gain across folds","C  Spatial permutation test"),
    adj=0,font=2,cex=1.09,col=ink)
  text(c(10,230), 407, c("D  Spatial validation folds","E  Palaeoclimatic gain by locality"),
    adj=0,font=2,cex=1.09,col=ink)
  text(10, 86, "F  Effect of excluding Tasmania", adj=0,font=2,cex=1.09,col=ink)
  par(las=1,tcl=-.24,mgp=c(1.54,.45,0))

  par(fig=fig(7,148,420,539),new=TRUE,mar=c(2.9,5.4,.7,.65))
  plot.new(); plot.window(c(.1207,.1304),c(.47,2.86))
  segments(.121,c(1,2),.130,c(1,2),col="#E1DBE9",lwd=.65)
  segments(rmse[1],2,rmse[2],1,lwd=1.2,col="#9F8EB9")
  points(rmse,2:1,pch=21,bg=c(ink,violet),col="white",cex=1.08,lwd=.55)
  text(rmse,2:1,sprintf("%.5f",rmse),pos=3,offset=.47,cex=.89,col=ink)
  text(.12555,2.72,sprintf("−%.1f%% RMSE",100*obs$observed_Delta_RMSE/rmse[1]),cex=.84,col=violet)
  axis(1,at=c(.122,.126,.130),cex.axis=.85,col=axis_col,lwd=.60)
  axis(2,at=2:1,labels=c("Present","+ Palaeoclimate"),tick=FALSE,cex.axis=.78,line=-.25)
  mtext("Mean held-out RMSE", side=1, line=1.48, cex=.90)

  par(fig=fig(151,298,420,539),new=TRUE,mar=c(2.9,3.65,.7,.9))
  v <- 1000*fold$observed_Delta_RMSE
  b <- barplot(v,col=fold_col,border=ink,ylim=c(-5.0,9.9),axes=FALSE,
    xlab="",ylab="",space=.44,lwd=.45)
  abline(h=0,col=ink,lwd=.82)
  axis(1,at=b,labels=fold$fold,cex.axis=.87,lwd=.60,col=axis_col)
  axis(2,at=c(-4,0,4,8),labels=c("−4","0","4","8"),cex.axis=.85,lwd=.60,col=axis_col)
  mtext("Spatial fold",1,line=1.48,cex=.90)
  mtext("ΔRMSE (×10⁻³)",2,line=2.05,cex=.87,las=0)

  par(fig=fig(301,446,420,539),new=TRUE,mar=c(2.9,3.4,.7,1.0))
  h <- hist(null$Delta_RMSE_null*1000,breaks=34,plot=FALSE)
  plot(h,axes=FALSE,ann=FALSE,ylim=c(0,max(h$counts)*1.22),col="#CBD5E1",border="white",lwd=.34)
  segments(obs$observed_Delta_RMSE*1000,0,obs$observed_Delta_RMSE*1000,max(h$counts)*1.07,col=amber,lwd=1.8)
  u <- par("usr")
  text(u[2]-.13,max(h$counts)*1.16,sprintf("P = %.3f", obs$p_value),adj=1,cex=.91,col=ink)
  axis(1,at=c(-3,0,3,6),labels=c("−3","0","3","6"),cex.axis=.85,lwd=.60,col=axis_col)
  axis(2,cex.axis=.85,lwd=.60,col=axis_col)
  mtext("Mean ΔRMSE (×10⁻³)",1,line=1.48,cex=.90)
  mtext("Permutations",2,line=2.0,cex=.87,las=0)

  for (panel in 1:2) {
    x0 <- if (panel==1) 9 else 229
    par(fig=fig(x0,x0+211,143,397),new=TRUE,mar=rep(0,4),xaxs="i",yaxs="i")
    plot.new()
    plot.window(bb[c("xmin","xmax")],bb[c("ymin","ymax")],asp=1)
    plot(land,add=TRUE,col="#F1EFEB",border=NA)
    plot(coast,add=TRUE,col="#9A94A1",lwd=.54)
    text(st_lab[,1],st_lab[,2],c("QLD","NSW","VIC","SA","TAS"),col="#615773",cex=.78)
    if (panel == 1) {
      points(loc_xy,pch=fold_pch[loc$fold],bg=fold_col[loc$fold],col=ink,cex=.80,lwd=.57)
      x <- par("usr")[1] + 140000; y <- bb["ymin"] + 170000
      segments(x,y,x+300000,y,col=ink,lwd=.85)
      segments(x+c(0,300000),y-24000,x+c(0,300000),y+24000,col=ink,lwd=.85)
      text(x+150000,y-69000,"300 km",cex=.69,col=ink)
    } else {
      j <- order(abs(z),decreasing=TRUE)
      points(gain_xy[j,,drop=FALSE],pch=21,bg=gain_col(z[j]),col=ink,cex=gain_cex(z[j]),lwd=.52)
    }
  }

  par(fig=fig(0,450,96,143),new=TRUE,mar=rep(0,4))
  plot.new(); plot.window(c(0,450),c(0,47))
  text(19,39,"Spatial fold",adj=0,font=2,cex=.91)
  xx <- seq(23,183,by=40)
  points(xx,rep(24,length(xx)),pch=fold_pch,bg=fold_col,col=ink,cex=1.0,lwd=.56)
  text(xx+8,24,1:5,adj=0,cex=.90)
  text(234,43,"Mean squared-error gain (×10⁻³)",adj=0,font=2,cex=.86)
  x <- seq(237,341,length.out=258)
  rect(x[-258],30,x[-1],36,col=pal,border=NA)
  rect(237,30,341,36,border="#9A94A1",lwd=.47)
  text(c(237,289,341),24,c(paste0("−",lmt*1000),"0",paste0("+",lmt*1000)),cex=.77)
  text(234,9,"Absolute gain (×10⁻³)",adj=0,font=2,cex=.83)
  mag <- c(.001,.003,.006)
  mx <- c(326,364,402)
  points(mx,rep(9,length(mx)),pch=21,bg="#F7F4F2",col=ink,cex=gain_cex(mag),lwd=.56)
  text(mx+8,9,mag*1000,adj=0,cex=.78)

  par(fig=fig(10,442,1,78),new=TRUE,mar=c(2.83,13.7,.40,.65))
  tv <- rev(tas$Delta_RMSE)*1000
  ys <- barplot(tv,horiz=TRUE,col=c("#9B66D8",violet,ink),border=NA,
    xlim=c(0,4.23),axes=FALSE,ylab="",xlab="",space=.42)
  text(tv+.09,ys,sprintf("%.3f",tv),adj=0,col=ink,cex=.88)
  axis(1,at=0:4,cex.axis=.86,lwd=.60,col=axis_col)
  axis(2,at=ys,labels=c("Without the Tasmania fold","Without Tasmanian pairs",
    "All held-out pairs"),tick=FALSE,line=-.15,cex.axis=.82)
  mtext("Pooled ΔRMSE (×10⁻³)",1,line=1.51,cex=.93)
}

stem <- file.path(out,"Figure5_palaeoclimate_final")
cairo_pdf(paste0(stem,".pdf"),width=page[1]/72,height=page[2]/72,pointsize=9,bg="white")
draw(); invisible(dev.off())
png(paste0(stem,".png"),width=page[1]/72,height=page[2]/72,units="in",
  res=600,type="cairo",pointsize=9,bg="white")
draw(); invisible(dev.off())
cat(sprintf("Equal-fold RMSE: %.9f vs %.9f; gain %.9f (%.2f%%); P = %.3f\n",
  rmse[1],rmse[2],obs$observed_Delta_RMSE,100*obs$observed_Delta_RMSE/rmse[1],obs$p_value))
cat("101 localities; counts by gain sign: "); print(table(sign(z)))
