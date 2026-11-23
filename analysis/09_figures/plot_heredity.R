#### Vector figures from frozen and matched-CV tables ####
suppressPackageStartupMessages({library(sf);library(ape)})
a <- commandArgs(TRUE); p <- a[1]; out <- a[2]; local <- a[3]
t <- file.path(out,"tables")
read <- function(n) read.csv(file.path(t,n))
ink <- "#251462"; violet <- "#6811DB"; gold <- "#946212"
original_palette <- read.csv(file.path(local,"results/Fig2_Nuclear_structure_revision/P4_palette.csv"))
palette <- setNames(original_palette$colour,original_palette$role)[paste0("K",1:4)]
coast <- st_read(file.path(local,"results/Lampropholis_Fig1_Geographic_context/data/abs_states_generalized.geojson"),quiet=TRUE)
coast <- st_geometry(st_transform(st_make_valid(coast),3577))
xy <- function(lon,lat) st_coordinates(st_transform(st_as_sf(data.frame(lon,lat),coords=c("lon","lat"),crs=4326),3577))
aa <- subset(read("random_fold_assignments.csv"),repetition==1)
points_map <- xy(aa$X,aa$Y)
bx <- range(points_map[,1])+c(-125000,130000); by <- range(points_map[,2])+c(-100000,110000)
setup <- function() par(family="sans",fg=ink,col.lab=ink,col.axis=ink,mgp=c(1.65,.44,0),tcl=-.22,las=1)
map <- function() {
  plot(NA,xlim=bx,ylim=by,asp=1,axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  plot(coast,add=TRUE,col="#F6F4F2",border="#B6AFC0",lwd=.32)
}
save <- function(stem,w,h,draw) {
  cairo_pdf(paste0(stem,".pdf"),w,h,pointsize=8.8,bg="white")
  draw(); invisible(dev.off())
  png(paste0(stem,".png"),w,h,units="in",res=600,type="cairo",pointsize=8.8,bg="white")
  draw(); invisible(dev.off())
}
national <- read("random_cv_summary.csv")
regional <- subset(read("regional_cv_summary.csv"),metric=="delta_gain")
region_names <- c("VIC, TAS & far SE NSW","Central–southern NSW","Northern NSW & SE QLD","Northern QLD","South Australia")
draw_cv <- function() {
  layout(rbind(c(1,2),c(3,3),c(4,4)),heights=c(1.32,1,1.05))
  setup()
  pch <- c(21,22,23,24,25); gr <- c("#FAFAFA","#DFDEE3","#BEBCC5","#9D99A9","#777281")
  for(k in 1:2) {
    par(mar=c(1.75,.7,2.2,.55))
    map()
    f <- if(k==1) aa$fold else aa$random_fold
    points(points_map,pch=pch[f],bg=gr[f],col="#332D40",lwd=.45,cex=.99)
    mtext(c("A   Original spatial folds","Random locality folds · repeat 1")[k],3,.60,adj=0,font=2)
    legend("bottomleft",inset=c(.01,-.09),xpd=NA,horiz=TRUE,legend=1:5,pch=pch,pt.bg=gr,
      col="#332D40",title="Fold",cex=.82,bty="n")
  }
  par(mar=c(3.02,11,2.25,.85))
  plot(c(.070,.137),c(.4,5.4),type="n",axes=FALSE,xlab="",ylab="")
  axis(1); axis(2,5:1,c("Spatial · original pairs","Random · all pairs","Spatial · matched pairs","Random · matched pairs",""),tick=FALSE,cex.axis=.85)
  mtext("B   National held-out prediction",3,.72,adj=0,font=2)
  mtext("RMSE (lower is better)",1,1.66)
  rows <- data.frame(target=c("all_original_spatial_pairs","all_test_pairs","matched_test_pairs","matched_test_pairs"),
    aggregation=c("equal_original_region","equal_random_fold","equal_original_region","equal_original_region"),
    design=c("spatial","random","spatial","random"))
  for(i in 1:4) for(j in 1:2) {
    z <- subset(national,target==rows$target[i] & aggregation==rows$aggregation[i] &
       design==rows$design[i] & metric==c("rmse_m2","rmse_m3")[j])
    y <- 6-i+c(.15,-.15)[j]
    col <- c(ink,gold)[j]; shape <- c(21,23)[j]
    segments(z$q025,y,z$q975,y,col=col,lwd=1.4)
    points(z$median,y,pch=shape,bg=col,col=col,cex=1.17,lwd=.55)
  }
  legend("bottomright",inset=c(.01,-.015),legend=c("M2 · contemporary","M3 · + PaleoMean"),pch=c(21,23),
         pt.bg=c(ink,gold),col=c(ink,gold),horiz=TRUE,bty="n",cex=.81)
  par(mar=c(3.75,11,2.05,1.1))
  plot(c(-.0105,.0065),c(.55,5.50),type="n",axes=FALSE,xlab="",ylab="")
  axis(1);axis(2,5:1,paste(1:5,region_names,sep=" · "),tick=FALSE,cex.axis=.79)
  abline(v=0,lty=2,col="#91879B",lwd=.8)
  regional <- regional[match(1:5,regional$region),]
  segments(regional$q025,5:1,regional$q975,5:1,col=violet,lwd=1.5)
  points(regional$median,5:1,pch=21,bg=violet,col="white",cex=1.2)
  text(.0064,5:1,paste0("n=",regional$n),cex=.73,adj=1,col=ink)
  mtext("C   Palaeoclimate gain on identical test pairs",3,.7,adj=0,font=2)
  mtext("Random − spatial gain in RMSE",1,1.61)
  mtext("Negative: smaller gain under random validation",1,2.69,cex=.79)
}
save(file.path(out,"supplementary/Figure_S_Random_vs_Spatial_CrossValidation"),6.9,7.9,draw_cv)

#### Main discussion synthesis: observed, predictive, and historical evidence ####
s5 <- file.path(p,"results/landscape/gdm_definitive")
s6 <- file.path(p,"results/landscape/paleoclimate")
cs <- read.csv(file.path(s5,"validation/spatial_group_null_summary.csv"))
hs <- read.csv(file.path(s6,"null_corrected/corrected_historical_null_summary.csv"))
fr <- read.csv(file.path(s6,"null_corrected/corrected_observed_fold_effects.csv"))
sc <- read.csv(file.path(p,"results/demography/moments/tables/stage8b_model_comparison.csv"))
tas <- read.csv(file.path(local,"results/Lampropholis_demographic_chronology/figures/Fig3_sample_tip_key.csv"))
mat <- as.matrix(tas[,names(palette)])
dom <- max.col(mat,ties.method="first")
eligible <- apply(mat,1,max)>=.70
gx <- xy(tas$lon,tas$lat)
gain <- c(cs$observed_mean_rmse_gain,hs$observed_Delta_RMSE)
pv <- c(cs$p_value,hs$p_value)
draw_main <- function() {
  layout(rbind(c(1,2),c(1,3)),widths=c(.92,1.16),heights=c(1.19,1))
  setup()
  par(mar=c(2.5,.5,2.3,.5)); map()
  points(gx[!eligible,,drop=FALSE],pch=21,bg="#C9C4D0",col="white",cex=.83,lwd=.35)
  points(gx[eligible,,drop=FALSE],pch=21,bg=palette[dom[eligible]],col="white",cex=.97,lwd=.38)
  mtext("A   Nuclear structure",3,.8,adj=0,font=2)
  legend("bottom",inset=c(0,.02),xpd=NA,pch=21,pt.bg=c(palette,"#C9C4D0"),col="#675E77",
    legend=c(names(palette),"Mixed Q"),horiz=TRUE,bty="n",cex=.84)
  par(mar=c(2.7,7.4,2.1,.55))
  plot(c(0,.010),c(.6,4.25),type="n",axes=FALSE,xlab="",ylab="")
  axis(1,at=c(0,.003,.006,.009),labels=c("0",".003",".006",".009"))
  axis(2,3:1,c("Current climate","Resistance","Palaeoclimate"),tick=FALSE,cex.axis=.85)
  abline(v=0,col="#B2AABF"); segments(0,3:1,gain,3:1,lwd=1,col="#B2AABF")
  points(gain,3:1,pch=c(21,22,23),bg=c("white",gold,violet),col=c(ink,gold,violet),cex=1.23)
  text(gain+.0004,3:1,sprintf("P = %.3f",pv),adj=0,cex=.81)
  text(0,4,"Geographic distance: dominant component",adj=0,cex=.82,col=ink)
  mtext("B   Spatially tested associations",3,.7,adj=0,font=2)
  mtext("Conditional gain in RMSE",1,1.7)
  par(mar=c(.45,.55,1.8,.55)); plot.new(); plot.window(c(0,1),c(0,1))
  mtext("C   Historical interpretation",3,.5,adj=0,font=2)
  n_sc <- sum(sc$model=="SC" & sc$rank==1)
  separator <- c(.69,.42)
  segments(.01,separator,.99,separator,col="#E3DDE9",lwd=.8)
  text(.01,c(.89,.62,.34),c("Palaeoclimate improvement was heterogeneous",
      "Isolation and renewed gene flow are plausible","Taxonomy and timing remain unresolved"),
      adj=0,font=2,cex=.91,col=ink)
  text(.01,c(.79,.51,.23,.12),c(sprintf("Positive in %d of 5 spatial folds",sum(fr$observed_Delta_RMSE>1e-10)),
      sprintf("Baseline secondary contact ranked first in %d contrasts",n_sc),
      "No validated nuclear node ages", "Bunya specimen identities require confirmation"),
      adj=0,cex=.83,col=ink)
}
save(file.path(out,"main/Figure_Heredity_scientific_synthesis"),7.2,4.45,draw_main)

#### Focused taxonomic supplement: original topology and recorded ancestry ####
tax <- file.path(local,"results/Lampropholis_demographic_chronology/qage")
v <- read.csv(file.path(tax,"tables/LAMPROPHOLIS_SAMPLE_SPECIES_CROSSWALK.csv"),colClasses=c(sample_id="character"))
ids <- c("411716",as.character(411828:411834),"411710")
tree <- read.tree(file.path(p,"treemix/svdquartets_ref5/svdquartets_ref5_170_random_final.tre"))
small <- keep.tip(tree,ids); small$edge.length <- NULL
v <- v[match(ids,v$sample_id),]
draw_bunya <- function() {
  layout(matrix(c(1,2,1,3),nrow=2,byrow=TRUE),widths=c(.94,1),heights=c(1.1,1))
  setup()
  par(mar=c(2,1,2.2,.8))
  plot(small,type="phylogram",use.edge.length=FALSE,direction="rightwards",cex=.9,
       tip.color=ifelse(small$tip.label=="411710","#62576F",ink),edge.color=ink,label.offset=.10,no.margin=FALSE)
  mtext("A   Induced nuclear topology",3,.66,adj=0,font=2)
  mtext("Bunya split: 87.4% in the full tree",1,.45,cex=.83)
  par(mar=c(1.3,.5,2.1,.7)); map()
  b <- tas[match(ids,tas$id),]
  pts <- unique(xy(b$lon[1:8],b$lat[1:8]))
  points(pts,pch=21,bg="#EBB201",col=ink,cex=1.12)
  t0 <- xy(b$lon[9],b$lat[9]); points(t0,pch=22,bg="white",col=ink,cex=.9)
  text(pts[1,1]-33000,pts[1,2]+52000,"Bunya (2 sites)",adj=1,cex=.78,col=ink)
  text(t0[1]-35000,t0[2]-16000,"411710 · NSW",adj=1,cex=.78,col=ink)
  mtext("B   Recorded localities",3,.62,adj=0,font=2)
  par(mar=c(2.15,4,2,.45))
  mat <- as.matrix(v[,names(palette)])
  barplot(t(mat[nrow(mat):1,]),horiz=TRUE,col=palette,border=NA,names.arg=rev(ids),las=1,cex.names=.78)
  mtext("C   Global sNMF proportions",3,.62,adj=0,font=2)
  legend("bottom",inset=c(0,-.28),xpd=NA,legend=names(palette),fill=palette,border=NA,horiz=TRUE,bty="n",cex=.80)
}
save(file.path(tax,"figures/Figure_S_Bunya_specimen_review"),6.4,5.45,draw_bunya)
