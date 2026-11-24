# Figure 7: demographic-model comparison and frequency-spectrum residuals
# Reconstruct the comparison of strict isolation, continuous migration, ancient migration and secondary contact from previously estimated folded site-frequency spectra and model fits.
# The fitted spectra and AIC rankings are inputs to plotting; this script does not optimize demographic parameters. Residual plots locate model mismatches rather than testing each frequency bin.
#
#### Read fitted spectra and model summaries without refitting
args <- commandArgs(TRUE)
stopifnot(length(args) == 3)
package <- args[1]; project <- args[2]; out <- args[3]
tab <- file.path(package,"03_tables")
plot_in <- file.path(package,"04_reproducibility/plot_inputs")
write_to <- file.path(out,"03_tables/moments")
fig_to <- file.path(out,"01_main_figures")
pair <- c("K1_K3","K1_K4","K2_K4")
labels <- gsub("_","–",pair)
palette <- read.csv(file.path(package,"04_reproducibility/notes/chapter2_palette.csv"))
idv <- setNames(palette$hex,palette$key)
K <- c(K1=idv[["violet"]],K2="#06ADF1",K3="#B8D6E2",K4=idv[["amber"]])
ink <- idv[["ink"]]
mod <- read.csv(file.path(tab,"moments/baseline_moments_model_comparison.csv"))
article <- read.csv(file.path(tab,"moments/Table_moments_demography_article.csv"))
stopifnot(nrow(mod)==12L,all(mod$AIC == 2*mod$k-2*mod$logLik |
                          abs(mod$AIC - (2*mod$k-2*mod$logLik)) < 1e-8))
for(type in c("SI","IM","AM","SC")) {
  z <- mod[mod$model==type,]
  stopifnot(max(abs(article[[paste0("delta_AIC_",type)]]-z$delta_AIC[match(pair,z$pair)]))<1e-8)
}
write.csv(article,file.path(write_to,"Table_moments_demography_article_R6_checked.csv"),row.names=FALSE,na="")

audit <- data.frame()
arrays <- list()
for (p in pair) {
  sources <- file.path(plot_in,"S_moments",paste0(p,"_SC_",c("observed","expected","residual"),".csv"))
  arr <- lapply(sources,function(f) as.matrix(read.csv(f,header=FALSE)))
  projection <- dim(arr[[1]])-1L
  f <- file.path(project,"results/demography/sfs/final",paste0(p,"_projected_folded_2d.fs"))
  fs <- readLines(f,n=3)
  head <- scan(text=fs[1],what="",quiet=TRUE)
  stopifnot(all(as.integer(head[1:2])==dim(arr[[1]])),head[3]=="folded",
            identical(head[4:5],strsplit(p,"_")[[1]]))
  mask <- matrix(as.integer(strsplit(trimws(fs[3])," +")[[1]]) != 0L,
                 nrow=projection[1]+1L,byrow=TRUE)
  mask[1,1] <- TRUE
  x <- scan(text=fs[2],quiet=TRUE)
  obs <- matrix(x,nrow=projection[1]+1L,byrow=TRUE)
  stopifnot(all(dim(obs)==dim(mask)),
            isTRUE(all.equal(obs[!mask],arr[[1]][!mask],tolerance=1e-10)),
            max(abs(arr[[1]]-arr[[2]]-arr[[3]]))<1e-8,
            abs(sum(arr[[1]][!mask])-sum(arr[[2]][!mask]))<1e-7)
  arrays[[p]] <- list(obs=arr[[1]],exp=arr[[2]],res=arr[[3]],mask=mask)
  audit <- rbind(audit,data.frame(contrast=p,projection=paste(projection,collapse="x"),
             n_unmasked=sum(!mask),obs_mass=sum(arr[[1]][!mask]),expected_mass=sum(arr[[2]][!mask]),
             obs_origin=arr[[1]][1,1],expected_origin=arr[[2]][1,1],
             min_residual=min(arr[[3]][!mask]),max_residual=max(arr[[3]][!mask])))
}
write.csv(audit,file.path(write_to,"Figure7_R6_SFS_lineage.csv"),row.names=FALSE)

#### Geographic order and matching membership
div <- read.csv(file.path(tab,"diversity/population_diversity.csv"))
fst <- read.csv(file.path(tab,"diversity/pairwise_fst.csv"),row.names=1)
q <- read.csv(file.path(tab,"treemix/TreeMix_population_Q_summary.csv"))
loc <- read.csv(file.path(package,"04_reproducibility/plot_inputs/03_phylogeny/Fig3_sample_tip_key.csv"))
stopifnot(all(c("pop","lat","lon") %in% names(loc)))
geo <- aggregate(cbind(lat,lon)~pop,loc,median)
geo <- geo[geo$pop %in% div$pop[div$n>=3 & div$pop!="unknown"],]
geo <- geo[order(-geo$lat,geo$lon,geo$pop),]
q <- q[match(geo$pop,q$population),]
fst <- as.matrix(fst[geo$pop,geo$pop])
stopifnot(nrow(geo)==14,nrow(q)==14,sum(q$n)==163,
          identical(q$population,geo$pop),all(q$n==div$n[match(q$population,div$pop)]),
          sum(is.finite(fst[lower.tri(fst)]))==91,isSymmetric(fst,check.attributes=FALSE),
          all(abs(rowSums(q[names(K)])-1)<0.001))

#### Export the exact cells, populations and model rows behind both figures
src <- list()
src[[1]] <- data.frame(kind="A_model",id=mod$pair,id2=mod$model,index1=NA,index2=NA,
                       value=mod$delta_AIC,observed=NA,expected=NA,masked=NA)
for(p in pair) {
  a <- arrays[[p]]
  src[[length(src)+1L]] <- data.frame(kind="A_SFS_cell",id=p,id2="SC",
      index1=as.vector(row(a$res))-1L,index2=as.vector(col(a$res))-1L,
      value=as.vector(a$res),observed=as.vector(a$obs),expected=as.vector(a$exp),masked=as.vector(a$mask))
}
g <- which(lower.tri(fst),arr.ind=TRUE)
src[[length(src)+1L]] <- data.frame(kind="B_population_pair",id=geo$pop[g[,1]],id2=geo$pop[g[,2]],
      index1=g[,1],index2=g[,2],value=fst[g],observed=NA,expected=NA,masked=NA)
for(k in names(K)) src[[length(src)+1L]] <- data.frame(kind="B_mean_Q",id=geo$pop,id2=k,
      index1=seq_len(nrow(q)),index2=q$n,value=q[[k]],observed=NA,expected=NA,masked=NA)
write.csv(do.call(rbind,src),file.path(write_to,"Figure7_R6_source_data.csv"),row.names=FALSE,na="")
write.csv(cbind(geo,q[setdiff(names(q),"population")]),file.path(write_to,"Figure7_R6_population_order.csv"),row.names=FALSE)

#### Figure A: ranking and signed residuals answer different questions
lim <- max(unlist(lapply(arrays,function(a) abs(a$res[!a$mask]))))
res_cols <- colorRampPalette(c(idv[["violet"]],idv[["paper"]],idv[["amber"]]))(257)
stretch <- function(x) asinh(x / .25)
draw_A <- function() {
  layout(matrix(c(1,1,1,2,3,4,5,5,5),nrow=3,byrow=TRUE),heights=c(1.02,1.15,.31))
  par(font=1,cex=.86,fg=ink,col.axis=ink,col.lab=ink,mar=c(3.5,6.8,1.8,1.0),mgp=c(2.05,.52,0),tcl=-.2)
  plot(NA,xlim=log1p(c(0,540)),ylim=c(.45,3.8),type="n",axes=FALSE,xlab="ΔAIC relative to SC [log(1 + ΔAIC) spacing]",ylab="")
  at <- c(0,5,10,25,50,100,250,500)
  axis(1,at=log1p(at),labels=at)
  for(i in 1:3){
    y <- 4-i
    points(log1p(c(0,400)),c(y,y),type="l",col=idv[["rule"]])
    b <- strsplit(pair[i],"_")[[1]]
    usr <- par("usr")
    points(usr[1]-.45,y,pch=21,bg=K[b[1]],col=ink,xpd=NA,cex=.9)
    points(usr[1]-.26,y,pch=21,bg=K[b[2]],col=ink,xpd=NA,cex=.9)
    text(usr[1]-.54,y,labels[i],pos=2,xpd=NA,offset=.1)
    value <- mod$delta_AIC[mod$pair==pair[i] & mod$model=="IM"]
    points(log1p(value),y,pch=21,bg=ink,col=ink,cex=1.25)
    text(log1p(value),y,sprintf("%.2f",value),pos=4,offset=.45,cex=.95)
  }
  mtext("A   Support relative to the best alternative",side=3,adj=0,cex=.87,font=2,line=.4)
  for(i in 1:3) {
    a <- arrays[[pair[i]]]; vals <- a$res; vals[a$mask] <- NA_real_
    lims <- dim(vals)-1
    par(mar=c(3.15,if(i==1)3.4 else 2.7,2.0,.5),mgp=c(1.92,.48,0))
    plot.new()
    plot.window(c(-.5,lims[1]+.5),c(-.5,lims[2]+.5),xaxs="i",yaxs="i")
    rect(-.5,-.5,lims[1]+.5,lims[2]+.5,col="#EAE7ED",border=NA)
    image(0:lims[1],0:lims[2],stretch(vals),col=res_cols,zlim=c(-1,1)*stretch(lim),useRaster=FALSE,add=TRUE)
    axis(1,at=pretty(0:lims[1],4)); axis(2,at=pretty(0:lims[2],4),las=1)
    b <- strsplit(pair[i],"_")[[1]]
    mtext(paste(b[1],"allele count"),1,line=1.65,cex=.73)
    mtext(paste(b[2],"allele count"),2,line=1.65,cex=.73)
    title(main=paste0(if(i==1)"B   " else "",labels[i]),cex.main=.97,line=.48)
    box(col=idv[["rule"]])
  }
  par(mar=c(1.6,4.9,.4,4.9),mgp=c(1,.18,0))
  plot.new(); plot.window(c(-stretch(lim),stretch(lim)),c(0,1),xaxs="i",yaxs="i")
  e <- seq(-stretch(lim),stretch(lim),length.out=258)
  rect(head(e,-1),.60,tail(e,-1),.90,col=res_cols,border=NA)
  ticks <- c(-500,-50,-5,0,5,50,500)
  axis(1,at=stretch(ticks),labels=ticks,pos=.60,cex.axis=.75)
  mtext("Observed − expected SFS mass [asinh spacing]",1,line=.35,cex=.72)
}

#### Figure B: geography sets the order, Q is a separate annotation
draw_B <- function() {
  par(cex=.83,fg=ink,col.axis=ink,col.lab=ink,mar=c(6.3,12.8,1.6,1.1),mgp=c(1.8,.43,0),tcl=-.18)
  plot.new(); plot.window(c(-.55,17.7),c(14.05,-.55),xaxs="i",yaxs="i")
  breaks <- seq(0,max(fst),length.out=257)
  cols <- colorRampPalette(c("#FFF7E2",idv[["amber"]],idv[["violet"]],ink))(256)
  for(i in seq_len(14)) {
    for(j in seq_len(14)) {
      val <- fst[i,j]
      fill <- if(i==j || !is.finite(val)) "#EAE7ED" else if(j>i) "#FFFFFF" else
        cols[findInterval(val,breaks,all.inside=TRUE)]
      rect(j-1.5,i-1.5,j-.5,i-.5,col=fill,border="white",lwd=.6)
      if(i==j) text(j-1,i-1,"–",col=ink,cex=.85)
    }
    left <- 14.1
    for(k in names(K)) {
      right <- left + q[[k]][i]*3.30
      rect(left,i-1.27,right,i-.73,col=K[k],border=NA); left <- right
    }
  }
  label <- paste0(seq_len(14),"  ",q$geographic_label,ifelse(q$population=="colossus","*", ""),"  (",q$n,")")
  text(-.66,0:13,label,adj=1,xpd=NA,cex=.89)
  text(0:13,13.80,seq_len(14),adj=.5,xpd=NA,cex=.95)
  text(15.8,-.45,"Mean ancestry (Q)",font=2,cex=.92,xpd=NA)
  for(k in 1:4){
    rect(14.06+(k-1)*.90,13.64,14.29+(k-1)*.90,13.95,col=K[k],border=ink,xpd=NA)
    text(14.35+(k-1)*.90,13.80,names(K)[k],adj=0,cex=.81,xpd=NA)
  }
  x <- seq(0,9,length.out=257)
  rect(head(x,-1),14.98,tail(x,-1),15.30,col=cols,border=NA,xpd=NA)
  ticks <- c(0,.15,.30,.45,.60)
  text(ticks/max(fst)*9,15.56,sprintf("%.2f",ticks),cex=.86,xpd=NA)
  text(4.5,16.22,"Pairwise nuclear differentiation (FST)",cex=.94,xpd=NA)
}

save <- function(tag,draw,height) {
  f <- file.path(fig_to,paste0("Figure7_candidate_",tag,"_R6"))
  cairo_pdf(paste0(f,".pdf"),width=174/25.4,height=height/25.4,family="sans",pointsize=10)
  draw(); dev.off()
  png(paste0(f,".png"),width=174,height=height,units="mm",res=600,type="cairo",pointsize=10)
  draw(); dev.off()
}
save("A",draw_A,128)
save("B",draw_B,152)
