#### Rebuild candidate A from archived summary values, not model fitting
args <- commandArgs(TRUE)
stopifnot(length(args)==1,dir.exists(args[1]))
pair <- c("K1_K3","K1_K4","K2_K4")
labels <- gsub("_","–",pair)
pal <- read.csv("design/idv/palette/chapter2_palette.csv")
idv <- setNames(pal$hex,pal$key)
K <- c(K1=idv[["violet"]],K2="#06ADF1",K3="#B8D6E2",K4=idv[["amber"]])
ink <- idv[["ink"]]
z <- read.csv("results/summary_tables/moments/Figure7_A_model_and_SFS_cells.csv")
i <- z$kind=="A_model"
mod <- data.frame(pair=z$id[i],model=z$id2[i],delta_AIC=z$value[i])
arrays <- setNames(vector("list",length(pair)),pair)
for(p in pair){
  x <- z[z$kind=="A_SFS_cell" & z$id==p,]
  n <- max(x$index1)+1L; m <- max(x$index2)+1L
  stopifnot(nrow(x)==n*m,all(x$index1==rep(0:(n-1L),m)),
    all(x$index2==rep(0:(m-1L),each=n)),x$masked[1],
    max(abs(x$observed-x$expected-x$value))<1e-8)
  arrays[[p]] <- list(res=matrix(x$value,n,m),mask=matrix(x$masked,n,m))
}
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
out <- file.path(args[1],"Figure_7_baseline_demography")
cairo_pdf(paste0(out,".pdf"),width=174/25.4,height=128/25.4,family="sans",pointsize=10)
draw_A()
dev.off()
png(paste0(out,".png"),width=174,height=128,units="mm",res=600,type="cairo",pointsize=10)
draw_A()
dev.off()
