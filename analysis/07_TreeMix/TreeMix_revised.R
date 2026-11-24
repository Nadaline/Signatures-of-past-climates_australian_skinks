# Figure 6: population relationships and admixture in TreeMix
# Visualize genetic relationships among the 14 geographic groups, with migration edges estimated from allele-frequency covariance. Edges indicate fitted covariance patterns, not observed movements or independently dated migration episodes.
# Preserve the original edge weights and geographic-group definitions; distinguish the plotted graph from evidence for exactly four historical migration events.
#
#### 1. TreeMix relationships among sampling regions ####
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
project <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = TRUE)
plotting <- normalizePath(args[3], mustWork = TRUE)
state_file <- normalizePath(args[4], mustWork = TRUE)
source(plotting)
tmdir <- file.path(project, "treemix")
stem <- file.path(tmdir, "outputs/lampro_m4_r1")
lik <- read.csv(file.path(tmdir, "likelihood_all_runs.csv"))
stopifnot(nrow(lik) == 90L, all(table(lik$m, lik$replicate) == 1L),
  setequal(lik$m, 0:8), setequal(lik$replicate, 1:10), all(is.finite(lik$loglik)),
  lik$loglik[match("lampro_m4_r1.llik", lik$file)] == max(lik$loglik[lik$m == 4]))
curve <- aggregate(loglik ~ m, lik, mean)
vert <- read.table(gzfile(paste0(stem, ".vertices.gz")), as.is = TRUE, quote = "", comment.char = "")
edges <- read.table(gzfile(paste0(stem, ".edges.gz")), as.is = TRUE, quote = "", comment.char = "")
newick <- readLines(gzfile(paste0(stem, ".treeout.gz")))
tip <- vert[vert$V5 == "TIP", c("V1","V2")]
root <- vert$V1[vert$V3 == "ROOT"]
mig <- which(edges$V5 == "MIG")
w <- edges$V4[mig]
stopifnot(nrow(tip) == 14L, length(root) == 1L, length(mig) == 4L, !anyDuplicated(tip$V2),
  identical(w, as.numeric(sub(" .*", "", newick[-1]))))

#### 2. Geographic group definitions and range-wide K = 4 ancestry ####
lab <- c(colossus = "Bunya Mountains", deli1b = "Sarina–Eungella", deli1c = "Rockhampton–Boyne Island",
  deli2 = "Kroombit Tops", deli4a = "Northern NSW ranges", deli4b = "Tweed–Nightcap",
  deli4c = "New England uplands", deli5b = "NSW mid-north ranges", deli7a = "ACT–Southern Tablelands",
  deli7b = "South Australia", deli8 = "Coolah–Liverpool Range", deli9a = "Sydney & surrounds",
  deli9c = "Illawarra–Southern Highlands", deli9d = "Vic.–Tas. & coastal NSW")
pal <- c(K1 = "#6811DB", K2 = "#06ADF1", K3 = "#B8D6E2", K4 = "#EBB201")
ink <- "#251462"
meta <- read.csv(file.path(project, "population_structure/pcoa_scores.csv"), colClasses = c(id = "character"))
site_meta <- read.csv(file.path(project, "Lampropholis_sample_metadata_audit/TableS1_Individual_metadata.csv"),
  colClasses = c(Sample_ID = "character"))[, c("Sample_ID","Analytical_locality","State_region",
  "Locality_description","Latitude","Longitude")]
states <- read.csv(state_file, colClasses = c(id = "character"))[, c("id","state")]
q <- read.csv(file.path(project, "tables/snmf_selected_qmatrix_k4.csv"), colClasses = c(id = "character"))
sizes <- read.csv(file.path(tmdir, "population_sample_sizes.csv"))
stopifnot(nrow(meta) == 170L, !anyDuplicated(meta$id), !anyDuplicated(site_meta$Sample_ID),
  !anyDuplicated(q$id), !anyDuplicated(states$id), setequal(meta$id, q$id),
  setequal(meta$id, site_meta$Sample_ID), setequal(meta$id, states$id),
  setequal(tip$V2, names(lab)), setequal(tip$V2, sizes$population))
header <- strsplit(readLines(gzfile(file.path(tmdir, "populations.curated.treemix.gz")), n = 1L), " ")[[1]]
stopifnot(setequal(header, tip$V2))
meta <- meta[meta$pop %in% tip$V2, ]
site_meta <- site_meta[match(meta$id, site_meta$Sample_ID), ]
stopifnot(all(meta$pop == site_meta$Analytical_locality),
  max(abs(meta$lon-site_meta$Longitude), abs(meta$lat-site_meta$Latitude)) < .00006)
Q <- as.matrix(q[match(meta$id, q$id), names(pal)])
stopifnot(max(abs(rowSums(Q)-1)) < 1e-5, all(Q >= 0 & Q <= 1))
meta$global_top_component <- names(pal)[max.col(Q, ties.method = "first")]
meta$state_from_coordinates <- states$state[match(meta$id, states$id)]
meta$state_as_recorded <- site_meta$State_region
meta$state_for_summary <- ifelse(meta$state_as_recorded == "ACT", "ACT", meta$state_from_coordinates)
meta$locality_as_recorded <- site_meta$Locality_description
ix <- split(seq_len(nrow(meta)), meta$pop)
stopifnot(nrow(meta) == 163L, all(lengths(ix)[sizes$population] == sizes$n))
means <- t(vapply(ix, function(i) colMeans(Q[i,,drop=FALSE]), numeric(4)))
stopifnot(max(abs(rowSums(means)-1)) < 1e-5)
q_summary <- data.frame(population = names(ix), geographic_label = unname(lab[names(ix)]),
  n = unname(lengths(ix)), states = vapply(ix, function(i)
    paste(sort(unique(meta$state_for_summary[i])), collapse = "; "), character(1)),
  means, row.names = NULL)
for (k in names(pal)) q_summary[[paste0("n_top_",k)]] <- vapply(ix, function(i)
  sum(meta$global_top_component[i] == k), integer(1))
write.csv(q_summary, file.path(out,"TreeMix_population_Q_summary.csv"), row.names=FALSE)
crosswalk <- data.frame(population = meta$pop, geographic_label = unname(lab[meta$pop]),
  population_n = unname(lengths(ix)[meta$pop]), sample_id = meta$id,
  locality_as_recorded = meta$locality_as_recorded, longitude = meta$lon, latitude = meta$lat,
  state_as_recorded = meta$state_as_recorded, state_from_coordinates = meta$state_from_coordinates,
  state_for_summary = meta$state_for_summary, global_top_component = meta$global_top_component, Q)
write.csv(crosswalk[order(crosswalk$population,crosswalk$sample_id),],
  file.path(out,"TreeMix_population_label_crosswalk.csv"), row.names=FALSE, na="")

#### 3. Original graph coordinates, unchanged drift scale and directed edges ####
# The previously saved TreeMix renderer applies the edge-weight factor twice to its
# display coordinates. Reproduce that convention; previously saved edge data stay intact.
e_plot <- edges
e_plot$V3 <- edges$V3 * edges$V4 * edges$V4
vert$x <- vert$y <- vert$ymin <- vert$ymax <- NA_real_
vert <- set_y_coords(vert)
vert <- set_x_coords(vert, e_plot)
invisible(capture.output(vert <- set_mig_coords(vert, e_plot)))
from <- match(edges$V1,vert$V1); to <- match(edges$V2,vert$V1)
stopifnot(all(is.finite(vert$x)), all(is.finite(vert$y)), !anyNA(from), !anyNA(to),
  vert$x[match(root,vert$V1)] == 0)
tv <- which(vert$V5 == "TIP")
line_col <- colorRampPalette(c("#BD927A","#AE6549","#84392B"))(256)
weight_range <- c(0.10,0.42)
stopifnot(all(w >= weight_range[1] & w <= weight_range[2]))
mcolour <- function(x) line_col[1L + round((x-weight_range[1])/diff(weight_range)*255)]
print(data.frame(from=edges$V1[mig],to=edges$V2[mig],weight=w,
  xfrom=vert$x[from[mig]],xto=vert$x[to[mig]],yfrom=vert$y[from[mig]],yto=vert$y[to[mig]]),
  digits=10,row.names=FALSE)

#### 4. Vector TreeMix plate ####
draw <- function() {
  par(family="sans",fg=ink,col.axis=ink,col.lab=ink,col.main=ink,las=1,tcl=-.24,mgp=c(1.52,.43,0))
  par(fig=c(0,1,.742,1),mar=c(3.12,5.0,2.75,2.35))
  plot(curve$m,curve$loglik,type="o",pch=19,cex=.59,lwd=1.17,col=ink,bty="l",
    axes=FALSE,xlab="Migration edges (m)",ylab="TreeMix log-likelihood",
    xlim=c(-.2,8.12),ylim=c(480,818))
  axis(1,at=0:8,cex.axis=.89); axis(2,at=c(500,600,700,800),cex.axis=.89)
  box(bty="l",col="#77708A",lwd=.65)
  points(4,curve$loglik[curve$m==4],pch=21,bg="#B45346",col="white",cex=1.25,lwd=.65)
  text(4.15,curve$loglik[curve$m==4]-21,"m = 4",adj=0,col="#84392B",cex=.91)
  mtext("A  TreeMix log-likelihood",side=3,at=-1.22,adj=0,line=1.15,font=2,cex=1.13)

  par(fig=c(0,1,.072,.743),new=TRUE,mar=c(3,1.95,2.95,1.08),xaxs="i",yaxs="i")
  plot.new()
  plot.window(xlim=c(-.001,.116),ylim=c(-.01,1.10))
  mtext("B  Population graph",side=3,adj=0,line=.97,font=2,cex=1.13)
  tx <- c(0,.02,.04,.06)
  axis(1,at=tx,labels=formatC(tx,digits=2,format="f"),cex.axis=.89,lwd=.65)
  mtext("Drift parameter",side=1,at=.03,line=1.65,cex=.94)
  guides <- vert$y[tv]
  segments(vert$x[tv],guides,.056,guides,col="#DDD9E3",lwd=.6,lty=3)
  a <- which(edges$V5=="NOT_MIG")
  segments(vert$x[from[a]],vert$y[from[a]],vert$x[to[a]],vert$y[to[a]],col=ink,lwd=.97)
  t <- seq(0,1,length.out=65)
  # Ease the arcs away from one another without changing either endpoint.
  bend <- c(-.0060,-.0006,.0039,.0049)
  for (j in seq_along(mig)) {
    a <- from[mig[j]]; b <- to[mig[j]]
    xx <- (1-t)*vert$x[a]+t*vert$x[b]+bend[j]*4*t*(1-t)
    yy <- (1-t)*vert$y[a]+t*vert$y[b]
    lines(xx,yy,col="white",lwd=2.60)
    lines(xx,yy,col=mcolour(w[j]),lwd=1.56)
    arrows(xx[62],yy[62],xx[65],yy[65],length=.067,angle=25,col=mcolour(w[j]),lwd=1.35)
  }
  points(vert$x[tv],guides,pch=19,cex=.40,col=ink)
  text(.071,guides,lab[vert$V2[tv]],adj=0,cex=.93,col=ink)
  bx <- .057
  text(bx,1.056,"Mean population Q (global sNMF, K = 4)",adj=0,cex=.85,font=2,col=ink)
  for (i in tv) {
    q_i <- means[vert$V2[i],]
    q_i <- q_i/sum(q_i)
    ends <- bx+c(0,cumsum(q_i))*.0120
    rect(ends[-5],vert$y[i]-.011,ends[-1],vert$y[i]+.011,col=pal,border=NA)
    rect(bx,vert$y[i]-.011,bx+.0120,vert$y[i]+.011,border="#8B819F",lwd=.38)
  }

  par(fig=c(0,1,0,.079),new=TRUE,mar=rep(0,4),xaxs="i",yaxs="i")
  plot.new(); plot.window(c(0,450),c(0,39))
  text(18,25,"Migration weight",adj=0,font=2,col=ink,cex=.92)
  cb <- seq(18,135,length.out=257)
  rect(cb[-257],14,cb[-1],20,col=line_col,border=NA)
  ticks <- c(.10,.20,.30,.40)
  ct <- 18+(ticks-weight_range[1])/diff(weight_range)*117
  segments(ct,12,ct,14,col=ink,lwd=.57)
  text(ct,7,formatC(ticks,digits=2,format="f"),col=ink,cex=.77)
  xx <- c(208,250,292,334)
  rect(xx,14,xx+8,21,col=pal,border=NA)
  text(xx+12,17.5,names(pal),adj=0,col=ink,cex=.89)
}
cairo_pdf(file.path(out,"TreeMix_revised.pdf"),width=450/72,height=498/72,pointsize=9,bg="white")
draw(); invisible(dev.off())
png(file.path(out,"TreeMix_revised.png"),width=450/72,height=498/72,units="in",
  type="cairo",res=600,pointsize=9,bg="white")
draw(); invisible(dev.off())
