# Figure 2: genetic structure across the sampled range
# Combine observed sampling locations with principal coordinate analysis (PCoA) and sNMF ancestry coefficients. These complementary summaries illustrate genetic variation without treating ancestry components as named populations or species.
# The plots use previously estimated ordination and sNMF results. Check that sample identifiers agree across coordinate, PCoA and ancestry tables before drawing the panels.
#
#### 1. Read sampling coordinates and global and nested ancestry estimates ####
library(sf)
library(terra)
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 3:4)
project <- normalizePath(args[1])
boundary <- normalizePath(args[2])
out <- normalizePath(args[3])
phase <- if (length(args) == 4L) args[4] else "final"
stopifnot(phase %in% c("preview", "final"))
pal <- read.csv(file.path(out, "P4_palette.csv"))
pal <- setNames(pal$colour, pal$role)
K <- paste0("K", 1:4)
coords <- read.csv(file.path(project, "population_structure", "pcoa_scores.csv"),
  colClasses = c(id = "character"))
arch <- readRDS(file.path(project, "..", "Reference", "pcoa.ref7.rds"))
stopifnot(nrow(coords) == 170L, !anyDuplicated(coords$id), nrow(arch$loadings) == 19352L,
  nrow(arch$scores) == 170L, setequal(coords$id, rownames(arch$scores)))
scores <- arch$scores[coords$id, 1:3, drop = FALSE]
stopifnot(min(diag(cor(scores[, 1:2], coords[, c("PCo1", "PCo2")]))) > .999999)
positive <- arch$eig[arch$eig > 0]
percent <- 100 * positive[1:3] / sum(positive)
stopifnot(identical(round(percent[1:2], 1), c(28.2, 10.8)))
qd <- read.csv(file.path(project, "tables", "snmf_selected_qmatrix_k4.csv"),
  colClasses = c(id = "character"))
stopifnot(!anyDuplicated(qd$id), setequal(coords$id, qd$id))
Q <- as.matrix(qd[match(coords$id, qd$id), K])
stopifnot(all(is.finite(Q)), all(Q >= 0 & Q <= 1), max(abs(rowSums(Q) - 1)) < 1e-5)
long <- read.csv(file.path(project, "population_structure", "snmf_k4_individual_ancestry_long.csv"),
  colClasses = c(id = "character"))
ij <- cbind(match(long$id, coords$id), match(long$ancestry, K))
stopifnot(nrow(long) == 680L, !anyNA(ij), !anyDuplicated(ij),
  max(abs(Q[ij] - long$q)) < 1e-8,
  all(long$longitude == coords$lon[ij[, 1]]), all(long$latitude == coords$lat[ij[, 1]]))
dom <- K[max.col(Q, ties.method = "first")]
hier <- file.path(project, "results", "structure", "hierarchical_snmf")
h <- read.csv(file.path(hier, "hierarchical_snmf_assignments.csv"), colClasses = c(id = "character"))
decision <- read.csv(file.path(hier, "hierarchical_snmf_final_hierarchy.csv"))
stopifnot(!anyDuplicated(h$id), identical(decision$selected_child_K[match(K, decision$parent)], c(1L, 1L, 3L, 2L)))
children <- list(K3 = c("K3a", "K3b", "K3c"), K4 = c("K4a", "K4b"))
nested_ids <- lapply(names(children), function(k) {
  ii <- which(dom == k & apply(Q, 1, max) >= .70)
  hh <- h[h$global_component == k & h$selected_child_K == length(children[[k]]), ]
  stopifnot(setequal(coords$id[ii], hh$id))
  z <- as.matrix(hh[, children[[k]]])
  stopifnot(all(is.finite(z)), all(z >= 0 & z <= 1), max(abs(rowSums(z) - 1)) < 1e-5,
    identical(children[[k]][max.col(z, ties.method = "first")], hh$child_component))
  ii
})
names(nested_ids) <- names(children)
stopifnot(identical(lengths(nested_ids), c(K3 = 56L, K4 = 45L)))

#### 2. Geographic ordering is independent of the ancestry coefficients ####
crs_map <- 9473
states <- st_read(boundary, quiet = TRUE)
states <- states[states$state_code_2016 %in% as.character(1:8), ]
states <- st_cast(suppressWarnings(st_collection_extract(
  st_make_valid(st_transform(states, crs_map)), "POLYGON")), "MULTIPOLYGON")
sample_geo <- st_transform(st_as_sf(coords, coords = c("lon", "lat"), crs = 4326), crs_map)
state_hit <- st_intersects(sample_geo, states)
stopifnot(!any(lengths(state_hit) > 1L))
nearest <- st_nearest_feature(sample_geo, states)
assigned <- vapply(state_hit, function(z) if (length(z)) z[1] else NA_integer_, integer(1))
offshore <- is.na(assigned)
gap <- rep(0, nrow(coords))
gap[offshore] <- as.numeric(st_distance(sample_geo[offshore, ], states[nearest[offshore], ], by_element = TRUE))
stopifnot(all(gap < 12000))
assigned[offshore] <- nearest[offshore]
region <- states$state_name_2016[assigned]
region[region == "Australian Capital Territory"] <- "New South Wales"
regions <- names(sort(tapply(coords$lat, region, median), decreasing = TRUE))
geo_order <- order(match(region, regions), -coords$lat, coords$lon, coords$id)
stopifnot(length(geo_order) == 170L)
sites <- aggregate(Q, coords[c("lon", "lat")], mean)
counts <- aggregate(id ~ lon + lat, coords, length)
stopifnot(nrow(sites) == 101L, identical(sites[c("lon", "lat")], counts[c("lon", "lat")]), sum(counts$id) == 170L)
site_xy <- st_coordinates(st_transform(st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), crs_map))

#### 3. check the plotted scores and the K2–K4 comparison ####
pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
separation <- function(pair, confident = FALSE) {
  eligible <- !confident | apply(Q, 1, max) >= .70
  a <- scores[dom == "K2" & eligible, pair, drop = FALSE]
  b <- scores[dom == "K4" & eligible, pair, drop = FALSE]
  cx <- colMeans(a); cy <- colMeans(b)
  d <- sqrt(sum((cx - cy)^2))
  spread <- sqrt((sum(sweep(a, 2, cx)^2) + sum(sweep(b, 2, cy)^2)) / (nrow(a) + nrow(b) - 2))
  data.frame(axes = paste(pair, collapse = "x"), subset = if (confident) "global Q >= 0.70" else "dominant Q",
    n_K2 = nrow(a), n_K4 = nrow(b), centroid_distance = d, pooled_within_RMS = spread,
    distance_per_within_RMS = d / spread)
}
sep <- do.call(rbind, c(lapply(pairs, separation), lapply(pairs, separation, confident = TRUE)))
write.csv(sep, file.path(out, "Fig2_K2_K4_projection_check.csv"), row.names = FALSE)
audit <- data.frame(id = coords$id, lon = coords$lon, lat = coords$lat,
  region = region, nearest_boundary_m = gap, Q, dominant = dom, scores)
write.csv(audit[geo_order, ], file.path(out, "Fig2_individual_order_and_scores.csv"), row.names = FALSE)
write.csv(data.frame(sites, individuals = counts$id), file.path(out, "Fig2_locality_mean_Q.csv"), row.names = FALSE)
print(data.frame(axis = 1:3, percent))
print(sep)
print(table(region))

#### 4. Common typography and layout helpers ####
struct <- pal["structural"]
ink <- "#3C3747"
lim <- lapply(1:3, function(j) range(scores[, j]) + c(-1.1, 1.1))
axis_label <- sprintf("PCoA %d (%.1f%%)", 1:3, percent)
axis_label[3] <- sprintf("PCoA 3 (%.2f%%)", percent[3])
header <- function(rect, txt, first = FALSE) {
  par(fig = rect, new = !first, mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(.009, .50, txt, adj = 0, col = struct, font = 2, cex = 1.02)
}
pcoa <- function(rect, pair, units = 4.15) {
  par(fig = rect, new = TRUE, mar = rep(0, 4))
  plot.new()
  width <- diff(rect[1:2]) * par("din")[1] * 72
  height <- diff(rect[3:4]) * par("din")[2] * 72
  xa <- lim[[pair[1]]]; ya <- lim[[pair[2]]]
  stopifnot(diff(xa)*units < width-33, diff(ya)*units < height-27)
  start_x <- max(28, (width - diff(xa) * units) / 2)
  start_y <- max(26, (height - diff(ya) * units) / 2)
  par(plt = c(start_x, start_x + diff(xa)*units, start_y, start_y + diff(ya)*units) / c(width,width,height,height))
  plot.window(xa, ya, xaxs = "i", yaxs = "i")
  abline(h = 0, v = 0, lty = 3, lwd = .52, col = "#D5CFD9")
  points(scores[, pair], pch = 21, bg = pal[dom], col = "#352A47", cex = .79, lwd = .47)
  axis(1, at = pretty(xa, 4), cex.axis = .88, col = struct, col.ticks = struct, col.axis = ink,
    lwd = .58, mgp = c(2, .34, 0), tck = -.035)
  axis(2, at = pretty(ya, 4), cex.axis = .88, col = struct, col.ticks = struct, col.axis = ink,
    lwd = .58, las = 1, mgp = c(2, .33, 0), tck = -.035)
  mtext(axis_label[pair[1]], side = 1, line = 1.38, col = struct, cex = .92)
  mtext(axis_label[pair[2]], side = 2, line = 1.75, col = struct, cex = .92)
}
global_key <- function(rect, txt = "sNMF (K = 4)") {
  par(fig = rect, new = TRUE, mar = rep(0, 4))
  plot.new(); plot.window(c(0,1), c(0,1), xaxs = "i", yaxs = "i")
  text(.015,.5,txt,adj=0,col=struct,font=2,cex=.92)
  at <- c(.38,.53,.68,.83)
  points(at, rep(.50,4), pch=21,bg=pal[K],col=struct,cex=.95,lwd=.55)
  text(at+.03,.50,K,adj=0,col=ink,cex=.92)
}
bars <- function(rect, mat, label, groups, state_names = FALSE, axis_left = TRUE) {
  par(fig = rect, new = TRUE, mar = c(if (state_names) 2.95 else 1.65,
    if (axis_left) 3.5 else .55, .25, .5), xaxs = "i", yaxs = "i")
  plot.new()
  n <- nrow(mat)
  plot.window(c(.5, n+.5), c(0,1))
  normalized <- mat / rowSums(mat)
  end <- t(apply(normalized, 1, cumsum))
  staircase <- as.vector(rbind((1:n)-.5, (1:n)+.5))
  for (i in seq_along(label)) {
    base <- if (i == 1L) rep(0, n) else end[, i-1L]
    polygon(c(staircase, rev(staircase)),
      c(rep(base, each = 2L), rev(rep(end[, i], each = 2L))), col = pal[label[i]], border = NA)
  }
  rr <- rle(groups)
  cuts <- cumsum(rr$lengths)
  if (length(cuts) > 1L) abline(v = head(cuts, -1L)+.5, col = "#FFFFFF", lwd = .92)
  if (axis_left) {
    axis(2, at = c(0,.5,1), labels = c("0","0.5","1"), las = 1, mgp = c(2,.34,0),
      tck = -.027, col = struct, col.axis = ink, lwd = .6, cex.axis = .86)
    mtext("Ancestry (Q)", side = 2, line = 1.62, cex = .88, col = struct)
  }
  if (state_names) {
    centers <- cuts - (rr$lengths-1)/2
    loc_label <- sub("New South Wales", "New South\nWales", rr$values)
    loc_label <- sub("South Australia", "South\nAustralia", loc_label)
    tier <- rr$values == "Tasmania"
    down <- -.08 - tier * (1.15 * strheight("M", cex = .87) + .025)
    segments(centers[tier], -.015, centers[tier], down[tier] + .036, col = "#9F93A6",
      xpd = NA, lwd = .51)
    text(centers, down, loc_label, xpd = NA, cex = .87, col = ink, adj = c(.5,1))
  } else {
    x <- (seq_along(label)-.58) / length(label) * n
    rect(x-1.4, -.33, x-.2, -.20, col = pal[label], border = NA, xpd = NA)
    text(x, -.266, label, adj = 0, xpd = NA, cex = .90, col = ink)
  }
}

#### 5. Four projected maps and the two proposed compositions ####
map_bounds <- c(260000, 2140000, -4920000, -2010000)
qcolour <- function(group, q) {
  base <- as.numeric(col2rgb(pal[group])) / 255
  white <- c(.992, .99, .998)
  q <- pmin(1, pmax(0, q))
  val <- outer(q, base) + outer(1-q, white)
  rgb(val[,1], val[,2], val[,3])
}
map_panel <- function(rect, group) {
  par(fig = rect, new = TRUE, mar = c(.12,.12,.15,.12))
  plot.new()
  plot.window(map_bounds[1:2], map_bounds[3:4], asp = 1, xaxs = "i", yaxs = "i")
  plot(st_geometry(states), add = TRUE, border = "#A9A1AF", lwd = .36, col = "#F1EEF1")
  points(site_xy, pch = 21, bg = qcolour(group, sites[[group]]), col = "#514260", cex = .60, lwd = .40)
  u <- par("usr")
  text(u[1]+.055*diff(u[1:2]), u[4]-.075*diff(u[3:4]), group, adj = 0, cex = .95, col = struct, font = 2)
}
maps <- function(x, y) {
  xb <- seq(x[1], x[2], length.out=3L)
  yb <- seq(y[1]+.055, y[2], length.out=3L)
  j <- 1L
  for (row in 2:1) for (col in 1:2) {
    map_panel(c(xb[col], xb[col+1], yb[row], yb[row+1]), K[j])
    j <- j+1L
  }
  par(fig = c(x, y[1], y[1]+.055), new=TRUE, mar=rep(0,4))
  plot.new(); plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  text(.013,.76,"Mean ancestry coefficient (Q)", adj=0, cex=.91, col=struct)
  px <- c(.20,.51,.83)
  text(px,.10,c("0","0.5","1"), cex=.85, col=ink)
  for (j in 1:4)
    points(px+(j-2.5)*.034, rep(.37,3), pch=21, bg=qcolour(K[j],c(0,.5,1)), col="#514260",cex=.51,lwd=.44)
}
nested <- function(y0,y1) {
  for (j in 1:2) {
    k <- names(children)[j]
    ii <- geo_order[geo_order %in% nested_ids[[k]]]
    mat <- as.matrix(h[match(coords$id[ii], h$id), children[[k]]])
    xb <- if (j == 1L) c(.025,.494) else c(.515,.978)
    header(c(xb,y1-.036,y1), sprintf("%s  ·  n = %d", k, length(ii)))
    bars(c(xb,y0,y1-.040), mat, children[[k]], region[ii], axis_left=j==1L)
  }
}
draw_full <- function() {
  par(family="sans", oma=rep(0,4))
  header(c(.032,.49,.962,.994),"A  PCoA 1 × 2", TRUE)
  header(c(.522,.978,.962,.994),"B  PCoA 1 × 3")
  header(c(.522,.978,.757,.789),"C  PCoA 2 × 3")
  pcoa(c(.030,.492,.708,.951), pairs[[1]])
  pcoa(c(.526,.980,.787,.954), pairs[[2]])
  pcoa(c(.526,.980,.577,.749), pairs[[3]])
  global_key(c(.068,.465,.671,.712))
  header(c(.027,.980,.550,.580),"D  Global sNMF ancestry")
  bars(c(.025,.980,.413,.549), Q[geo_order,], K, region[geo_order], TRUE)
  header(c(.027,.482,.379,.407),"E  Spatial distribution")
  maps(c(.031,.484), c(.026,.378))
  header(c(.515,.978,.379,.407),"F  Nested sNMF")
  for (j in 1:2) {
    k <- names(children)[j]
    ii <- geo_order[geo_order %in% nested_ids[[k]]]
    ends <- c(.208,.024)[j]
    header(c(.516,.977,ends+.136,ends+.169), sprintf("%s  ·  n = %d", k, length(ii)))
    bars(c(.513,.977,ends,ends+.135), as.matrix(h[match(coords$id[ii],h$id),children[[k]]]),
      children[[k]],region[ii])
  }
}
draw_main <- function() {
  par(family="sans", oma=rep(0,4))
  header(c(.032,.542,.954,.995),"A  Nuclear PCoA", TRUE)
  pcoa(c(.028,.544,.618,.942), pairs[[1]])
  global_key(c(.061,.503,.583,.626))
  header(c(.565,.982,.954,.995),"B  Spatial distribution of ancestry")
  maps(c(.559,.978), c(.604,.953))
  header(c(.027,.978,.546,.581),"C  Global sNMF ancestry")
  bars(c(.023,.980,.359,.548), Q[geo_order,], K, region[geo_order], TRUE)
  header(c(.027,.978,.300,.347),"D  Nested sNMF (global Q ≥ 0.70)")
  nested(.032,.300)
}
draw_supp <- function() {
  par(family="sans", oma=rep(0,4))
  header(c(.024,.479,.911,.982),"A  PCoA axes 1 and 3", TRUE)
  header(c(.514,.982,.911,.982),"B  PCoA axes 2 and 3")
  pcoa(c(.026,.482,.17,.903), pairs[[2]])
  pcoa(c(.512,.982,.17,.903), pairs[[3]])
  global_key(c(.278,.746,.026,.14))
}

#### 6. Two 450 pt proofs; retain the main figure and the two-plane supplement ####
proof <- file.path(out, paste0("Fig2_composition_", 1:2, ".png"))
height <- c(850, 490)
draw <- list(draw_full, draw_main)
for (i in 1:2) {
  png(proof[i], width = 450/72, height = height[i]/72, units = "in", res = 144,
    type = "cairo", pointsize = 9, bg = "white")
  draw[[i]](); invisible(dev.off())
}
if (phase == "final") {
  cairo_pdf(file.path(out, "Fig2_Nuclear_structure.pdf"), width = 450/72, height = height[2]/72,
    pointsize = 9, bg = "white")
  draw_main(); invisible(dev.off())
  png(file.path(out, "Fig2_Nuclear_structure.png"), width = 450/72, height = height[2]/72,
    units = "in", res = 600, type = "cairo", pointsize = 9, bg = "white")
  draw_main(); invisible(dev.off())
  cairo_pdf(file.path(out, "FigS_PCoA_three_axes.pdf"), width = 450/72, height = 200/72,
    pointsize = 9, bg = "white")
  draw_supp(); invisible(dev.off())
  total <- max(height) + 20
  cairo_pdf(file.path(out, "Fig2_composition_comparison.pdf"), width = 900/72, height = total/72,
    pointsize = 9, bg = "white")
  par(mar=rep(0,4)); plot.new(); plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  for (i in 1:2) {
    r <- suppressWarnings(rast(proof[i]))
    z <- values(r, mat=TRUE)[,1:3]/255
    px <- matrix(rgb(z[,1],z[,2],z[,3]), nrow=nrow(r),byrow=TRUE)
    rasterImage(as.raster(px[nrow(px):1,,drop=FALSE]), (i-1)/2,
      (max(height)-height[i])/total, i/2, max(height)/total)
  }
  text(c(.017,.517), .986, c("1  All panels together", "2  Main figure + PCoA supplement"),
    adj=0,col=struct,font=2,cex=1.1)
  invisible(dev.off())
}
