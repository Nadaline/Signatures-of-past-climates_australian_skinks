#### 1. Final SVDQuartets tree, sampling and nested ancestry ####
library(ape)
library(sf)
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L || (length(args)==4L && args[4] %in% c("--bootstrap","--discussion")))
bootstrap_view <- length(args) == 4L
discussion_view <- bootstrap_view && args[4]=="--discussion"
project <- normalizePath(args[1], mustWork = TRUE)
boundary <- normalizePath(args[2], mustWork = TRUE)
out <- normalizePath(args[3], mustWork = TRUE)
tree_file <- file.path(project, "treemix/svdquartets_ref5/svdquartets_ref5_170_random_final.tre")
tr <- read.tree(tree_file)
n <- length(tr$tip.label)
stopifnot(inherits(tr, "phylo"), n == 170L, !anyDuplicated(tr$tip.label),
  length(tr$edge.length) == nrow(tr$edge), all(tr$edge.length >= 0 & tr$edge.length <= 100))
read_ids <- function(path) read.csv(file.path(project, path), colClasses = c(id = "character"))
meta <- read_ids("population_structure/pcoa_scores.csv")
hier <- read_ids("results/structure/hierarchical_snmf/hierarchical_snmf_assignments.csv")
q <- read_ids("tables/snmf_selected_qmatrix_k4.csv")
stopifnot(!anyDuplicated(meta$id), !anyDuplicated(hier$id), !anyDuplicated(q$id),
  setequal(meta$id, tr$tip.label), setequal(q$id, tr$tip.label), all(hier$id %in% tr$tip.label))
meta <- meta[match(tr$tip.label, meta$id), c("id", "pop", "lon", "lat")]
gl <- readRDS(file.path(project, "data/derived/lampro_landscape_background.rds"))
ll <- as.data.frame(gl@other$latlon)[match(meta$id, gl@ind.names), ]
stopifnot(setequal(gl@ind.names, meta$id), identical(as.numeric(meta$lon), as.numeric(ll$lon)),
  identical(as.numeric(meta$lat), as.numeric(ll$lat)))
rm(gl, ll)
macro <- c(K1 = "#6811DB", K2 = "#06ADF1", K3 = "#B8D6E2", K4 = "#EBB201")
nested <- c(K1 = macro[["K1"]], K2 = macro[["K2"]], K3a = macro[["K3"]],
  K3b = "#7FA6B5", K3c = "#416C80", K4a = macro[["K4"]], K4b = "#AA7006")
Q <- as.matrix(q[match(meta$id, q$id), names(macro)])
stopifnot(max(abs(rowSums(Q) - 1)) < 1e-5, all(Q >= 0 & Q <= 1))
meta$nested_assignment <- hier$child_component[match(meta$id, hier$id)]
meta$display_group <- sub("^(K[12])a$", "\\1", meta$nested_assignment)
q_max <- unname(apply(Q, 1, max))
paired <- q_max < .70
stopifnot(setequal(na.omit(meta$display_group), names(nested)),
  identical(paired, is.na(meta$nested_assignment)), sum(paired) == 15L)
ii <- which(paired)
top <- t(apply(Q[ii, , drop = FALSE], 1, function(z) order(-z, seq_along(z))[1:2]))
halves <- t(apply(top, 1, sort))
meta$display_group[ii] <- paste(names(macro)[halves[, 1]], names(macro)[halves[, 2]], sep = "–")
meta$fill_left <- unname(nested[meta$display_group])
meta$fill_right <- meta$fill_left
meta$fill_left[ii] <- unname(macro[halves[, 1]])
meta$fill_right[ii] <- unname(macro[halves[, 2]])
pair_groups <- sort(unique(meta$display_group[ii]))
k3 <- meta$display_group %in% c("K3a", "K3b", "K3c")
k3_link <- c(K3a = "#578195", K3b = "#507488", K3c = "#426B81")
link_fill <- meta$fill_left
link_fill[k3] <- unname(k3_link[meta$display_group[k3]])
link_fill[paired] <- "#7E758B"
stopifnot(!anyNA(meta$fill_left), !anyNA(meta$fill_right))

#### 2. East-coast map from the existing basemap and exact coordinates ####
crs_map <- 3577
pts <- st_transform(st_as_sf(meta, coords = c("lon", "lat"), crs = 4326, remove = FALSE), crs_map)
xy <- st_coordinates(pts)
states <- st_read(boundary, quiet = TRUE)
states <- states[states$state_code_2016 %in% as.character(c(1, 2, 3, 4, 6)), ]
states <- st_cast(suppressWarnings(st_collection_extract(
  st_make_valid(st_transform(states, crs_map)), "POLYGON")), "MULTIPOLYGON")
state_names <- setNames(c("NSW", "VIC", "QLD", "SA", "TAS"), c(1, 2, 3, 4, 6))
meta$state <- unname(state_names[as.character(states$state_code_2016[st_nearest_feature(pts, states)])])
stopifnot(setequal(unique(meta$state), unname(state_names)))
box <- c(xmin = min(xy[, 1]) - 145000, ymin = min(xy[, 2]) - 125000,
  xmax = max(xy[, 1]) + 115000, ymax = max(xy[, 2]) + 115000)
clip <- st_as_sfc(st_bbox(box, crs = st_crs(crs_map)))
land <- suppressWarnings(st_intersection(st_geometry(states), clip))
coast <- suppressWarnings(st_intersection(st_boundary(st_geometry(states)), clip))
site <- match(paste(meta$lon, meta$lat), unique(paste(meta$lon, meta$lat)))
sites <- unname(split(seq_len(n), site))
stopifnot(length(sites) == 101L, all(lengths(sites) >= 1))
north_k3 <- which(meta$state == "QLD" & meta$display_group == "K3a")
north_k3 <- north_k3[which.max(meta$lat[north_k3])]
stopifnot(length(north_k3) == 1L, meta$id[north_k3] == "411871",
  length(sites[[site[north_k3]]]) == 1L)

#### 3. Rotate child blocks to minimize crossings without changing the tree ####
root <- setdiff(tr$edge[, 1], tr$edge[, 2])
stopifnot(length(root) == 1L)
children <- split(tr$edge[, 2], tr$edge[, 1])
nt <- n + tr$Nnode
desc <- vector("list", nt)
depth <- integer(nt)
walk <- function(node, kids, level = 0L) {
  depth[node] <<- level
  v <- if (node <= n) node else unlist(lapply(kids[[as.character(node)]], walk,
    kids = kids, level = level + 1L), use.names = FALSE)
  desc[[node]] <<- v
  v
}
original <- walk(root, children)
target <- -xy[, 2]
cross <- function(order) sum(outer(target[order], target[order], ">")[upper.tri(matrix(FALSE, n, n))])
ladder <- ape::ladderize(tr)
ladder_order <- walk(root, split(ladder$edge[, 2], ladder$edge[, 1]))
original <- walk(root, children)
opt <- children
internal <- order(depth, decreasing = TRUE)
internal <- internal[internal > n]
for (node in internal) {
  k <- as.character(node)
  kid <- children[[k]]
  m <- length(kid)
  stopifnot(m >= 2L, m <= 12L)
  # Between-child crossing cost depends on their order, not internal rotations.
  cost <- matrix(0, m, m)
  for (i in seq_len(m)) for (j in setdiff(seq_len(m), i))
    cost[i, j] <- sum(outer(target[desc[[kid[i]]]], target[desc[[kid[j]]]], ">"))
  # Subset dynamic programming: append the final child to the best preceding subset.
  bits <- bitwShiftL(1L, seq_len(m) - 1L)
  best <- c(0, rep(Inf, 2^m - 1))
  orders <- vector("list", 2^m)
  for (s in seq_len(2^m - 1)) {
    present <- which(bitwAnd(bits, s) != 0L)
    for (j in rev(present)) {
      prev <- s - bits[j]
      delta <- best[prev + 1L] + sum(cost[setdiff(present, j), j])
      if (delta < best[s + 1L]) {
        best[s + 1L] <- delta
        orders[[s + 1L]] <- c(orders[[prev + 1L]], j)
      }
    }
  }
  opt[[k]] <- kid[orders[[2^m]]]
}
tip_order <- walk(root, opt)
stopifnot(setequal(tip_order, seq_len(n)), cross(tip_order) <= min(cross(original), cross(ladder_order)))
rotation_count <- sum(!vapply(names(children), function(k) identical(children[[k]], opt[[k]]), logical(1)))
# Reorder existing edge rows with their support values; do not reroot or edit a branch.
edge_key <- paste(tr$edge[, 1], tr$edge[, 2])
edge_walk <- function(node) {
  kid <- opt[[as.character(node)]]
  if (!length(kid)) return(integer())
  unlist(lapply(kid, function(ch) c(match(paste(node, ch), edge_key), edge_walk(ch))), use.names = FALSE)
}
ei <- edge_walk(root)
rot <- tr
rot$edge <- tr$edge[ei, , drop = FALSE]
rot$edge.length <- tr$edge.length[ei]
attr(rot, "order") <- "cladewise"
stopifnot(!anyNA(ei), identical(sort(ei), seq_len(nrow(tr$edge))),
  identical(rot$edge[order(ei), , drop = FALSE], tr$edge),
  identical(rot$edge.length[order(ei)], tr$edge.length),
  identical(rot$tip.label, tr$tip.label), identical(rot$node.label, tr$node.label),
  identical(is.rooted(rot), is.rooted(tr)))
tas <- which(meta$state == "TAS")
metrics <- data.frame(orientation = c("Source Newick", "Source Newick flipped",
  "Ladderized", "Ladderized flipped", "Crossing-minimized"),
  crossings = vapply(list(original, rev(original), ladder_order, rev(ladder_order), tip_order), cross, numeric(1)),
  Tasmania_median_tip_from_top = vapply(list(original, rev(original), ladder_order, rev(ladder_order), tip_order),
    function(o) median(match(tas, o)), numeric(1)))
write.csv(metrics, file.path(out, "Fig3_orientation_audit.csv"), row.names = FALSE)

#### 4. Tree and map on one vector canvas ####
page <- if (discussion_view) c(660,495) else c(450, if (bootstrap_view) 622 else 495)
if (bootstrap_view) {
  sc_source <- file.path(project,"results/demography/moments/tables/stage8b_final_model_maxima.csv")
  fits <- read.csv(sc_source)
  contrasts <- c("K1_K3","K1_K4","K2_K4")
  sc <- subset(fits,model=="SC" & pair %in% contrasts)
  sc <- sc[match(contrasts,sc$pair),]
  sizes <- do.call(rbind,lapply(strsplit(sc$final_best_parameters,";",fixed=TRUE),as.numeric))
  stopifnot(identical(sc$pair,contrasts),ncol(sizes)==5L,all(is.finite(sizes)),all(sizes>0))
  windows <- data.frame(pair=contrasts,nu1=sizes[,1],nu2=sizes[,2],
    T_iso=sizes[,3],T_sc=sizes[,4],contact_start=sizes[,3]/rowSums(sizes[,3:4]),
    isolation_time_at_upper_bound=sizes[,3]>=10-1e-8,source=basename(sc_source))
  write.csv(windows,file.path(out,"Fig3_pairwise_SC_size_profiles.csv"),row.names=FALSE)
}
ink <- "#251462"
box_x <- c(242, 441)
box_y <- c(41, 467)
sc <- min(diff(box_x)/(box["xmax"] - box["xmin"]), diff(box_y)/(box["ymax"] - box["ymin"]))
tx <- function(x) mean(box_x) + (x - mean(box[c("xmin", "xmax")])) * sc
ty <- function(y) mean(box_y) + (y - mean(box[c("ymin", "ymax")])) * sc
map_x <- tx(xy[, 1]); map_y <- ty(xy[, 2])
stopifnot(all(map_x > box_x[1]), all(map_x < box_x[2]))
tip_x <- 214
entry_x <- box_x[1]
# Topological height only: every actual terminal edge reaches the same X.
height <- integer(nt)
for (node in internal) height[node] <- 1L + max(height[opt[[as.character(node)]]])
node_x <- 12 + (tip_x - 12) * depth / (depth + height)
node_y <- numeric(nt)
node_y[tip_order] <- seq(461, 48, length.out = n)
for (node in internal) node_y[node] <- mean(range(node_y[opt[[as.character(node)]]]))
parent <- tr$edge[, 1]; child <- tr$edge[, 2]
stopifnot(all(node_x[seq_len(n)] == tip_x), all(node_x[child] > node_x[parent]))
support <- tr$edge.length
if (discussion_view) {
  groups <- colnames(Q)[max.col(Q,ties.method="first")]
  split_details <- list()
  split_audit <- do.call(rbind,lapply(windows$pair,function(pair) {
    g <- strsplit(pair,"_",fixed=TRUE)[[1]]
    left <- which(q_max>=.70 & groups==g[1])
    right <- which(q_max>=.70 & groups==g[2])
    stopifnot(length(left)>=2L,length(right)>=2L)
    edge <- which(child>n)
    wrong <- vapply(desc[child[edge]],function(tips) {
      l <- sum(left %in% tips); r <- sum(right %in% tips)
      min(length(left)-l+r,length(right)-r+l)
    },integer(1))
    exact <- edge[wrong==0L]
    top <- exact[support[exact]==max(c(support[exact],-Inf))]
    guide <- if(length(top)) top[which.max(node_x[child[top]]-node_x[parent[top]])] else NA_integer_
    if(length(exact)) split_details[[pair]] <<- data.frame(pair=pair,parent=parent[exact],
      child=child[exact],bootstrap=support[exact],example_marked=exact %in% guide)
    data.frame(pair=pair,n_left=length(left),n_right=length(right),
      min_mismatched_tips=min(wrong),n_compatible_splits=length(exact),
      max_bootstrap_exact=if(length(exact)) max(support[exact]) else NA_real_,
      example_edge=guide)
  }))
  write.csv(split_audit,file.path(out,"Fig3_discussion_split_compatibility.csv"),row.names=FALSE)
  if(length(split_details)) write.csv(do.call(rbind,split_details),
    file.path(out,"Fig3_discussion_compatible_edges.csv"),row.names=FALSE)
}
eligible <- which(child > n & support >= 95 & lengths(desc)[child] >= 12)
support_labels <- integer()
for (j in eligible[order(lengths(desc)[child[eligible]], decreasing = TRUE)]) {
  if (all(abs(node_y[child[j]] - node_y[child[support_labels]]) >= 10))
    support_labels <- c(support_labels, j)
}
if (bootstrap_view) {
  support_labels <- which(child > n & support >= 70)
  support_dots <- which(child > n & support >= 95)
  label_right <- node_x[child[support_labels]] < 29
  label_x <- node_x[child[support_labels]] + ifelse(label_right, 4, -3)
  label_y <- node_y[child[support_labels]] + 4
  # Stagger nearby labels across the incoming split; never move the dot or branch.
  for (i in seq_along(support_labels)) {
    previous <- seq_len(i - 1L)
    collision <- abs(label_x[i] - label_x[previous]) < 23 &
      abs(label_y[i] - label_y[previous]) < 8
    if (any(collision)) label_y[i] <- node_y[child[support_labels[i]]] - 8
  }
}
label_text <- sub("\\.0$", "", format(support[support_labels], trim = TRUE, scientific = FALSE))
tas_order <- match(tas, tip_order)
meta$tip_position_top_to_bottom <- match(seq_len(n), tip_order)
meta$locality_id <- sprintf("L%03d", site)
write.csv(cbind(meta[tip_order, c("tip_position_top_to_bottom", "id", "locality_id", "pop", "state",
  "lon", "lat", "nested_assignment", "display_group", "fill_left", "fill_right")], Q[tip_order, , drop = FALSE]),
  file.path(out, "Fig3_sample_tip_key.csv"), row.names = FALSE, na = "")
pair_table <- data.frame(sample_id = meta$id[ii], locality_id = meta$locality_id[ii],
  tip_position_top_to_bottom = meta$tip_position_top_to_bottom[ii], max_Q = q_max[ii],
  first_component = names(macro)[top[, 1]], first_Q = Q[cbind(ii, top[, 1])],
  second_component = names(macro)[top[, 2]], second_Q = Q[cbind(ii, top[, 2])],
  displayed_pair = meta$display_group[ii],
  plot_representation = ifelse(lengths(sites)[site[ii]] == 1,
    "50:50 split terminal and locality symbols",
    "50:50 split terminal and equal halves within the individual's locality sector"))
write.csv(pair_table[order(pair_table$tip_position_top_to_bottom), ],
  file.path(out, "Fig3_top_two_global_Q.csv"), row.names = FALSE)
write.csv(data.frame(parent_node = parent, child_node = child,
  sample_id = ifelse(child <= n, tr$tip.label[child], NA),
  bootstrap_percent = support, tip_count = lengths(desc)[child],
  label_on_figure = seq_along(child) %in% support_labels),
  file.path(out, "Fig3_branch_support.csv"), row.names = FALSE, na = "")
if (bootstrap_view) {
  write.csv(data.frame(parent_node = parent[child > n], child_node = child[child > n],
    bootstrap_percent = support[child > n], tips_on_descendant_side = lengths(desc)[child[child > n]],
    label_on_figure = which(child > n) %in% support_labels,
    dot_on_figure = which(child > n) %in% support_dots),
    file.path(out, "Fig3_bootstrap_internal_splits.csv"), row.names = FALSE)
}

# Plotting regions are adjusted only by a common linear scale on both map axes.
map_on_page <- function(g) {
  st_crs(g) <- NA
  g * sc + c(mean(box_x) - mean(box[c("xmin", "xmax")]) * sc,
    mean(box_y) - mean(box[c("ymin", "ymax")]) * sc)
}
map_land <- map_on_page(land)
map_coast <- map_on_page(coast)

draw_symbol <- function(x, y, radius, colours, weights = rep(1, length(colours)), stroke = .56) {
  if (length(colours) == 1L) {
    symbols(x, y, circles = radius, inches = FALSE, fg = "white", bg = colours, lwd = stroke, add = TRUE)
  } else {
    a <- pi/2 + c(0, cumsum(weights)/sum(weights)*2*pi)
    for (j in seq_along(colours)) {
      t <- seq(a[j], a[j+1L], length.out = 30)
      polygon(x + c(0, radius*cos(t)), y + c(0, radius*sin(t)), col = colours[j], border = NA)
    }
    symbols(x, y, circles = radius, inches = FALSE, fg = "white", bg = NA, lwd = stroke, add = TRUE)
  }
}
draw_bootstrap_dot <- function(x, y) {
  symbols(x, y, circles = rep(2.05, length(x)), inches = FALSE,
    fg = "white", bg = ink, lwd = .55, add = TRUE)
}
draw_size_profiles <- function() {
  text(12,-12,"C  Relative population sizes under baseline secondary contact",
    adj=0,font=2,col=ink,cex=1.11)
  text(12,-24,"Independent pairwise fits; time is relative to each split, not tied to the tree above.",
    adj=0,col=ink,cex=.83)
  for (i in seq_len(nrow(windows))) {
    a <- 12+(i-1)*147; split <- a+30; present <- a+126
    grp <- strsplit(windows$pair[i],"_",fixed=TRUE)[[1]]
    sizes <- unlist(windows[i,c("nu1","nu2")],use.names=FALSE)
    y <- c(-52,-76); h <- 1.28*sizes
    contact <- split+(present-split)*windows$contact_start[i]
    rect(contact,-88,present,-42,col="#F3EEF8",border=NA)
    rect(a,-64-1.28,split,-64+1.28,col="#C1B8CE",border=NA)
    for (j in 1:2) {
      polygon(c(split,split+6,split+6,split),c(-64+1.28,y[j]+h[j],y[j]-h[j],-64-1.28),
        border=NA,col=adjustcolor(macro[grp[j]],.50))
      rect(split+6,y[j]-h[j],present,y[j]+h[j],border=NA,col=adjustcolor(macro[grp[j]],.50))
    }
    segments(contact,-88,contact,-42,col=ink,lty=2,lwd=.52)
    text(c(a+4,a+4),y,grp,adj=0,col=ink,font=2,cex=.86)
    text(c(present+4,present+4),y,sprintf("%.2f",sizes),adj=0,col=ink,cex=.76)
    text(c(split,present),-96,c("Split","Present"),cex=.72,col=ink)
    if (windows$isolation_time_at_upper_bound[i]) text(contact,-41,"*",col=ink,cex=.92)
  }
  text(12,-111,"Band width: Ne relative to each fitted ancestor (Ne = 1). Dashed: onset of contact.",
    adj=0,col=ink,cex=.79)
  text(12,-121,"* Isolation duration reached the fitting limit for K1–K3.",adj=0,col=ink,cex=.76)
}
draw_discussion_profiles <- function() {
  text(467,485,"Pairwise demographic contrasts",adj=0,font=2,col=ink,cex=1.15)
  text(467,472,"Baseline secondary-contact fits",adj=0,col=ink,cex=.88)
  for (i in seq_len(nrow(windows))) {
    mid <- c(398,267,136)[i]
    pair <- windows[i,]
    audit <- split_audit[i,]
    g <- strsplit(pair$pair,"_",fixed=TRUE)[[1]]
    sizes <- c(pair$nu1,pair$nu2)
    split <- 515; present <- 618
    onset <- split+(present-split)*pair$contact_start
    rect(463,mid-63,651,mid+49,border=NA,col="#F9F7FB")
    for(j in 1:2) {
      x <- c(475,507)[j]
      symbols(x,mid+35,circles=2.8,inches=FALSE,fg="white",bg=macro[g[j]],lwd=.45,add=TRUE)
      text(x+6,mid+35,g[j],adj=0,col=ink,font=2,cex=.97)
    }
    y <- mid+c(7,-17); h <- 1.35*sizes
    rect(onset,mid-28,present,mid+20,border=NA,col="#EFEBF5")
    rect(470,mid-5-1.35,split,mid-5+1.35,border=NA,col="#B8ADCA")
    for(j in 1:2) {
      fill <- adjustcolor(macro[g[j]],.62)
      polygon(c(split,split+6,split+6,split),
        c(mid-5+1.35,y[j]+h[j],y[j]-h[j],mid-5-1.35),border=NA,col=fill)
      rect(split+6,y[j]-h[j],present,y[j]+h[j],border=NA,col=fill)
      text(624,y[j],sprintf("%.2f",sizes[j]),adj=0,col=ink,cex=.86)
    }
    segments(onset,mid-28,onset,mid+20,lty=2,col=ink,lwd=.55)
    if(pair$isolation_time_at_upper_bound) text(onset,mid+24,"*",col=ink)
    text(c(split,present),mid-36,c("Split","Present"),col=ink,cex=.75)
    note <- if(audit$n_compatible_splits) {
      sprintf("%s  ·  %d compatible splits (max. BS %g%%)", LETTERS[i],audit$n_compatible_splits,audit$max_bootstrap_exact)
    } else sprintf("Tree: no split separates both groups")
    text(469,mid-47,note,adj=0,col=ink,cex=.78)
    text(469,mid-56,sprintf("Compared %d + %d tips with ancestry Q ≥0.70",audit$n_left,audit$n_right),
      adj=0,col=ink,cex=.72)
  }
  text(467,c(61,50,39,28,17),c("Width: effective size relative to each ancestor (Ne = 1).",
    "Dashed: onset of contact. * K1–K3 isolation hit a bound.",
    "A–C mark examples, not unique demographic splits.",
    "Time progresses to the right within each fit only.",
    "No model times are assigned to cladogram branches."),adj=0,col=ink,cex=.77)
}
draw <- function() {
  par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", lend = "round")
  plot.new()
  plot.window(c(0, page[1]), c(495-page[2],495))
  title <- if(discussion_view) c("SVDQuartets cladogram","Sampling localities") else
    c("A  SVDQuartets cladogram", "B  Sampling localities")
  text(c(12, 244), if (bootstrap_view) 485 else 482, title,
    adj = 0, col = ink, font = 2, cex = 1.15)
  if (bootstrap_view) {
    draw_bootstrap_dot(15, 472)
    text(20, 472, "Bootstrap ≥95%; labels show ≥70%", adj = 0, col = ink, cex = .88)
  }
  plot(map_land, add = TRUE, col = "#F3F0EC", border = NA)
  plot(map_coast, add = TRUE, col = "#A19CA9", lwd = .48)
  label_coords <- st_coordinates(st_transform(st_as_sf(data.frame(
    lon = c(146.2, 146.4, 143.4, 137.0, 149.3), lat = c(-22.1, -32.5, -36.1, -32.1, -43.65)),
    coords = c("lon", "lat"), crs = 4326), crs_map))
  lab_x <- tx(label_coords[, 1]); lab_y <- ty(label_coords[, 2])
  keep <- lab_x > box_x[1] + 7 & lab_x < box_x[2] - 7 & lab_y > box_y[1]
  text(lab_x[keep], lab_y[keep], c("QLD", "NSW", "VIC", "SA", "Tasmania")[keep],
    col = "#625A75", cex = .76)
  t <- seq(0, 1, length.out = 35)
  # The curves all use the same monotone X parameter up to the map entry;
  # afterwards each one is level. Crossing counts in section 3 are exact
  # pairwise intersections for distinct map latitudes, excluding shared endpoints.
  bx <- (1-t)^3 * tip_x + 3*(1-t)^2*t*(tip_x+9) +
    3*(1-t)*t^2*(entry_x-9) + t^3 * entry_x
  w <- t^2 * (3 - 2*t)
  link_order <- c(tip_order[!k3[tip_order]], setdiff(tip_order[k3[tip_order]], north_k3), north_k3)
  for (i in link_order) {
    x <- c(bx, map_x[i])
    y <- c(node_y[i]*(1-w) + map_y[i]*w, map_y[i])
    if (i == north_k3) lines(x, y, col = "white", lwd = 1.23)
    cl <- adjustcolor(link_fill[i], alpha.f = if (i == north_k3) .95 else if (k3[i]) .52 else .25)
    lines(x, y, col = cl, lwd = if (i == north_k3) .90 else if (k3[i]) .53 else .43)
  }
  segments(node_x[parent], node_y[parent], node_x[parent], node_y[child], col = ink, lwd = .53)
  segments(node_x[parent], node_y[child], node_x[child], node_y[child],
    col = ifelse(child <= n, adjustcolor(ink, alpha.f = .76), ink),
    lwd = ifelse(child <= n, .42, .53))
  # The tip marks sit just beyond the purple branch ends, so the ancestry
  # colours read independently. Their diameter stays below the tip spacing.
  sample_x <- tip_x + if (bootstrap_view) 1.15 else 0
  symbols(rep(sample_x, sum(!paired)), node_y[which(!paired)],
    circles = rep(if (bootstrap_view) 1.15 else 1.12, sum(!paired)),
    inches = FALSE, fg = "white", bg = meta$fill_left[!paired],
    lwd = if (bootstrap_view) .52 else .36, add = TRUE)
  for (i in ii) draw_symbol(sample_x, node_y[i], if (bootstrap_view) 1.15 else 1.20,
    c(meta$fill_left[i], meta$fill_right[i]), stroke = if (bootstrap_view) .52 else .29)
  if (!bootstrap_view) symbols(tip_x, node_y[north_k3], circles = 1.19, inches = FALSE,
    fg = k3_link[["K3a"]], bg = meta$fill_left[north_k3], lwd = .46, add = TRUE)
  if (length(support_labels)) {
    sx <- if (bootstrap_view) label_x else node_x[child[support_labels]] - 3
    sy <- if (bootstrap_view) label_y else node_y[child[support_labels]] + 3.4
    size <- if (bootstrap_view) .84 else .66
    tw <- strwidth(label_text, cex = size, font = if (bootstrap_view) 2 else 1)
    th <- strheight(label_text, cex = size)
    right <- if (bootstrap_view) label_right else rep(FALSE, length(sx))
    rect(sx - ifelse(right, 0, tw) - .6, sy - .45,
      sx + ifelse(right, tw, 0) + .65, sy + th + .6, border = NA, col = "white")
    for (i in seq_along(sx)) text(sx[i], sy[i], label_text[i],
      adj = c(ifelse(right[i], 0, 1), 0), col = ink, cex = size, font = if (bootstrap_view) 2 else 1)
  }
  if (bootstrap_view && length(support_dots)) {
    draw_bootstrap_dot(node_x[child[support_dots]], node_y[child[support_dots]])
  }
  if (discussion_view) {
    e <- na.omit(unique(split_audit$example_edge))
    for(j in e) {
      lx <- max(22,mean(node_x[c(parent[j],child[j])]))
      ly <- node_y[child[j]]-5
      lab <- paste(LETTERS[which(split_audit$example_edge==j)],collapse="/")
      tw <- strwidth(lab,cex=.77,font=2)/2+2
      rect(lx-tw,ly-3,lx+tw,ly+3,border=NA,col="white")
      text(lx,ly,lab,col=ink,font=2,cex=.77)
    }
  }
  for (members in sites) {
    i <- members[1]
    profiles <- table(meta$display_group[members])
    reps <- members[match(names(profiles), meta$display_group[members])]
    radius <- if (bootstrap_view) {
      if (length(members) == 1L) 2.55 else 3.35
    } else if (length(members) == 1L) 2.29 else 3.14
    colours <- character()
    weights <- numeric()
    for (j in seq_along(reps)) {
      r <- reps[j]
      colours <- c(colours, if (paired[r]) c(meta$fill_left[r], meta$fill_right[r]) else meta$fill_left[r])
      weights <- c(weights, if (paired[r]) rep(profiles[j]/2, 2L) else profiles[j])
    }
    draw_symbol(map_x[i], map_y[i], radius, colours, weights, stroke = if (bootstrap_view) .82 else .56)
  }
  if (!bootstrap_view) symbols(map_x[north_k3], map_y[north_k3], circles = 1.94, inches = FALSE,
    fg = k3_link[["K3a"]], bg = NA, lwd = .44, add = TRUE)
  text(12, 29, "Nested sNMF ancestry", adj = 0, font = 2, col = ink, cex = .96)
  key_x <- c(15, 42, 69, 99, 129, 159, 189)
  symbols(key_x, rep(14, 7), circles = rep(2.65, 7), inches = FALSE,
    fg = "white", bg = nested, lwd = .5, add = TRUE)
  text(key_x + 5, 14, names(nested), adj = 0, col = ink, cex = .94)
  text(244, 29, "Top two global Q components", adj = 0, font = 2, col = ink, cex = .96)
  pair_x <- seq(248, by = 56, length.out = length(pair_groups))
  for (j in seq_along(pair_groups)) {
    r <- ii[match(pair_groups[j], meta$display_group[ii])]
    draw_symbol(pair_x[j], 14, 3.05, c(meta$fill_left[r], meta$fill_right[r]))
  }
  text(pair_x + 6, 14, pair_groups, adj = 0, col = ink, cex = .94)
  if (discussion_view) draw_discussion_profiles() else if (bootstrap_view) draw_size_profiles()
}
stem <- file.path(out, if (discussion_view) "Fig3_SVDQuartets_demography_discussion" else
  if (bootstrap_view) "Fig3_SVDQuartets_geography_bootstrap" else "Fig3_SVDQuartets_geography_v4")
cairo_pdf(paste0(stem, ".pdf"), width = page[1]/72, height = page[2]/72, pointsize = 8.5, bg = "white")
draw(); invisible(dev.off())
png(paste0(stem, ".png"), width = page[1]/72, height = page[2]/72, units = "in",
  res = 600, pointsize = 8.5, type = "cairo", bg = "white")
draw(); invisible(dev.off())
cat(sprintf("Tree: %d tips, %d edges, rooted=%s; changed child order at %d nodes.\n",
  n, nrow(tr$edge), is.rooted(tr), rotation_count))
print(metrics, row.names = FALSE)
cat(sprintf("Tasmania: %d tips, ranks %d–%d; %d bootstrap labels printed.\n",
  length(tas), min(tas_order), max(tas_order), length(support_labels)))
print(table(meta$display_group))
