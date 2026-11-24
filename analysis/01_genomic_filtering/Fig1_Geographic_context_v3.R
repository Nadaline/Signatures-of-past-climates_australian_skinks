# 1. Load sampling records, boundaries and vegetation data ----
library(sf)
library(terra)
library(jsonlite)
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 4:5)
sample_file <- normalizePath(args[1])
previous_data <- normalizePath(args[2])
out <- normalizePath(args[3])
v2 <- normalizePath(args[4], mustWork = TRUE)
v2_preview <- file.path(v2, "Fig1_Geographic_context_v2_preview_2_ilustração.png")
if (!file.exists(v2_preview)) v2_preview <- file.path(v2, "Fig1_Geographic_context_v2_preview_2.png")
if (length(args) == 5L) v2_preview <- normalizePath(args[5], mustWork = TRUE)
data_dir <- file.path(v2, "data")
# Read individual identifiers and geographic coordinates.
sample <- read.csv(sample_file)[, c("id", "lon", "lat")]
stopifnot(nrow(sample) == 170L, !anyDuplicated(sample$id),
  all(is.finite(sample$lon)), all(is.finite(sample$lat)))
loc <- aggregate(id ~ lon + lat, sample, length)
names(loc)[3] <- "n"
stopifnot(nrow(loc) == 101L, sum(loc$n) == 170L, max(loc$n) == 6L)
loc$size <- cut(loc$n, c(0, 1, 3, 6), labels = c("1", "2–3", "4–6"))
point_size <- c(.72, .94, 1.14)
crs_map <- 9473
# Project sampling coordinates into the map coordinate system.
xy <- function(lon, lat) unname(st_coordinates(st_transform(
  st_as_sf(data.frame(lon, lat), coords = c("lon", "lat"), crs = 4326), crs_map)))
point_xy <- xy(loc$lon, loc$lat)
states <- st_read(file.path(previous_data, "abs_states_generalized.geojson"), quiet = TRUE)
states <- states[states$state_code_2016 %in% as.character(1:8), ]
states <- st_make_valid(st_transform(states, crs_map))
states <- st_cast(suppressWarnings(st_collection_extract(states, "POLYGON")), "MULTIPOLYGON")
land <- vect(states)
dem <- rast(file.path(previous_data, "ga_east_dem.tif"))
panel <- list(A = c(1270000, 2245000, -3360000, -1960000),
              B = c(285000, 2145000, -4960000, -3100000))
placement <- list(A = c(.012, .405, .202, .926), B = c(.423, .989, .202, .926))
contains <- function(box) point_xy[, 1] >= box[1] & point_xy[, 1] <= box[2] &
  point_xy[, 2] >= box[3] & point_xy[, 2] <= box[4]
loc$panel_A <- contains(panel$A)
loc$panel_B <- contains(panel$B)
stopifnot(all(loc$panel_A | loc$panel_B), any(loc$panel_A & loc$panel_B))
write.csv(loc, file.path(out, "Fig1_panel_sampling_audit.csv"), row.names = FALSE)

# Group NVIS vegetation classes by vegetation type.
nvmeta <- fromJSON(file.path(data_dir, "nvis_ext_mvs_metadata.json"), simplifyVector = FALSE)
stopifnot(identical(nvmeta$metadata[[which(names(nvmeta$metadata) == "")]]$RASTER_DATASET,
  "NVIS7_0_AUST_EXT_MVS_ALB"))
lookup <- do.call(rbind, lapply(nvmeta$rat$row, function(r)
  data.frame(code = as.integer(r$f[[1]]), official_name = r$f[[3]])))
category <- list("Rainforest" = c(1L, 2L, 6L),
  "Wet eucalypt forest / woodland" = c(3L, 59L, 60L),
  "Dry forest / woodland" = c(13L, 27L, 62L, 66L),
  "Arid / semi-arid vegetation" = c(20L, 22L, 23L, 31L, 33L, 45L, 51L, 52L))
ids <- as.integer(unlist(category, use.names = FALSE))
stopifnot(!anyDuplicated(ids), all(ids %in% lookup$code))
classes <- data.frame(code = ids, display_group = rep(names(category), lengths(category)))
classes$official_name <- lookup$official_name[match(classes$code, lookup$code)]
stopifnot(all(grepl("rainforest", classes$official_name[classes$code %in% category[[1]]], ignore.case = TRUE)),
  grepl("wet sclerophyll", classes$official_name[classes$code == 3L], fixed = TRUE))
write.csv(classes, file.path(out, "Fig1_NVIS_MVS_classes.csv"), row.names = FALSE)

# 2. Prepare elevation and vegetation layers ----
structure_col <- "#251462"
ocean <- "#F1F6F8"
height_colour <- c("#F4EFE4", "#E1E3D7", "#B7C2B6")
veg_colour <- c("#24734B", "#68A98B", "#C6B18C", "#E8D6B9")
veg_alpha <- c(.93, .88, .76, .81)
gap_col <- "#754326"
mt_col <- "#255348"
ink <- "#393849"

terrain_panel <- function(k) {
  suffix <- if (k == "A") "north" else "south"
  map <- rast(file.path(data_dir, paste0("nvis_ext_mvs_", suffix, ".tif")))
  mask_nvis <- rast(file.path(data_dir, paste0("nvis_ext_mvs_", suffix, "_valid.tif")))
  base <- rast(map)
  stopifnot(same.crs(map, paste0("EPSG:", crs_map)),
    compareGeom(map, mask_nvis, stopOnError = FALSE), all(res(base) == 1800))
  z <- mask(project(dem, base, method = "bilinear"), land)
  sa <- terrain(z, v = c("slope", "aspect"), unit = "radians")
  sh <- terra::shade(sa[["slope"]], sa[["aspect"]], angle = 42, direction = 315)
  v <- values(map, mat = FALSE)
  v[values(mask_nvis, mat = FALSE) != 255] <- NA
  e <- values(z, mat = FALSE)
  h <- values(sh, mat = FALSE)
  landmask <- is.finite(values(rasterize(land, base, field = 1), mat = FALSE))
  stopifnot(all(unique(v[is.finite(v)]) %in% c(0L, 255L, lookup$code)), sum(is.finite(e)) > 1000L)
  cat <- match(classes$display_group[match(v, classes$code)], names(category))
  cat[!landmask] <- NA_integer_
  stopifnot(all(c(1L, 2L) %in% cat))
  list(grid = base, h = h, e = e, cat = cat, land = landmask)
}
terrain_layers <- setNames(lapply(names(panel), terrain_panel), names(panel))
write.csv(do.call(rbind, lapply(names(panel), function(k) {
  t <- table(factor(terrain_layers[[k]]$cat, levels = 1:4))
  data.frame(panel = k, group = names(category), display_cells = as.integer(t))
})), file.path(out, "Fig1_NVIS_display_coverage.csv"), row.names = FALSE)

make_background <- function(x) {
  valid <- x$land
  e <- x$e
  e[!is.finite(e)] <- 0
  bands <- cut(e[valid], c(-Inf, 250, 650, Inf), labels = FALSE, right = FALSE)
  c0 <- t(col2rgb(height_colour[bands])) / 255
  h <- x$h[valid]
  h[!is.finite(h)] <- .68
  strength <- .12
  c0 <- pmin(pmax(c0 * (1 + strength * (h - .65)), 0), 1)
  group <- x$cat[valid]
  mapped <- !is.na(group)
  if (any(mapped)) {
    weight <- veg_alpha[group[mapped]]
    cover <- t(col2rgb(veg_colour[group[mapped]])) / 255
    c0[mapped, ] <- c0[mapped, , drop = FALSE] * (1 - weight) + cover * weight
  }
  all <- rep("#FFFFFF00", length(e))
  all[valid] <- rgb(c0[, 1], c0[, 2], c0[, 3])
  as.raster(matrix(all, nrow = nrow(x$grid), byrow = TRUE))
}

# 3. Draw geographic features and sampling localities ----
gaps <- rbind(Burdekin = c(147.32, -19.77), StLawrence = c(149.75, -22.34),
  Hunter = c(151.03, -32.56))
uplands <- xy(c(149.06, 150.95, 150.07), c(-23.73, -24.42, -31.75))
upland_number <- c(2L, 1L, 3L)
colnames(uplands) <- c("x", "y")
qpos <- xy(152.35, -28.22)
eastg <- xy(148.3, -37.51)
highlands <- xy(148.7, -35.55)
# Place the basin label near Mildura.
murray <- xy(142.16, -34.18)

label <- function(x, y, txt, col = ink, cex = .82, adj = .5, font = 1, srt = 0) {
  u <- par("usr")
  pin <- par("pin")
  dx <- diff(u[1:2]) / (72 * pin[1]) * .51
  dy <- diff(u[3:4]) / (72 * pin[2]) * .51
  for (a in seq(0, 7) * pi / 4)
    text(x + dx * cos(a), y + dy * sin(a), txt, cex = cex, col = "#FFFCF5",
      adj = adj, srt = srt, font = font)
  text(x, y, txt, cex = cex, col = col, adj = adj, srt = srt, font = font)
}
draw_gap <- function(i, tx, ty, txt, adj = 0) {
  x <- xy(gaps[i, 1], gaps[i, 2])
  segments(x[1] - 32000, x[2] - 6000, x[1] + 31000, x[2] + 15000,
    lty = 2, lwd = 1.4, col = gap_col)
  label(tx, ty, txt, col = gap_col, adj = adj, font = 2, cex = .83)
}
draw_isolate <- function(i, shiftx, shifty) {
  p <- uplands[i, ]
  c <- p + c(shiftx, shifty)
  segments(p[1], p[2], c[1], c[2], lwd = .48, col = mt_col)
  points(c[1], c[2], pch = 24, bg = "#96BFAC", col = mt_col, cex = .85, lwd = .70)
  label(c[1] - 25000, c[2], upland_number[i], col = mt_col, cex = .86, font = 2, adj = 1)
}
draw_scale <- function(b, width) {
  x <- b[1] + .06 * diff(b[1:2])
  y <- b[3] + .07 * diff(b[3:4])
  d <- .014 * diff(b[3:4])
  rect(x, y, x + width / 2, y + d, col = structure_col, border = structure_col, lwd = .55)
  rect(x + width / 2, y, x + width, y + d, col = "white", border = structure_col, lwd = .55)
  text(x, y - 1.25*d, "0", cex = .75, col = ink)
  text(x + width, y - 1.25*d, paste(width / 1000, "km"), cex = .75, col = ink, adj = 1)
}
draw_panel <- function(k) {
  par(fig = placement[[k]], new = k != "A", mar = rep(0, 4), xaxs = "i", yaxs = "i")
  b <- panel[[k]]
  x <- terrain_layers[[k]]
  plot.new()
  plot.window(b[1:2], b[3:4], asp = 1)
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col = ocean, border = NA)
  plot(st_geometry(states), col = height_colour[1], border = NA, add = TRUE)
  b_r <- as.vector(ext(x$grid))
  rasterImage(make_background(x), b_r[1], b_r[3], b_r[2], b_r[4], interpolate = FALSE)
  plot(st_geometry(states), col = NA, border = "#7D817F", lwd = .48, add = TRUE)
  box(col = "#8F8B99", lwd = .55)
  if (k == "A") {
    label(1390000, -2905000, "QLD", cex = .92, font = 2)
    draw_gap(1, 1654000, -2173000, "Burdekin Gap")
    draw_gap(2, 1875000, -2488000, "St Lawrence\nGap")
    label(1545000, -2026000, "Wet Tropics", col = mt_col, adj = 1, font = 3)
    label(1645000, -3050000, "Great Dividing Range", col = mt_col, srt = 64, font = 3)
    label(1920000, -3288000, "McPherson /\nMain Range", col = mt_col, adj = 1, cex = .81)
    segments(qpos[1], qpos[2], 1925000, -3271000, col = mt_col, lwd = .48)
  } else {
    label(503000, -3630000, "SA", font = 2)
    label(1110000, -3620000, "NSW", font = 2)
    label(918000, -4170000, "VIC", font = 2)
    label(957000, -4640000, "TAS", font = 2)
    label(murray[1], murray[2], "Murray Basin", font = 3)
    draw_gap(3, 1710000, -3675000, "Hunter Valley", adj = 1)
    label(1380000, -3959000, "Southern NSW\nhighlands", col = mt_col, adj = 1, cex = .79)
    segments(highlands[1], highlands[2], 1390000, -3968500, col = mt_col, lwd = .48)
    label(1518000, -4301000, "East Gippsland", cex = .83, adj = 0)
    segments(eastg[1], eastg[2], 1518000, -4281000, col = "#606876", lwd = .48)
    label(1360000, -4451000, "Bass Strait", font = 3)
  }
  shown <- which(if (k == "A") loc$panel_A else loc$panel_B)
  shown <- shown[order(loc$n[shown], decreasing = TRUE)]
  points(point_xy[shown, , drop = FALSE], pch = 21, cex = point_size[as.integer(loc$size[shown])],
    col = "white", bg = structure_col, lwd = .48)
  if (k == "A") {
    draw_isolate(1, -42000, -15000); draw_isolate(2, -39000, -25000)
  } else draw_isolate(3, -58000, 18000)
  draw_scale(b, if (k == "A") 150000 else 250000)
  u
}

# Assemble the map panels and inset.
draw_figure <- function() {
  par(family = "sans", oma = rep(0, 4))
  a <- draw_panel("A")
  b <- draw_panel("B")
  par(fig = c(0, 1, .93, 1), new = TRUE, mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(.014, .46, "A  Northern eastern Australia", adj = 0, col = structure_col, cex = 1.01, font = 2)
  text(.424, .46, "B  South-eastern Australia", adj = 0, col = structure_col, cex = 1.01, font = 2)
  par(fig = c(.012, .207, .017, .161), new = TRUE, mar = rep(0, 4))
  plot.new()
  bb <- st_bbox(states)
  plot.window(c(bb$xmin, bb$xmax) + c(-1, 1)*300000,
    c(bb$ymin, bb$ymax) + c(-210000, 320000), asp = 1)
  plot(st_geometry(states), col = "#E6E1EA", border = "#898491", lwd = .35, add = TRUE)
  rect(a[1], a[3], a[2], a[4], border = structure_col, lwd = .94,
    col = adjustcolor("#C6B7DF", .24))
  rect(b[1], b[3], b[2], b[4], border = structure_col, lwd = .94, lty = 2,
    col = adjustcolor("#A69DB5", .14))
  text(mean(a[1:2]), mean(a[3:4]), "A", cex = .78, font = 2, col = structure_col)
  text(mean(b[1:2]), mean(b[3:4]), "B", cex = .78, font = 2, col = structure_col)
  par(fig = c(.012, .207, .156, .193), new = TRUE, mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(.50, .48, sprintf("%d localities shown in both", sum(loc$panel_A & loc$panel_B)),
    cex = .70, col = structure_col)

  par(fig = c(.224, .989, .007, .186), new = TRUE, mar = rep(0, 4))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(0, .94, "Sampling locality", adj = 0, cex = .83, font = 2, col = structure_col)
  u <- c(.030, .136, .257)
  points(u, rep(.75, 3), pch = 21, cex = point_size, bg = structure_col, col = "white")
  text(u + .026, .75, levels(loc$size), adj = 0, cex = .79, col = ink)
  text(0, .60, "Individuals per point", adj = 0, cex = .72, col = ink)
  segments(.005, .41, .062, .41, lty = 2, lwd = 1.25, col = gap_col)
  text(.079, .41, "Dry habitat gap", adj = 0, cex = .79, col = ink)
  points(.03, .24, pch = 24, bg = "#96BFAC", col = mt_col, cex = .86)
  text(.079, .24, "Isolated upland habitat", adj = 0, cex = .79, col = ink)
  text(0, .066, "1 Kroombit Tops  ·  2 Blackdown Tableland  ·  3 Coolah Tops", adj = 0, cex = .73, col = ink)
  text(.43, .94, "Vegetation (NVIS)", adj = 0, cex = .83, font = 2, col = structure_col)
  ys <- c(.75, .58, .41, .24)
  blend <- veg_alpha
  swatch <- t(col2rgb(veg_colour)) / 255 * blend + (1 - blend) * rep(.93, 4)
  rect(.43, ys-.057, .468, ys+.053, col = rgb(swatch[, 1], swatch[, 2], swatch[, 3]), border = NA)
  text(.480, ys, names(category), adj = 0, cex = .78, col = ink)
  text(.827, .94, "Elevation (m)", adj = 0, cex = .83, font = 2, col = structure_col)
  rect(.827, ys[1:3]-.057, .865, ys[1:3]+.053, col = height_colour, border = "#B7B7AD", lwd = .35)
  text(.880, ys[1:3], c("<250", "250–649", "≥650"), adj = 0, cex = .78, col = ink)
}

# 4. Export the sampling map ----
stem <- file.path(out, "Fig1_Geographic_context_v3")
proof <- paste0(stem, "_print_preview.png")
png(proof, width = 450/72, height = 350/72, units = "in",
  res = 144, type = "cairo", pointsize = 8.7, bg = "white")
draw_figure(); invisible(dev.off())
cairo_pdf(paste0(stem, ".pdf"), width = 450/72, height = 350/72, pointsize = 8.7, bg = "white")
draw_figure(); invisible(dev.off())
png(paste0(stem, ".png"), width = 450/72, height = 350/72, units = "in",
  res = 600, type = "cairo", pointsize = 8.7, bg = "white")
draw_figure(); invisible(dev.off())
cairo_pdf(paste0(stem, "_comparison.pdf"), width = 2*450/72, height = 365/72, pointsize = 8.7, bg = "white")
par(mar = rep(0, 4))
plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
comparison <- c(v2_preview, proof)
for (i in seq_along(comparison)) {
  preview <- rast(comparison[i])
  pix <- values(preview, mat = TRUE)[, 1:3, drop = FALSE] / 255
  image <- matrix(rgb(pix[, 1], pix[, 2], pix[, 3]), nrow = nrow(preview), byrow = TRUE)
# terra stores unreferenced PNG rows in reverse display order.
  image <- as.raster(image[nrow(image):1, , drop = FALSE])
  rasterImage(image, (i-1)/2, 0, i/2, 350/365, interpolate = FALSE)
}
text(c(.015, .515), .983, c("Figure 1 v2", "Figure 1 v3"),
  adj = 0, col = structure_col, font = 2, cex = 1.15)
invisible(dev.off())
