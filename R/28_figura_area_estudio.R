# ============================================================
# Study-area figure: (a) location of Peru, (b) relief and fossil
# presence cells by era, with the two enlarged basins, (c) Pisco-Ica
# basin, (d) Bagua basin, (e) elevation and distance to roads of the
# presence cells, which shows the sampling bias the paper corrects.
# Output: salidas/figuras/fig_area_estudio_en.png
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(sf); library(ggplot2); library(grid)})
RUTA <- path.expand("~/Desktop/dinos")
sf_use_s2(FALSE)
LL <- "EPSG:4326"

# ---- Data ------------------------------------------------------
peru  <- st_as_sf(readRDS(file.path(RUTA, "datos/crudos/gadm/gadm41_PER_0_pk.rds")) |> unwrap())
deps  <- st_as_sf(readRDS(file.path(RUTA, "datos/crudos/gadm/gadm41_PER_1_pk.rds")) |> unwrap())
world <- st_as_sf(readRDS(file.path(RUTA, "datos/crudos/gadm/gadm36_adm0_r5_pk.rds")) |> unwrap())
elev  <- rast(file.path(RUTA, "datos/crudos/elevation/PER_elv_msk.tif"))

g   <- rast(file.path(RUTA, "datos/procesados/grilla_2km.tif"))
occ <- fread(file.path(RUTA, "datos/procesados/ocurrencias_con_celda.csv"))
occ[, era := fifelse(era %in% c("Cenozoico","Mesozoico","Paleozoico"), era, "Sin edad")]
# Era of a cell: the most frequent era among its dated occurrences
celdas <- occ[, .(n = .N, era = { e <- era[era != "Sin edad"]
                                  if (length(e)) names(which.max(table(e))) else "Sin edad" }), by = celda]
xy <- xyFromCell(g, celdas$celda)
ll <- project(vect(xy, crs = crs(g)), LL)
celdas[, c("lon","lat") := as.data.table(crds(ll))]
tm <- fread(file.path(RUTA, "datos/procesados/tabla_modelado.csv"), select = c("celda","elevacion","dist_vias_km"))
celdas <- merge(celdas, tm, by = "celda")
ERAS <- c(Cenozoico = "Cenozoic", Mesozoico = "Mesozoic", Paleozoico = "Paleozoic or older", "Sin edad" = "No age")
celdas[, Era := factor(ERAS[era], levels = ERAS)]
COL <- setNames(c("#E69F00", "#0072B2", "#CC79A7", "#8C8C8C"), ERAS)
med_terr <- median(tm$dist_vias_km); med_pres <- median(celdas$dist_vias_km)

# Relief: hillshade plus hypsometric tint, aggregated for plotting
el <- aggregate(elev, 3, mean, na.rm = TRUE)
hs <- shade(terrain(el, "slope", unit = "radians"), terrain(el, "aspect", unit = "radians"), 40, 315)
rel <- as.data.frame(c(el, hs), xy = TRUE, na.rm = TRUE); names(rel) <- c("x","y","elev","hs")

rios <- tryCatch({
  f <- file.path(tempdir(), "HydroRIVERS_v10_sa.gdb")
  if (!dir.exists(f)) unzip(file.path(RUTA, "datos/crudos/HydroRIVERS_v10_sa.gdb.zip"), exdir = tempdir())
  r <- st_read(list.files(tempdir(), "HydroRIVERS_v10_sa.gdb$", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1],
               query = "SELECT * FROM HydroRIVERS_v10_sa WHERE UPLAND_SKM >= 300", quiet = TRUE)
  suppressWarnings(st_intersection(st_make_valid(r), st_union(peru)))
}, error = function(e) NULL)

# ---- Helpers -----------------------------------------------------
TEMA <- theme_bw(base_size = 8) + theme(panel.grid = element_blank(), axis.title = element_blank(),
  panel.border = element_rect(colour = "grey25", linewidth = 0.4), plot.margin = margin(2, 2, 2, 2),
  legend.position = "none", panel.background = element_rect(fill = "#F2F1EF"))
barra <- function(x0, y0, km, lat) {           # scale bar in degrees at a given latitude
  d <- km / (111.32 * cos(lat * pi / 180)); s <- d / 4
  list(annotate("rect", xmin = x0 + s * 0:3, xmax = x0 + s * 1:4, ymin = y0, ymax = y0 + d * 0.035,
                fill = rep(c("black","white"), 2), colour = "black", linewidth = 0.25),
       annotate("text", x = x0 + d / 2, y = y0 + d * 0.035, label = paste(km, "km"), vjust = -0.5, size = 2.3))
}
norte <- function(x, y, h) list(
  annotate("polygon", x = c(x, x - h * 0.28, x, x + h * 0.28), y = c(y + h, y, y + h * 0.3, y),
           fill = "black", colour = "black", linewidth = 0.2),
  annotate("text", x = x, y = y + h * 1.25, label = "N", size = 2.5, fontface = "bold"))
caja <- function(b) annotate("rect", xmin = b[1], xmax = b[2], ymin = b[3], ymax = b[4],
                             fill = NA, colour = "#C0142B", linewidth = 0.6)
BC <- c(-76.60, -74.90, -15.60, -13.50)   # Pisco-Ica
BD <- c(-78.90, -77.90, -6.30, -5.00)     # Bagua
puntos <- function(d, escala) geom_point(data = d, aes(lon, lat, colour = Era, size = n), alpha = 0.8)

# ---- (a) Globe -----------------------------------------------------
ORTO <- "+proj=ortho +lat_0=-12 +lon_0=-68"
# Visible hemisphere cut on the sphere (s2) before projecting, so that
# countries on the far side do not produce invalid geometries
sf_use_s2(TRUE)
tapa <- st_buffer(st_sfc(st_point(c(-68, -12)), crs = 4326), 9.6e6)
w_v  <- suppressWarnings(st_intersection(st_make_valid(world), tapa))
sf_use_s2(FALSE)
w_o  <- st_make_valid(st_transform(w_v, ORTO))
circ <- st_buffer(st_sfc(st_point(c(0, 0)), crs = ORTO), 6371000)
pe_o <- st_transform(peru, ORTO)
bb   <- st_as_sfc(st_bbox(pe_o)) |> st_buffer(150000)
pa <- ggplot() + geom_sf(data = circ, fill = "#D6E6F0", colour = "grey30", linewidth = 0.3) +
  geom_sf(data = w_o, fill = "#F4F2EE", colour = "grey60", linewidth = 0.15) +
  geom_sf(data = pe_o, fill = "#E7B08A", colour = "grey30", linewidth = 0.2) +
  geom_sf(data = bb, fill = NA, colour = "#C0142B", linewidth = 0.6) +
  coord_sf(crs = ORTO, datum = NA) + theme_void()

# ---- (b) Peru ------------------------------------------------------
pb <- ggplot() +
  geom_raster(data = rel, aes(x, y, fill = elev), alpha = 1) +
  scale_fill_gradientn(colours = c("#2F7D4B", "#A7C47F", "#E9E1B8", "#C9A27A", "#9A6B4F", "#F4F4F4"),
                       values = scales::rescale(c(0, 500, 1500, 3000, 4500, 6500)), guide = "none") +
  geom_raster(data = rel, aes(x, y, alpha = hs), fill = "black") +
  scale_alpha(range = c(0.35, 0), guide = "none") +
  { if (!is.null(rios)) geom_sf(data = rios[rios$UPLAND_SKM >= 20000, ], colour = "#7FA9C9", linewidth = 0.25) } +
  geom_sf(data = deps, fill = NA, colour = "grey35", linewidth = 0.15, linetype = "dotted") +
  geom_sf(data = peru, fill = NA, colour = "black", linewidth = 0.45) +
  geom_point(data = celdas, aes(lon, lat, colour = Era), size = 0.8, alpha = 0.95, stroke = 0) +
  scale_colour_manual(values = COL) +
  caja(BC) + caja(BD) +
  annotate("text", x = BC[1] - 0.15, y = BC[4] + 0.25, label = "c", fontface = "bold", size = 3) +
  annotate("text", x = BD[1] - 0.15, y = BD[4] + 0.25, label = "d", fontface = "bold", size = 3) +
  annotate("label", x = -73.2, y = -7.2, label = "PERU", size = 2.6, fontface = "bold", linewidth = 0, fill = "white", alpha = 0.7) +
  barra(-80.9, -17.9, 300, -17) + norte(-69.3, -1.6, 0.9) +
  coord_sf(xlim = c(-81.5, -68.5), ylim = c(-18.5, 0.2), expand = FALSE) + TEMA

# ---- (c) Pisco-Ica and (d) Bagua ------------------------------------
zoom <- function(b, lab, bar_km, rio_min = 2000, xbreaks = waiver()) {
  d <- celdas[lon >= b[1] & lon <= b[2] & lat >= b[3] & lat <= b[4]]
  ggplot() + geom_sf(data = peru, fill = "#F4F2EE", colour = "grey40", linewidth = 0.3) +
    geom_sf(data = deps, fill = NA, colour = "grey45", linewidth = 0.2, linetype = "dotted") +
    { if (!is.null(rios)) geom_sf(data = rios[rios$UPLAND_SKM >= rio_min, ], colour = "#9DBCD4", linewidth = 0.3) } +
    puntos(d) + scale_colour_manual(values = COL) +
    scale_size_area(max_size = 8, limits = c(1, max(celdas$n)), breaks = c(1, 10, 100)) +
    annotate("label", x = b[1] + 0.03 * (b[2] - b[1]), y = b[4] - 0.04 * (b[4] - b[3]), label = lab,
             hjust = 0, vjust = 1, size = 2.6, fontface = "bold", linewidth = 0, fill = "white") +
    barra(b[2] - 0.40 * (b[2] - b[1]), b[3] + 0.05 * (b[4] - b[3]), bar_km, mean(b[3:4])) +
    scale_x_continuous(breaks = xbreaks) +
    coord_sf(xlim = b[1:2], ylim = b[3:4], expand = FALSE) + TEMA +
    theme(panel.background = element_rect(fill = "#D6E6F0"), panel.border = element_rect(colour = "#C0142B", linewidth = 0.8))
}
pc <- zoom(BC, "Pisco–Ica basin", 40) + norte(BC[2] - 0.12, BC[4] - 0.35, 0.18)
pd <- zoom(BD, "Bagua basin", 20, rio_min = 500, xbreaks = c(-78.8, -78.4, -78.0))

# ---- Legend ----------------------------------------------------------
leyenda <- cowplot_free <- ggplot(celdas, aes(lon, lat, colour = Era, size = n)) + geom_point() +
  scale_colour_manual(values = COL, name = "Era of the cell") +
  scale_size_area(max_size = 8, limits = c(1, max(celdas$n)), breaks = c(1, 10, 100), name = "Occurrences\nper cell") +
  guides(colour = guide_legend(override.aes = list(size = 3), order = 1)) +
  theme_void(base_size = 8) + theme(legend.position = "right", legend.box = "vertical")
leg <- ggplotGrob(leyenda)
leg <- leg$grobs[grep("guide-box", sapply(leg$grobs, function(x) x$name))]
leg <- leg[!sapply(leg, inherits, "zeroGrob")][[1]]

# ---- (e) Elevation vs distance to roads ------------------------------
pe <- ggplot(celdas, aes(dist_vias_km, elevacion, colour = Era, size = n)) +
  geom_vline(xintercept = med_terr, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = med_pres, linetype = "dashed", colour = "#C0142B", linewidth = 0.35) +
  annotate("text", x = med_terr * 1.08, y = 5250, label = sprintf("territory median\n%.1f km", med_terr), hjust = 0, size = 2.3, colour = "grey30") +
  annotate("text", x = med_pres * 1.08, y = 5250, label = sprintf("presence-cell median\n%.1f km", med_pres), hjust = 0, size = 2.3, colour = "#C0142B") +
  geom_point(alpha = 0.7) + scale_colour_manual(values = COL) +
  scale_size_area(max_size = 8, limits = c(1, max(celdas$n)), guide = "none") +
  scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 5, 10, 25, 50, 100)) +
  scale_y_continuous(limits = c(NA, 5600), breaks = seq(0, 5000, 1000)) +
  labs(x = "Distance from the presence cell to the nearest road (km, square-root scale)", y = "Elevation (m a.s.l.)") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), legend.position = "none",
                                  panel.border = element_rect(colour = "grey25", linewidth = 0.4))

# ---- Compose ---------------------------------------------------------
out <- file.path(RUTA, "salidas/figuras/fig_area_estudio_en.png")
png(out, width = 180, height = 200, units = "mm", res = 300)
grid.newpage()
vp <- function(x, y, w, h) viewport(x = x, y = y, width = w, height = h, just = c("left","bottom"))
etq <- function(l, x, y) grid.text(l, x, y, gp = gpar(fontface = "bold", fontsize = 11), just = c("left","top"))
print(pa, vp = vp(0.00, 0.79, 0.30, 0.20)); etq("a", 0.01, 0.995)
print(pb, vp = vp(0.00, 0.30, 0.45, 0.50)); etq("b", 0.01, 0.80)
print(pc, vp = vp(0.47, 0.60, 0.52, 0.39)); etq("c", 0.46, 0.995)
print(pd, vp = vp(0.47, 0.30, 0.25, 0.29)); etq("d", 0.46, 0.595)
pushViewport(vp(0.74, 0.30, 0.25, 0.29)); grid.draw(leg); popViewport()
print(pe, vp = vp(0.00, 0.00, 1.00, 0.29)); etq("e", 0.01, 0.295)
dev.off()
cat("Figure written:", out, "| presence cells:", nrow(celdas), "| median distance to roads: presences",
    round(med_pres, 1), "km, territory", round(med_terr, 1), "km\n")
