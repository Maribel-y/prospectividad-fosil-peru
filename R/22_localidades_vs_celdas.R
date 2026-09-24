# ============================================================
# Reconciliation of the two counts of the inventory:
#   666 independent localities  (complete-linkage clusters, 2 km)
#   710 presence cells          (occupied cells of the 2 km grid)
# Both come from the same occurrences but answer different questions,
# so neither is an error. The script traces every step from one to the
# other and writes the table used in the supplementary material.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"

occ  <- fread(file.path(RUTA, "datos/procesados/ocurrencias_final.csv"))       # 11 140
occc <- fread(file.path(RUTA, "datos/procesados/ocurrencias_con_celda.csv"))   # with valid cell
g    <- rast(file.path(RUTA, "datos/procesados/grilla_2km.tif"))

cluster <- function(dt, h) {
  u <- unique(dt[, .(lon, lat)])
  xy <- crds(project(vect(u, geom = c("lon","lat"), crs = "EPSG:4326"), PROJ))
  u[, loc := cutree(hclust(dist(xy), method = "complete"), h = h)]
  merge(dt, u, by = c("lon","lat"))
}

# ---- 1. Localities over the full inventory ----------------------
A <- cluster(occ, 2000)
n_loc <- uniqueN(A$loc)

# ---- 2. Localities over the occurrences that received a cell ----
B <- cluster(occc, 2000)
n_loc_celda <- uniqueN(B$loc)

# ---- 3. How localities map onto cells ----------------------------
# A cluster has a maximum diameter of 2 km, so it can straddle up to
# four cells of 2 km; conversely, two clusters can share one cell.
lc <- unique(B[, .(loc, celda)])
por_loc  <- lc[, .(n_celdas = .N), by = loc]
por_celda <- lc[, .(n_loc = .N), by = celda]

n_celdas <- uniqueN(occc$celda)
tab <- data.table(
  paso = c("Occurrences in the inventory",
           "Unique coordinates",
           "Independent localities (complete linkage, 2 km)",
           "Occurrences with a valid land cell (after snapping)",
           "Independent localities among those occurrences",
           "Presence cells (2 km grid)",
           "Localities split across 2 cells",
           "Localities split across 3-4 cells",
           "Cells shared by 2 or more localities",
           "Net difference cells - localities"),
  valor = c(nrow(occ), nrow(unique(occ[, .(lon, lat)])), n_loc,
            nrow(occc), n_loc_celda, n_celdas,
            sum(por_loc$n_celdas == 2), sum(por_loc$n_celdas >= 3),
            sum(por_celda$n_loc >= 2), n_celdas - n_loc_celda))
print(tab)

# Identity: cells = localities + extra cells from split localities - cells merged
extra <- sum(por_loc$n_celdas - 1); fusion <- sum(por_celda$n_loc - 1)
cat("\nCheck:", n_loc_celda, "+", extra, "(split) -", fusion, "(shared) =",
    n_loc_celda + extra - fusion, "| presence cells:", n_celdas, "\n")

# Same reconciliation for the Mesozoic stratum
Bm <- cluster(occc[era == "Mesozoico"], 2000)
cat("Mesozoic: localities", uniqueN(Bm$loc), "| cells", uniqueN(Bm$celda), "\n")

# Coordinates with degree precision (integer lon and lat), which the text
# says were excluded: quantify how many remain
deg <- occ[lon == round(lon) & lat == round(lat)]
cat("Occurrences with integer-degree coordinates:", nrow(deg),
    "| unique:", nrow(unique(deg[, .(lon, lat)])), "\n")

fwrite(rbind(tab, data.table(paso = c("Split cells added", "Cells merged",
                                      "Mesozoic localities", "Mesozoic presence cells",
                                      "Occurrences at integer-degree coordinates"),
                             valor = c(extra, fusion, uniqueN(Bm$loc), uniqueN(Bm$celda), nrow(deg)))),
       file.path(RUTA, "salidas/tablas/tabla_localidades_vs_celdas.csv"))
