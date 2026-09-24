# ============================================================
# Fase 2 (parte C) — Covariables de accesibilidad
# Modelan el sesgo de muestreo: dónde se ha podido buscar fósiles.
# Fuentes: OSM (vías, centros poblados) e HydroRIVERS (ríos).
# ============================================================
suppressPackageStartupMessages({library(terra); library(geodata)})
RUTA <- path.expand("~/Desktop/dinos")
g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
tierra <- g$tierra

dist_a <- function(v, nombre) {
  v <- project(v, crs(g))
  r <- rasterize(v, tierra, field = 1)
  d <- distance(r)                    # metros al elemento más próximo
  d <- mask(d, tierra) / 1000         # en km
  names(d) <- nombre
  d
}

cat("Descargando red vial OSM...\n")
vias <- osm("PER", "highways", path = file.path(RUTA,"datos/crudos"))
cat("  segmentos:", nrow(vias), "\n")
d_via <- dist_a(vias, "dist_vias_km")

cat("Descargando centros poblados OSM...\n")
pobl <- osm("PER", "places", path = file.path(RUTA,"datos/crudos"))
cat("  lugares:", nrow(pobl), "\n")
d_pob <- dist_a(pobl, "dist_poblados_km")

# Ríos: en la Amazonía son la vía de acceso principal, así que pesan
# en el sesgo tanto como las carreteras.
f_rios <- file.path(RUTA,"datos/crudos/HydroRIVERS_v10_sa.gdb.zip")
if (!file.exists(f_rios)) {
  download.file("https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_sa.gdb.zip",
                f_rios, mode = "wb", quiet = TRUE)
}
rios <- tryCatch({
  v <- vect(paste0("/vsizip/", f_rios))
  v[v$ORD_FLOW <= 5]     # solo cursos de orden mayor (navegables)
}, error = function(e) { cat("  ríos no disponibles:", conditionMessage(e), "\n"); NULL })

capas <- c(d_via, d_pob)
if (!is.null(rios)) {
  cat("  tramos de río:", nrow(rios), "\n")
  capas <- c(capas, dist_a(rios, "dist_rios_km"))
}

writeRaster(capas, file.path(RUTA,"datos/procesados/covar_acceso_2km.tif"), overwrite = TRUE)
print(summary(values(capas), digits = 3))
cat("\nCovariables de accesibilidad escritas.\n")
