# ============================================================
# Fase 2 (parte A) — Covariables topográficas sin Google Earth Engine
# SRTM agregado vía geodata; pendiente, rugosidad y curvatura a 2 km.
# ============================================================
suppressPackageStartupMessages({library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"

g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
dem <- geodata::elevation_30s("PER", path=file.path(RUTA,"datos/crudos"))
cat("DEM origen:", paste(dim(dem), collapse=" x "), "| res:", round(res(dem)[1]*111,2), "km aprox\n")

pend <- terrain(dem, "slope",   unit="degrees")
rug  <- terrain(dem, "TRI")     # índice de rugosidad del terreno
asp  <- terrain(dem, "aspect",  unit="degrees")

# Proyectar y agregar a la grilla de análisis
capas <- c(dem, pend, rug, asp)
names(capas) <- c("elevacion","pendiente","rugosidad","orientacion")
capas_p <- project(capas, g$tierra, method="bilinear")
capas_p <- mask(capas_p, g$tierra)

writeRaster(capas_p, file.path(RUTA,"datos/procesados/covar_topo_2km.tif"), overwrite=TRUE)
print(summary(values(capas_p), digits=3))
cat("\nCovariables topográficas escritas a 2 km.\n")
