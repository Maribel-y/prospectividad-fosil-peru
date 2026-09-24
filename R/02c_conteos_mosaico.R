# ============================================================
# Reduce un mosaico de NDVI p90 de Sentinel-2 a ~20 m (guion 02b) a
# conteos por celda de la grilla nacional de 2 km:
#   validos  píxeles con dato que no son agua, nieve, urbano ni manglar
#   d015, d020, d025  píxeles válidos con NDVI p90 < 0,15 / 0,20 / 0,25
# Los conteos (y no las fracciones) permiten sumar los mosaicos vecinos
# en las celdas que caen sobre el borde entre dos mosaicos.
# Uso: Rscript R/02c_conteos_mosaico.R entrada.tif salida.tif
# ============================================================
suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0)
arg  <- commandArgs(trailingOnly = TRUE)
RUTA <- path.expand("~/Desktop/dinos")
g    <- rast(file.path(RUTA, "datos/procesados/grilla_2km.tif"))

ndvi <- rast(arg[1]); NAflag(ndvi) <- 255
wc <- vrt(list.files(file.path(RUTA, "datos/crudos/worldcover"), "_74m\\.tif$", full.names = TRUE),
          file.path(tempdir(), "wc.vrt"), overwrite = TRUE)
wc <- crop(wc, ndvi)
wc <- resample(disagg(wc, 4), ndvi, method = "near")   # malla de 74 m -> 20 m (factor 4 exacto)

dn <- function(u) round((u + 1) / 0.008)              # NDVI = DN * 0,008 - 1
valido <- !is.na(ndvi) & !(wc %in% c(0, 50, 70, 80, 95))
x <- c(valido,
       valido & ndvi < dn(0.15),
       valido & ndvi < dn(0.20),
       valido & ndvi < dn(0.25))
names(x) <- c("validos", "d015", "d020", "d025")
x <- aggregate(x, 10, fun = "sum", na.rm = TRUE)    # ~185 m antes de proyectar
y <- project(x, g$tierra, method = "sum")
y <- trim(y)
writeRaster(y, arg[2], overwrite = TRUE, datatype = "FLT4S")
