# ============================================================
# Covariable de exposición de roca a partir de Sentinel-2.
#
# Fuente principal: percentil 90 anual de NDVI de Sentinel-2 (2021),
# compuesto de ESA WorldCover a ~20 m (guiones 02b y 02c_conteos_mosaico).
# Un píxel con NDVI p90 < 0,20 no alcanza cobertura vegetal apreciable
# en ningún momento del año: es superficie persistentemente desnuda
# (roca, regolito o sedimento suelto). A diferencia de la clase
# "desnudo / vegetación escasa" de WorldCover, el criterio distingue el
# afloramiento andino rodeado de pajonal estacional, que en época de
# lluvias sí reverdece, y a 20 m resuelve afloramientos de decenas de
# metros que a 74 m quedan diluidos en píxeles mixtos.
#
# Del denominador se excluyen agua (80), nieve y hielo (70), zonas
# urbanas (50) y manglar (95) según WorldCover 2021.
#
# Capas de salida (fracción 0-1 por celda de 2 km):
#   exposicion      NDVI p90 < 0,20   (covariable del modelo)
#   exposicion_015  NDVI p90 < 0,15   (sensibilidad al umbral)
#   exposicion_025  NDVI p90 < 0,25   (sensibilidad al umbral)
#   desnudo_wc      clase 60 de WorldCover (alternativa por clasificación)
# ============================================================
suppressPackageStartupMessages(library(terra))
RUTA <- path.expand("~/Desktop/dinos")
g    <- rast(file.path(RUTA, "datos/procesados/grilla_2km.tif"))

# ---- Sentinel-2: suma de conteos de todos los mosaicos ----------
f <- list.files(file.path(RUTA, "datos/crudos/sentinel2_conteos"), "^conteos_.*\\.tif$", full.names = TRUE)
cat("mosaicos con conteos:", length(f), "\n")
tot <- NULL
for (x in f) {
  y <- extend(rast(x), g$tierra)
  y <- resample(y, g$tierra, method = "near")        # misma malla: solo alinea extensiones
  y <- subst(y, NA, 0)
  tot <- if (is.null(tot)) y else tot + y
}
fr <- function(n) ifel(tot$validos > 0, tot[[n]] / tot$validos, 0)
expo <- c(fr("d020"), fr("d015"), fr("d025"))
names(expo) <- c("exposicion", "exposicion_015", "exposicion_025")
cat("celdas sin ningún píxel válido (lagos, glaciares):",
    global(tot$validos == 0 & !is.na(g$tierra), "sum", na.rm = TRUE)[[1]], "\n")

# ---- WorldCover clase 60, como alternativa ----------------------
tmp <- file.path(tempdir(), "wc_bin"); dir.create(tmp, showWarnings = FALSE)
bins <- vapply(list.files(file.path(RUTA, "datos/crudos/worldcover"), "_74m\\.tif$", full.names = TRUE),
  function(f) {
    out <- file.path(tmp, sub("_74m", "_bare", basename(f)))
    r <- rast(f); b <- ifel(r == 60, 1, ifel(r %in% c(0, 50, 70, 80, 95), NA, 0))
    writeRaster(b, out, overwrite = TRUE, datatype = "INT1U"); out
  }, "")
wc <- project(vrt(bins, file.path(tmp, "bare.vrt"), overwrite = TRUE), g$tierra, method = "average")
wc <- ifel(is.na(wc), 0, wc)
names(wc) <- "desnudo_wc"

expo <- mask(c(expo, wc), g$tierra)
writeRaster(expo, file.path(RUTA, "datos/procesados/covar_exposicion_2km.tif"), overwrite = TRUE)
for (n in names(expo))
  cat(sprintf("%-15s celdas %d | media %.3f\n", n,
              global(!is.na(expo[[n]]), "sum")[[1]], global(expo[[n]], "mean", na.rm = TRUE)[[1]]))
