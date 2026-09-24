# ============================================================
# Limpieza de la capa geológica combinada.
#
# Al escribir un GeoTIFF multicapa con capas categóricas, las celdas sin
# dato quedaron con el centinela entero de GDAL (-2147483648) en vez de
# NoData, y R las leía como una categoría más (inflaba la cobertura al
# 99,8 %). Aquí se convierten en ausencia real y el resultado se guarda
# en formato nativo de terra (.grd), que conserva categorías y NA.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
CENTINELA <- -2147483648

r <- rast(file.path(RUTA,"datos/procesados/covar_geologia_2km.tif"))
names(r) <- c("edad_media","era_geol","ambiente","unidad","fuente_geol")
g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))

for (n in names(r)) {
  v <- values(r[[n]])
  v[v == CENTINELA] <- NA
  values(r[[n]]) <- v
}
r <- mask(r, g$tierra)

sal <- file.path(RUTA,"datos/procesados/covar_geologia_2km.grd")
writeRaster(r, sal, overwrite = TRUE)

ter <- sum(!is.na(values(g$tierra)))
cuenta <- function(n) sum(!is.na(values(r[[n]])))
tab <- data.table(
  variable = c("Unidad geológica","Era geológica","Edad numérica (Ma)",
               "Ambiente sedimentario","  fuente 1:50 000","  fuente 1:100 000"),
  celdas = c(cuenta("unidad"), cuenta("era_geol"), cuenta("edad_media"), cuenta("ambiente"),
             sum(values(r$fuente_geol) == 1, na.rm = TRUE),
             sum(values(r$fuente_geol) == 2, na.rm = TRUE)))
tab[, pct := round(100*celdas/ter, 1)]
cat("===== COBERTURA GEOLÓGICA REAL (", format(ter, big.mark=" "), "celdas ) =====\n")
print(tab)
cat("\n-- ambiente --\n"); print(freq(r$ambiente))
cat("\n-- era --\n");      print(freq(r$era_geol))
fwrite(tab, file.path(RUTA,"salidas/tablas/tabla_cobertura_geologica.csv"))
