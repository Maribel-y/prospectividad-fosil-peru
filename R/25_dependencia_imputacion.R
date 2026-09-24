# ============================================================
# Do the priority areas depend disproportionately on imputed geology?
# Cells outside the 1:50 000 map take unit from the 1:100 000 map and
# age and environment from a dictionary of unit names (script 05).
# For each product we compare the share of 1:100 000 cells, and of cells
# with age or environment missing, against the whole territory; a ratio
# above 1 means the product leans on the poorer geological source.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
geo <- rast(file.path(RUTA, "datos/procesados/covar_geologia_2km.grd"))
P <- fread(file.path(RUTA, "datos/procesados/predicciones.csv"))
P[, fuente := values(geo$fuente_geol)[celda]]
P[, sin_edad := is.na(values(geo$edad_media)[celda])]
P[, sin_amb := is.na(values(geo$ambiente)[celda])]

top <- function(v, q) v >= quantile(v, q)
rob <- fread(file.path(RUTA, "salidas/tablas/tabla_priorizacion_robusta.csv"))$celda
pn  <- fread(file.path(RUTA, "salidas/tablas/tabla_priorizacion_nacional.csv"))$celda
pm  <- fread(file.path(RUTA, "salidas/tablas/tabla_priorizacion_mesozoico.csv"))$celda
conj <- list("Territory" = rep(TRUE, nrow(P)),
             "Presence cells" = P$pres == 1,
             "Top decile, neutralised national" = top(P$prospcorr_nacional, 0.9),
             "Top 1 %, neutralised national" = top(P$prospcorr_nacional, 0.99),
             "Single-projection priority set (national)" = P$celda %in% pn,
             "Robust priority set (national)" = P$celda %in% rob,
             "Mesozoic priority set" = P$celda %in% pm)
ref <- mean(P$fuente == 2, na.rm = TRUE)
tab <- rbindlist(lapply(names(conj), function(n) { s <- P[conj[[n]]]
  data.table(conjunto = n, celdas = nrow(s),
             pct_100k = 100*mean(s$fuente == 2, na.rm = TRUE),
             razon_100k = mean(s$fuente == 2, na.rm = TRUE)/ref,
             pct_sin_edad = 100*mean(s$sin_edad), pct_sin_ambiente = 100*mean(s$sin_amb)) }))
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
fwrite(tab, file.path(RUTA, "salidas/tablas/tabla_dependencia_imputacion.csv"))
