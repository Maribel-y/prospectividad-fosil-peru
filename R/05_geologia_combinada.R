# ============================================================
# Fase 2 (parte B-2) — Capa geológica nacional combinada
#
# El Mapa Integrado 1:50 000 (INGEMMET, ago-2025) es el más rico en
# atributos —edad numérica, litología, ambiente sedimentario— pero solo
# cubre el 61 % del territorio. La Carta Geológica 1:100 000 cubre el
# país entero, con menos atributos. Se combinan: 50k como capa primaria
# y 100k como relleno, imputando edad y ambiente por nombre de unidad
# a partir del diccionario que aporta el propio 50k.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
geo50 <- vect(file.path(RUTA,"datos/procesados/geologia_50k.gpkg"))

# El gpkg se escribió antes de derivar los atributos: se recalculan aquí
# con el mismo criterio del script 03 (dominio DGR_AMBIETE_SED del servicio).
if (!"edad_media" %in% names(geo50)) {
  x <- as.data.table(values(geo50))
  x[, edad_media := (as.numeric(E_MAX_MA) + as.numeric(E_MIN_MA))/2]
  x[, amb_txt := trimws(as.character(AMBIENTE_S))]
  x[, ambiente := fifelse(amb_txt %in% c("1","Continental"), "Continental",
                  fifelse(amb_txt %in% c("2","De transición","Transicional"), "Transicional",
                  fifelse(amb_txt %in% c("3","Marino"), "Marino", NA_character_)))]
  x[, era := fifelse(is.na(edad_media), NA_character_,
             fifelse(edad_media < 66, "Cenozoico",
             fifelse(edad_media < 252, "Mesozoico", "Paleozoico_o_mas")))]
  geo50$edad_media <- x$edad_media
  geo50$ambiente   <- x$ambiente
  geo50$era        <- x$era
  writeVector(geo50, file.path(RUTA,"datos/procesados/geologia_50k.gpkg"), overwrite = TRUE)
}

norm <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[áàä]","a",x); x <- gsub("[éèë]","e",x); x <- gsub("[íìï]","i",x)
  x <- gsub("[óòö]","o",x); x <- gsub("[úùü]","u",x)
  gsub("\\s+"," ", gsub("[^a-z0-9 ]"," ", x))
}

# ---- Diccionario de unidades a partir del 50k ----------------
d50 <- as.data.table(values(geo50))
d50[, u := norm(UNIDAD)]
dic <- d50[!is.na(edad_media), .(
  edad_dic = median(edad_media, na.rm = TRUE),
  amb_dic  = if (all(is.na(ambiente))) NA_character_
             else names(sort(table(ambiente[!is.na(ambiente)]), decreasing = TRUE))[1],
  n = .N), by = u]
cat("Diccionario de unidades:", nrow(dic), "entradas\n")

# ---- Carta 1:100 000 ------------------------------------------
arch <- list.files(file.path(RUTA,"datos/crudos/geologia_100k"), "\\.geojson$", full.names = TRUE)
cat("Lotes 100k:", length(arch), "\n")
geo100 <- do.call(rbind, lapply(arch, function(f) tryCatch(vect(f), error = function(e) NULL)))
geo100 <- project(geo100, crs(g))
cat("Polígonos 100k:", nrow(geo100), "\n")

d100 <- as.data.table(values(geo100))
d100[, u := norm(UNIDAD)]
d100 <- merge(d100, dic[, .(u, edad_dic, amb_dic)], by = "u", all.x = TRUE, sort = FALSE)
cobertura_dic <- 100 * sum(!is.na(d100$edad_dic)) / nrow(d100)
cat(sprintf("Polígonos 100k con edad imputada por diccionario: %.1f %%\n", cobertura_dic))

geo100$edad_media <- d100$edad_dic
geo100$ambiente   <- d100$amb_dic
geo100$era <- fifelse(is.na(d100$edad_dic), NA_character_,
              fifelse(d100$edad_dic < 66,  "Cenozoico",
              fifelse(d100$edad_dic < 252, "Mesozoico", "Paleozoico_o_mas")))
geo100$fuente_geol <- "100k"
geo50$fuente_geol  <- "50k"

# ---- Rasterizado por prioridad --------------------------------
# terra no fusiona los diccionarios de categorías al hacer cover(), así que
# las capas categóricas se codifican a enteros con niveles COMPARTIDOS entre
# ambas cartas; las etiquetas se reponen al final.
categ <- c("era","ambiente","UNIDAD")
niveles <- lapply(categ, function(f) {
  v <- unique(c(as.character(values(geo50)[[f]]), as.character(values(geo100)[[f]])))
  sort(v[!is.na(v)])
})
names(niveles) <- categ
for (f in categ) {
  geo50[[paste0(f,"_cod")]]  <- match(as.character(values(geo50)[[f]]),  niveles[[f]])
  geo100[[paste0(f,"_cod")]] <- match(as.character(values(geo100)[[f]]), niveles[[f]])
}

campos <- c("edad_media", paste0(categ, "_cod"))
r50  <- rast(lapply(campos, function(f) rasterize(geo50,  g$tierra, field = f)))
r100 <- rast(lapply(campos, function(f) rasterize(geo100, g$tierra, field = f)))
names(r50) <- names(r100) <- campos

comb <- cover(r50, r100)          # 50k manda; 100k rellena
for (f in categ) {
  i <- which(campos == paste0(f,"_cod"))
  levels(comb[[i]]) <- data.frame(ID = seq_along(niveles[[f]]), etiqueta = niveles[[f]])
}
geo50$fuente_cod <- 1L; geo100$fuente_cod <- 2L
proc <- rasterize(geo50, g$tierra, field = "fuente_cod")
proc <- cover(proc, rasterize(geo100, g$tierra, field = "fuente_cod"))
levels(proc) <- data.frame(ID = 1:2, etiqueta = c("50k","100k"))
names(proc) <- "fuente_geol"

capas <- mask(c(comb, proc), g$tierra)
names(capas) <- c("edad_media","era_geol","ambiente","unidad","fuente_geol")
writeRaster(capas, file.path(RUTA,"datos/procesados/covar_geologia_2km.tif"), overwrite = TRUE)

ter <- sum(!is.na(values(g$tierra)))
tab <- data.table(
  capa = c("Solo 1:50 000", "Combinada 50k + 100k", "Con edad asignada", "Con ambiente asignado"),
  celdas = c(sum(!is.na(values(r50$edad_media))), sum(!is.na(values(capas$unidad))),
             sum(!is.na(values(capas$edad_media))), sum(!is.na(values(capas$ambiente)))))
tab[, pct := round(100*celdas/ter, 1)]
cat("\n===== COBERTURA GEOLÓGICA (", ter, "celdas terrestres ) =====\n"); print(tab)
fwrite(tab, file.path(RUTA,"salidas/tablas/tabla_cobertura_geologica.csv"))
cat("\nCapa geológica combinada escrita.\n")
