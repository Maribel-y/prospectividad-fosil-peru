# ============================================================
# Fase 1 — Unidad de análisis: grilla nacional de 2 km
# Asigna la variable respuesta y proyecta (snap) las costeras.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})

RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
RES  <- 2000   # 2 km, justificado por la precisión posicional del inventario

# ---- Territorio y grilla ------------------------------------
peru  <- geodata::gadm("PER", level = 0, path = file.path(RUTA, "datos/crudos"))
peru_p <- project(peru, PROJ)
plantilla <- rast(ext(peru_p), resolution = RES, crs = PROJ)
mascara <- rasterize(peru_p, plantilla, field = 1)
names(mascara) <- "tierra"

n_tot <- ncell(mascara); n_ter <- sum(!is.na(values(mascara)))
cat("Grilla:", nrow(mascara), "x", ncol(mascara), "celdas de", RES/1000, "km\n")
cat("Celdas terrestres:", n_ter, "(", round(n_ter*RES^2/1e6), "km2 )\n\n")

# ---- Presencias ---------------------------------------------
occ <- fread(file.path(RUTA, "datos/procesados/ocurrencias_final.csv"))
pts <- project(vect(occ, geom = c("lon","lat"), crs = "EPSG:4326"), PROJ)

celda <- cellFromXY(mascara, crds(pts))
en_tierra <- !is.na(celda) & !is.na(values(mascara)[celda])
cat("Ocurrencias en celda terrestre directa:", sum(en_tierra), "\n")

# Snap: las costeras caen en celda marina o fuera; se reasignan a la
# celda terrestre más próxima dentro de 2 km. Si no la hay, se excluyen.
fuera <- which(!en_tierra)
if (length(fuera)) {
  tierra_xy <- xyFromCell(mascara, which(!is.na(values(mascara))))
  ids_tierra <- which(!is.na(values(mascara)))
  for (i in fuera) {
    p <- crds(pts)[i, ]
    d <- sqrt((tierra_xy[,1]-p[1])^2 + (tierra_xy[,2]-p[2])^2)
    j <- which.min(d)
    if (d[j] <= RES) celda[i] <- ids_tierra[j] else celda[i] <- NA
  }
}
recuperadas <- sum(!is.na(celda)) - sum(en_tierra)
cat("Recuperadas por snap (<= 2 km):", recuperadas, "\n")
cat("Excluidas definitivamente:", sum(is.na(celda)), "\n\n")

occ[, celda := celda]
occ_val <- occ[!is.na(celda)]

# ---- Variable respuesta --------------------------------------
resp <- occ_val[, .(ocurrencias = .N,
                    eras = paste(sort(unique(era[era != "" & !is.na(era)])), collapse = "|")),
                by = celda]
cat("===== UNIDAD DE ANÁLISIS =====\n")
res <- data.table(
  concepto = c("Celdas terrestres (fondo potencial)", "Celdas con presencia",
               "Prevalencia (%)", "Celdas con presencia mesozoica"),
  valor = c(n_ter, nrow(resp), round(100*nrow(resp)/n_ter, 3),
            nrow(occ_val[era == "Mesozoico", .N, by = celda])))
print(res)

pres <- rast(mascara); values(pres) <- NA
pres[resp$celda] <- 1
names(pres) <- "presencia"

writeRaster(c(mascara, pres), file.path(RUTA, "datos/procesados/grilla_2km.tif"),
            overwrite = TRUE)
fwrite(occ_val, file.path(RUTA, "datos/procesados/ocurrencias_con_celda.csv"))
fwrite(resp,    file.path(RUTA, "datos/procesados/celdas_presencia.csv"))
fwrite(res,     file.path(RUTA, "salidas/tablas/tabla_unidad_analisis.csv"))

# ---- Figura 1: inventario ------------------------------------
png(file.path(RUTA, "salidas/figuras/fig01_inventario.png"), width = 1600, height = 2000, res = 200)
plot(peru_p, col = "grey94", border = "grey60", axes = FALSE,
     main = "Localidades fósiles del Perú\n(inventario PBDB + GBIF, corte 18-ago-2026)")
xy <- xyFromCell(mascara, resp$celda)
mes <- occ_val[era == "Mesozoico", unique(celda)]
points(xy, pch = 16, cex = 0.35, col = "#4A6FA5")
points(xyFromCell(mascara, mes), pch = 16, cex = 0.5, col = "#C1440E")
legend("bottomleft", bty = "n", pch = 16, col = c("#4A6FA5", "#C1440E"),
       legend = c(paste0("Todas las localidades (n = ", nrow(resp), ")"),
                  paste0("Mesozoico (n = ", length(mes), ")")), cex = 0.8)
dev.off()
cat("\nFigura 1 escrita.\n")
