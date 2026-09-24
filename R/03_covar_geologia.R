# ============================================================
# Fase 2 (parte B) — Covariables geológicas a 2 km
# Fuente: Mapa Geológico Integrado 1:50 000 (INGEMMET, ago-2025),
# capa Litología, descargada vía API REST de GEOCATMIN.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(sf)})
RUTA <- path.expand("~/Desktop/dinos")

g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))

# ---- Unir los lotes descargados ------------------------------
arch <- list.files(file.path(RUTA,"datos/crudos/geologia_50k"), "\\.geojson$", full.names=TRUE)
cat("Lotes:", length(arch), "\n")
lst <- lapply(arch, function(f) tryCatch(vect(f), error=function(e) NULL))
lst <- lst[!sapply(lst, is.null)]
geo <- do.call(rbind, lst)
cat("Polígonos unidos:", nrow(geo), "\n")
# Nota: NO se usa makeValid() — en terra 1.9.34 colapsa la capa entera
# (20 000 polígonos -> 18) dejando los atributos descolgados.
nv <- sum(!is.valid(geo))
cat("Geometrías inválidas:", nv, sprintf("(%.3f %%)\n", 100*nv/nrow(geo)))
gv  <- project(geo, crs(g))
writeVector(gv, file.path(RUTA,"datos/procesados/geologia_50k.gpkg"), overwrite=TRUE)

# ---- Limpieza de atributos -----------------------------------
d <- as.data.table(values(gv))
d[, edad_media := (as.numeric(E_MAX_MA) + as.numeric(E_MIN_MA))/2]
d[, era := fifelse(is.na(edad_media), NA_character_,
           fifelse(edad_media < 66, "Cenozoico",
           fifelse(edad_media < 252, "Mesozoico", "Paleozoico_o_mas")))]
# AMBIENTE_S mezcla texto y códigos del dominio DGR_AMBIETE_SED del servicio:
# 1 = Continental, 2 = De transición, 3 = Marino. Se decodifica antes de usarlo.
d[, amb_txt := trimws(as.character(AMBIENTE_S))]
d[, ambiente := fifelse(amb_txt %in% c("1","Continental"),   "Continental",
                fifelse(amb_txt %in% c("2","De transición","Transicional"), "Transicional",
                fifelse(amb_txt %in% c("3","Marino"),        "Marino", NA_character_)))]
d[, sedimentaria := grepl("[Ss]edimen", TIPO_UNIDAD) | !is.na(MEDIO_SEDI)]
gv$edad_media   <- d$edad_media
gv$era          <- d$era
gv$ambiente     <- d$ambiente
gv$sedimentaria <- d$sedimentaria

cat("\n-- unidades distintas:", uniqueN(d$UNIDAD), "\n")
cat("-- eras:\n"); print(d[, .N, by=era][order(-N)])
cat("-- ambientes:\n"); print(d[, .N, by=ambiente][order(-N)])

# ---- Rasterizar a la grilla ----------------------------------
r_edad <- rasterize(gv, g$tierra, field="edad_media", fun="mean")
r_era  <- rasterize(gv, g$tierra, field="era")
r_amb  <- rasterize(gv, g$tierra, field="ambiente")
r_sed  <- rasterize(gv, g$tierra, field="sedimentaria")
r_uni  <- rasterize(gv, g$tierra, field="UNIDAD")
r_lito <- rasterize(gv, g$tierra, field="LITOLOGIA")

capas <- c(r_edad, r_era, r_amb, r_sed, r_uni, r_lito)
names(capas) <- c("edad_media","era_geol","ambiente","sedimentaria","unidad","litologia")
capas <- mask(capas, g$tierra)
writeRaster(capas, file.path(RUTA,"datos/procesados/covar_geologia_2km.tif"), overwrite=TRUE)

cob <- sum(!is.na(values(capas$edad_media)))
cat("\nCeldas con geología asignada:", cob, "de", sum(!is.na(values(g$tierra))),
    sprintf("(%.1f %%)\n", 100*cob/sum(!is.na(values(g$tierra)))))
cat("Covariables geológicas escritas.\n")
