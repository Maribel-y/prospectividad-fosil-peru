# Recupera localidades costeras aplicando un buffer marino al límite nacional
# y recalcula el N efectivo definitivo.
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"

occ <- fread(file.path(RUTA,"datos/procesados/ocurrencias_unificadas.csv"))
fuera <- fread(file.path(RUTA,"salidas/tablas/tabla_descartadas_diagnostico.csv"))

# Criterio: recuperar lo que esté a <= 5 km del territorio (costa + frontera),
# descartar sondeos oceánicos profundos y errores de georreferencia.
rec <- fuera[dist_km <= 5, .(fuente, id=as.character(id), lon, lat,
                            taxon=NA_character_, clase=NA_character_,
                            max_ma, formacion)]
rec[, era := fifelse(is.na(max_ma), NA_character_,
             fifelse(max_ma < 66, "Cenozoico",
             fifelse(max_ma < 252, "Mesozoico", "Paleozoico")))]
cat("Recuperadas (<= 5 km):", nrow(rec), "ocurrencias\n")
cat("Descartadas definitivamente (> 5 km):", nrow(fuera[dist_km > 5]), "\n\n")

occ2 <- rbindlist(list(occ, rec), fill=TRUE)

# Coordenadas con precisión de grado (~111 km): lon y lat enteras. Son
# centroides de un grado, no localidades, y no pueden asignarse a una
# celda de 2 km. Incluye las 224 de PBDB con latlng_precision = "degrees"
# y los registros de GBIF equivalentes.
grado <- occ2[lon == round(lon) & lat == round(lat)]
cat("Excluidas por precisión de grado:", nrow(grado), "ocurrencias en",
    nrow(unique(grado[, .(lon, lat)])), "coordenadas\n\n")
occ  <- occ[!(lon == round(lon) & lat == round(lat))]
occ2 <- occ2[!(lon == round(lon) & lat == round(lat))]

dedup <- function(dt, d) {
  u <- unique(dt[, .(lon, lat)])
  if (nrow(u) < 2) return(nrow(u))
  xy <- crds(project(vect(u, geom=c("lon","lat"), crs="EPSG:4326"), PROJ))
  length(unique(cutree(hclust(dist(xy), method="complete"), h=d)))
}

cat("===== N EFECTIVO DEFINITIVO (buffer costero 5 km) =====\n")
comp <- data.table(
  version = c("Solo territorio estricto", "Con costeras recuperadas"),
  ocurrencias = c(nrow(occ), nrow(occ2)),
  N_1km = c(dedup(occ,1000), dedup(occ2,1000)),
  N_2km = c(dedup(occ,2000), dedup(occ2,2000))
)
print(comp)

cat("\n===== POR ERA, VERSIÓN DEFINITIVA =====\n")
pe <- occ2[!is.na(era), .(ocurrencias=.N), by=era]
pe[, N_2km := sapply(era, function(e) dedup(occ2[era==e], 2000))]
print(pe[order(-ocurrencias)])

fwrite(occ2, file.path(RUTA,"datos/procesados/ocurrencias_final.csv"))
fwrite(comp, file.path(RUTA,"salidas/tablas/tabla_censo_definitivo.csv"))
fwrite(pe,   file.path(RUTA,"salidas/tablas/tabla_censo_por_era_definitivo.csv"))
