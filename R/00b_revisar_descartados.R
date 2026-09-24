# Diagnóstico de las ocurrencias que caen fuera del límite terrestre del Perú
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")

pbdb <- fread(file.path(RUTA,"datos/crudos/pbdb_peru_2026-08-18.csv"))
gbif <- fread(file.path(RUTA,"datos/crudos/gbif_peru.csv"))
occ <- rbindlist(list(
  data.table(fuente="PBDB", id=as.character(pbdb$occurrence_no),
             lon=as.numeric(pbdb$lng), lat=as.numeric(pbdb$lat),
             sitio=as.character(pbdb$county), precision=as.character(pbdb$latlng_precision),
             formacion=as.character(pbdb$formation), max_ma=as.numeric(pbdb$max_ma)),
  data.table(fuente="GBIF", id=as.character(gbif$key),
             lon=as.numeric(gbif$decimalLongitude), lat=as.numeric(gbif$decimalLatitude),
             sitio=as.character(gbif$locality), precision=NA_character_,
             formacion=NA_character_, max_ma=NA_real_)
), fill=TRUE)[!is.na(lon) & !is.na(lat)]

peru <- geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos"))
pts  <- vect(occ, geom=c("lon","lat"), crs="EPSG:4326")
fuera <- occ[!is.related(pts, peru, "intersects")]
cat("Descartadas:", nrow(fuera), "\n\n")

pf <- vect(fuera, geom=c("lon","lat"), crs="EPSG:4326")
fuera[, dist_km := round(as.numeric(distance(pf, peru))/1000, 2)]

# ¿En qué país caen?
mundo <- geodata::world(path=file.path(RUTA,"datos/crudos"))
ext <- extract(mundo, pf)
fuera[, pais := ext[[2]]]

cat("===== POR DISTANCIA AL TERRITORIO =====\n")
fuera[, rango := cut(dist_km, c(-1,1,5,20,100,500,Inf),
       labels=c("< 1 km (borde)","1-5 km","5-20 km","20-100 km","100-500 km","> 500 km"))]
print(fuera[, .(n=.N), by=rango][order(rango)])

cat("\n===== POR PAÍS (NA = en el mar) =====\n")
print(fuera[, .(n=.N, dist_media_km=round(mean(dist_km),1)), by=pais][order(-n)])

cat("\n===== POR FUENTE =====\n")
print(fuera[, .(n=.N, mediana_km=round(median(dist_km),1)), by=fuente])

cat("\n===== LOS 15 MÁS LEJANOS =====\n")
print(head(fuera[order(-dist_km), .(fuente, lon, lat, dist_km, pais, sitio=substr(sitio,1,28))], 15))

cat("\n===== RECUPERABLES: en el mar y a menos de 20 km de la costa =====\n")
rec <- fuera[is.na(pais) & dist_km <= 20]
cat("n =", nrow(rec), " | localidades únicas:", nrow(unique(rec[,.(lon,lat)])), "\n")
print(head(rec[order(-dist_km), .(fuente, lon, lat, dist_km, formacion, sitio=substr(sitio,1,25))], 12))

fwrite(fuera, file.path(RUTA,"salidas/tablas/tabla_descartadas_diagnostico.csv"))
