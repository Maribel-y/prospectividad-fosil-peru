# ============================================================
# Fase 0 — Censo de viabilidad e inventario
# Proyecto: Prospectividad paleontológica del Perú
# Salida: N efectivo de localidades independientes
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(terra); library(rgbif)
})

RUTA  <- path.expand("~/Desktop/dinos")
CORTE <- Sys.Date()
PROJ  <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"

# ---- 1. PBDB ------------------------------------------------
f_pbdb <- file.path(RUTA, "datos/crudos", "pbdb_peru_2026-08-18.csv")
pbdb <- fread(f_pbdb)
pbdb_std <- data.table(
  fuente    = "PBDB",
  id        = as.character(pbdb$occurrence_no),
  lon       = as.numeric(pbdb$lng),
  lat       = as.numeric(pbdb$lat),
  taxon     = as.character(pbdb$accepted_name),
  clase     = as.character(pbdb[["class"]]),
  max_ma    = as.numeric(pbdb$max_ma),
  formacion = as.character(pbdb$formation)
)[!is.na(lon) & !is.na(lat)]
message("PBDB con coordenadas: ", nrow(pbdb_std))

# ---- 2. GBIF ------------------------------------------------
f_gbif <- file.path(RUTA, "datos/crudos", "gbif_peru.csv")
if (!file.exists(f_gbif)) {
  n <- occ_search(country = "PE", basisOfRecord = "FOSSIL_SPECIMEN",
                  hasCoordinate = TRUE, limit = 1)$meta$count
  message("GBIF anuncia ", n, " registros con coordenadas; descargando...")
  pags <- lapply(seq(0, n, by = 300), function(off) {
    d <- occ_search(country = "PE", basisOfRecord = "FOSSIL_SPECIMEN",
                    hasCoordinate = TRUE, limit = 300, start = off)$data
    if (is.null(d)) NULL else as.data.table(d)
  })
  gbif <- rbindlist(pags, fill = TRUE)
  fwrite(gbif, f_gbif)
} else gbif <- fread(f_gbif)
message("GBIF descargado: ", nrow(gbif))

col <- function(dt, nm) if (nm %in% names(dt)) as.character(dt[[nm]]) else NA_character_
gbif_std <- data.table(
  fuente    = "GBIF",
  id        = col(gbif, "key"),
  lon       = as.numeric(gbif$decimalLongitude),
  lat       = as.numeric(gbif$decimalLatitude),
  taxon     = col(gbif, "scientificName"),
  clase     = col(gbif, "class"),
  max_ma    = NA_real_,
  formacion = col(gbif, "formation")
)[!is.na(lon) & !is.na(lat)]

# ---- 3. Unificación ----------------------------------------
occ <- rbindlist(list(pbdb_std, gbif_std), fill = TRUE)
occ[, era := fifelse(is.na(max_ma), NA_character_,
              fifelse(max_ma < 66,  "Cenozoico",
              fifelse(max_ma < 252, "Mesozoico", "Paleozoico")))]
message("Ocurrencias unificadas: ", nrow(occ))

# ---- 4. Filtro al territorio nacional ------------------------
peru <- geodata::gadm("PER", level = 0, path = file.path(RUTA, "datos/crudos"))
pts  <- vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326")
dentro <- is.related(pts, peru, "intersects")
occ <- occ[dentro]
message("Dentro del Perú: ", nrow(occ), "  (descartadas: ", sum(!dentro), ")")

# ---- 5. Deduplicación espacial -------------------------------
# El N efectivo son localidades independientes, no ocurrencias:
# muchos fósiles comparten un mismo sitio de colecta.
dedup <- function(dt, dist_m) {
  u <- unique(dt[, .(lon, lat)])
  if (nrow(u) < 2) return(nrow(u))
  xy <- crds(project(vect(u, geom = c("lon","lat"), crs = "EPSG:4326"), PROJ))
  length(unique(cutree(hclust(dist(xy), method = "complete"), h = dist_m)))
}

message("Deduplicando...")
N1 <- dedup(occ, 1000); N2 <- dedup(occ, 2000)

censo <- data.table(
  etapa = c("Ocurrencias brutas (PBDB+GBIF)", "Con coordenadas", "Dentro del Perú",
            "Coordenadas únicas", "N efectivo (1 km)", "N efectivo (2 km)"),
  n = c(nrow(pbdb) + nrow(gbif), nrow(pbdb_std) + nrow(gbif_std), nrow(occ),
        nrow(unique(occ[, .(lon, lat)])), N1, N2)
)
cat("\n===== CENSO DE VIABILIDAD  (corte:", format(CORTE), ") =====\n"); print(censo)

# ---- 6. N efectivo por estrato -------------------------------
por_era <- occ[!is.na(era), .(ocurrencias = .N), by = era]
por_era[, N_2km := sapply(era, function(e) dedup(occ[era == e], 2000))]
setorder(por_era, -ocurrencias)
cat("\n===== POR ERA (solo PBDB tiene edad) =====\n"); print(por_era)

por_fuente <- occ[, .(ocurrencias = .N, N_2km = dedup(.SD, 2000)), by = fuente]
cat("\n===== POR FUENTE =====\n"); print(por_fuente)

# ---- 7. Guardar ----------------------------------------------
fwrite(occ,        file.path(RUTA, "datos/procesados", "ocurrencias_unificadas.csv"))
fwrite(censo,      file.path(RUTA, "salidas/tablas", "tabla_censo_viabilidad.csv"))
fwrite(por_era,    file.path(RUTA, "salidas/tablas", "tabla_censo_por_era.csv"))
fwrite(por_fuente, file.path(RUTA, "salidas/tablas", "tabla_censo_por_fuente.csv"))
cat("\nCenso cerrado con fecha de corte:", format(CORTE), "\n")
