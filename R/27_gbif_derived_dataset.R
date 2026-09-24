# ============================================================
# Table for registering a GBIF derived dataset (DOI for the GBIF data).
# GBIF asks for the number of records used from each source dataset
# (datasetKey, count). The GBIF records kept in the final inventory are
# identified by coordinates: every coordinate retained keeps all the raw
# records located there, so the match is exact.
# Upload the output at https://www.gbif.org/derived-dataset/register
# ============================================================
suppressPackageStartupMessages(library(data.table))
RUTA <- path.expand("~/Desktop/dinos")
raw <- fread(file.path(RUTA, "datos/crudos/gbif_peru.csv"),
             select = c("datasetKey", "decimalLongitude", "decimalLatitude"))
fin <- unique(fread(file.path(RUTA, "datos/procesados/ocurrencias_final.csv"))[fuente == "GBIF", .(lon, lat)])
usados <- merge(raw, fin, by.x = c("decimalLongitude", "decimalLatitude"), by.y = c("lon", "lat"))
tab <- usados[, .(count = .N), by = datasetKey][order(-count)]
cat("GBIF records used:", sum(tab$count), "of", nrow(raw), "| source datasets:", nrow(tab), "\n")
fwrite(tab, file.path(RUTA, "salidas/gbif_derived_dataset.csv"))
