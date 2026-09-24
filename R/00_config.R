# ============================================================
# Configuración común del proyecto
#
# El límite nacional se define aquí y no en cada guion, de modo que
# cambiar de fuente cartográfica solo requiera editar este archivo.
# ============================================================

RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
RES  <- 2000        # lado de celda en metros

# Fuente del límite nacional: "IGN" (oficial peruano) o "GADM" (alternativa global).
# Para usar el del IGN, coloque el archivo en datos/crudos/limite_ign/ y
# ajuste FUENTE_LIMITE. El resto de los guiones no necesita cambios.
FUENTE_LIMITE <- "GADM"
RUTA_IGN      <- file.path(RUTA, "datos/crudos/limite_ign/limite_nacional_ign.shp")

limite_nacional <- function(proyectar = TRUE) {
  if (FUENTE_LIMITE == "IGN") {
    if (!file.exists(RUTA_IGN))
      stop("No se encuentra el límite del IGN en: ", RUTA_IGN,
           "\nDescárguelo del portal IDEP (https://portalgeo.idep.gob.pe) ",
           "o cambie FUENTE_LIMITE a \"GADM\".")
    v <- terra::vect(RUTA_IGN)
    v <- terra::aggregate(v)          # disolver divisiones internas
  } else {
    v <- geodata::gadm("PER", level = 0, path = file.path(RUTA, "datos/crudos"))
  }
  if (proyectar) v <- terra::project(v, PROJ)
  v
}

# Documenta en las salidas qué fuente se usó, para que las cifras sean trazables
cita_limite <- function() {
  if (FUENTE_LIMITE == "IGN")
    "Instituto Geográfico Nacional del Perú (IGN)"
  else
    "Global Administrative Areas (GADM), versión 4.1"
}
