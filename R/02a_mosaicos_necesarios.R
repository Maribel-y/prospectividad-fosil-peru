# Lists the 1-degree Sentinel-2 tiles that contain at least one land cell
# of the 2 km grid, so that 02b does not download tiles outside Peru.
suppressMessages(library(terra))
RUTA <- path.expand("~/Desktop/dinos")
g  <- rast(file.path(RUTA, "datos/procesados/grilla_2km.tif"))$tierra
xy <- crds(project(as.points(g), "EPSG:4326"))
t  <- ifelse(xy[,2] > 0, sprintf("N%02dW%03d", floor(xy[,2]), ceiling(-xy[,1])),
                         sprintf("S%02dW%03d", ceiling(-xy[,2]), ceiling(-xy[,1])))
writeLines(sort(unique(t)), file.path(RUTA, "datos/crudos/sentinel2_conteos/necesarios.txt"))
cat(length(unique(t)), "tiles needed\n")
