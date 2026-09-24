# Figura 1 del artículo — inventario nacional de localidades fósiles
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"

g    <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
occ  <- fread(file.path(RUTA,"datos/procesados/ocurrencias_con_celda.csv"))
peru <- project(geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos")), PROJ)

celdas <- unique(occ$celda); mes <- unique(occ[era=="Mesozoico", celda])
xy <- xyFromCell(g, celdas); xym <- xyFromCell(g, mes)
e  <- ext(peru); dx <- xmax(e)-xmin(e); dy <- ymax(e)-ymin(e)
# margen extra abajo para la leyenda y arriba para el título
e2 <- ext(xmin(e)-0.04*dx, xmax(e)+0.04*dx, ymin(e)-0.14*dy, ymax(e)+0.10*dy)

png(file.path(RUTA,"salidas/figuras/fig01_inventario.png"), width=1700, height=2150, res=210)
par(mar=c(0.5,0.5,0.5,0.5))
plot(e2, col=NA, border=NA, axes=FALSE, xlab="", ylab="")
plot(peru, col="grey95", border="grey55", add=TRUE)
points(xy,  pch=16, cex=0.40, col=adjustcolor("#3A5F8F",0.75))
points(xym, pch=16, cex=0.55, col="#C1440E")

text(xmin(e)+0.5*dx, ymax(e)+0.075*dy, "Localidades fósiles del Perú", cex=1.4, font=2)
text(xmin(e)+0.5*dx, ymax(e)+0.035*dy,
     "Inventario PBDB + GBIF · corte 18-ago-2026 · grilla de 2 km", cex=0.85, col="grey30")
legend(xmin(e)-0.02*dx, ymin(e)+0.02*dy, bty="n", pch=16, pt.cex=1.3, cex=0.88, y.intersp=1.3,
       col=c(adjustcolor("#3A5F8F",0.75), "#C1440E"),
       legend=c(sprintf("Todas las localidades (n = %d celdas)", length(celdas)),
                sprintf("Con registro mesozoico (n = %d)", length(mes))))
text(xmin(e)-0.02*dx, ymin(e)-0.115*dy, adj=0, cex=0.75, col="grey35",
     sprintf("Prevalencia: %.2f %% de %s celdas terrestres de 2 km",
             100*length(celdas)/sum(!is.na(values(g$tierra))),
             format(sum(!is.na(values(g$tierra))), big.mark=" ")))
dev.off()
cat("figura regenerada\n")
