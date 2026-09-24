# Figura 2 — favorabilidad geológica base (modelo de línea base)
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
r    <- rast(file.path(RUTA,"datos/procesados/favorabilidad_base.tif"))
D    <- fread(file.path(RUTA,"datos/procesados/tabla_analisis.csv"))
peru <- project(geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos")), PROJ)

# escala en percentiles: la probabilidad absoluta es muy baja (prevalencia 0,22 %)
q <- quantile(values(r), seq(0,1,0.1), na.rm=TRUE)
rc <- classify(r, unique(q), include.lowest=TRUE)

e <- ext(peru); dx <- xmax(e)-xmin(e); dy <- ymax(e)-ymin(e)
png(file.path(RUTA,"salidas/figuras/fig02_favorabilidad_base.png"), width=1700, height=2150, res=210)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ext(xmin(e)-0.04*dx, xmax(e)+0.16*dx, ymin(e)-0.12*dy, ymax(e)+0.10*dy),
     col=NA, border=NA, axes=FALSE, xlab="", ylab="")
plot(peru, col="grey96", border=NA, add=TRUE)
plot(rc, add=TRUE, col=hcl.colors(10, "YlOrRd", rev=TRUE), legend=FALSE)
plot(peru, border="grey45", add=TRUE)
points(D[pres==1, .(x,y)], pch=1, cex=0.28, col="#00000090", lwd=0.6)

text(xmin(e)+0.5*dx, ymax(e)+0.075*dy, "Favorabilidad geológica base", cex=1.4, font=2)
text(xmin(e)+0.5*dx, ymax(e)+0.035*dy,
     "Pesos de Evidencia + regresión logística · AUC = 0,76 con bloqueo espacial", cex=0.82, col="grey30")
leg <- hcl.colors(10, "YlOrRd", rev=TRUE)
yb <- seq(ymin(e)+0.30*dy, ymin(e)+0.62*dy, length.out=11)
rect(xmax(e)+0.02*dx, yb[-11], xmax(e)+0.07*dx, yb[-1], col=leg, border=NA)
text(xmax(e)+0.085*dx, yb[c(1,6,11)], c("bajo","medio","alto"), adj=0, cex=0.75)
text(xmax(e)+0.02*dx, ymin(e)+0.66*dy, "Decil de\nfavorabilidad", adj=0, cex=0.78, font=2)
legend(xmin(e)-0.02*dx, ymin(e)+0.03*dy, bty="n", pch=1, pt.cex=0.9, cex=0.8,
       col="black", legend="Localidad fósil conocida")
dev.off()
cat("figura 2 escrita\n")
