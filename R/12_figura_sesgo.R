# Figura 7 — evidencia visual de la tesis: el mismo modelo, proyectado con
# y sin las covariables de accesibilidad, recomienda buscar en sitios distintos.
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
peru <- project(geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos")), PROJ)
r <- rast(file.path(RUTA,"datos/procesados/prospectividad_nacional.tif"))

png(file.path(RUTA,"salidas/figuras/fig07_efecto_sesgo.png"), width=2200, height=1500, res=200)
par(mfrow=c(1,2), mar=c(1,0.5,4,0.5), oma=c(4,0,3,0))
for (i in seq_along(c("prospectividad","prospectividad_sin_sesgo"))) {
  capa <- c("prospectividad","prospectividad_sin_sesgo")[i]
  tit  <- c("(a) Proyección observacional","(b) Sesgo de acceso neutralizado")[i]
  sub  <- c("las covariables de accesibilidad actúan en la predicción",
            "fijadas al valor mediano de las celdas con presencia")[i]
  q  <- unique(quantile(values(r[[capa]]), seq(0,1,0.1), na.rm=TRUE))
  rc <- classify(r[[capa]], q, include.lowest=TRUE)
  plot(peru, col="grey96", border=NA, axes=FALSE, mar=c(1,0.5,4,0.5))
  plot(rc, add=TRUE, col=hcl.colors(length(q)-1,"YlOrRd",rev=TRUE), legend=FALSE)
  plot(peru, border="grey45", add=TRUE)
  title(main=tit, line=1.6, cex.main=1.15)
  mtext(sub, side=3, line=0.4, cex=0.78, col="grey35")
}
mtext("El mismo modelo, dos proyecciones", outer=TRUE, line=0.8, cex=1.4, font=2)
mtext(paste("Correlación de Spearman con la distancia a vías: -0,61 en (a) frente a -0,22 en (b);",
            "con la distancia a poblados: -0,53 frente a -0,19.",
            "\nEl 41,4 % del decil más prospectivo cambia de lugar entre ambas proyecciones,",
            "y el AUC validado del mapa corregido es 0,061 menor."),
      side=1, outer=TRUE, line=1.1, cex=0.78, col="grey30")
dev.off()
cat("figura 7 escrita\n")
