# ============================================================
# Figuras 3 a 6 del artículo
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
peru <- project(geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos")), PROJ)
D <- fread(file.path(RUTA,"datos/procesados/predicciones.csv"))
e <- ext(peru); dx <- xmax(e)-xmin(e); dy <- ymax(e)-ymin(e)

mapa <- function(r, titulo, subtitulo, paleta, etiquetas, archivo,
                 puntos = NULL, leyenda_puntos = NULL) {
  q  <- unique(quantile(values(r), seq(0,1,0.1), na.rm=TRUE))
  rc <- classify(r, q, include.lowest=TRUE)
  col <- hcl.colors(length(q)-1, paleta, rev = TRUE)   # siempre: oscuro = valor alto
  png(file.path(RUTA,"salidas/figuras",archivo), width=1700, height=2150, res=210)
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(ext(xmin(e)-0.04*dx, xmax(e)+0.18*dx, ymin(e)-0.12*dy, ymax(e)+0.10*dy),
       col=NA, border=NA, axes=FALSE, xlab="", ylab="")
  plot(peru, col="grey96", border=NA, add=TRUE)
  plot(rc, add=TRUE, col=col, legend=FALSE)
  plot(peru, border="grey45", add=TRUE)
  if (!is.null(puntos)) points(puntos, pch=1, cex=0.3, col="#00000095", lwd=0.6)
  text(xmin(e)+0.5*dx, ymax(e)+0.075*dy, titulo, cex=1.35, font=2)
  text(xmin(e)+0.5*dx, ymax(e)+0.035*dy, subtitulo, cex=0.8, col="grey30")
  yb <- seq(ymin(e)+0.30*dy, ymin(e)+0.62*dy, length.out=length(col)+1)
  rect(xmax(e)+0.03*dx, yb[-length(yb)], xmax(e)+0.08*dx, yb[-1], col=col, border=NA)
  text(xmax(e)+0.095*dx, yb[c(1, round(length(yb)/2), length(yb))], etiquetas, adj=0, cex=0.75)
  text(xmax(e)+0.03*dx, ymin(e)+0.66*dy, "Decil", adj=0, cex=0.78, font=2)
  if (!is.null(leyenda_puntos))
    legend(xmin(e)-0.02*dx, ymin(e)+0.03*dy, bty="n", pch=1, pt.cex=0.9, cex=0.8,
           col="black", legend=leyenda_puntos)
  dev.off()
}

for (cfg in list(
  list(suf="nacional",  capa="prospectividad_sin_sesgo", pal="YlOrRd",
       tit="Prospectividad paleontológica del Perú",
       sub="Gradient Boosting · AUC 0,859 con bloqueo espacial · sesgo de acceso neutralizado",
       etq=c("bajo","medio","alto"), arch="fig03_prospectividad_nacional.png",
       pts=as.matrix(D[pres==1, .(x,y)]), lp="Localidad fósil conocida"),
  list(suf="nacional",  capa="incertidumbre", pal="Blues",
       tit="Incertidumbre de la predicción",
       sub="Desviación estándar entre los cinco modelos de la validación espacial",
       etq=c("baja","media","alta"), arch="fig04_incertidumbre_nacional.png", pts=NULL, lp=NULL),
  list(suf="mesozoico", capa="prospectividad_sin_sesgo", pal="YlOrRd",
       tit="Prospectividad del Mesozoico",
       sub="Modelo estratificado con target-group background · AUC 0,834",
       etq=c("bajo","medio","alto"), arch="fig05_prospectividad_mesozoico.png",
       pts=as.matrix(D[pres_mes==1, .(x,y)]), lp="Localidad mesozoica conocida"))) {
  r <- rast(file.path(RUTA, paste0("datos/procesados/prospectividad_", cfg$suf, ".tif")))[[cfg$capa]]
  mapa(r, cfg$tit, cfg$sub, cfg$pal, cfg$etq, cfg$arch, cfg$pts, cfg$lp)
  cat("escrita:", cfg$arch, "\n")
}

# ---- Figura 6: priorización de campaña ------------------------
r <- rast(file.path(RUTA,"datos/procesados/prospectividad_nacional.tif"))
prio <- fread(file.path(RUTA,"salidas/tablas/tabla_priorizacion_nacional.csv"))
pm   <- fread(file.path(RUTA,"salidas/tablas/tabla_priorizacion_mesozoico.csv"))
png(file.path(RUTA,"salidas/figuras/fig06_priorizacion.png"), width=1700, height=2150, res=210)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ext(xmin(e)-0.04*dx, xmax(e)+0.10*dx, ymin(e)-0.14*dy, ymax(e)+0.10*dy),
     col=NA, border=NA, axes=FALSE, xlab="", ylab="")
plot(peru, col="grey96", border="grey55", add=TRUE)
points(prio[, .(x,y)], pch=15, cex=0.30, col="#C1440E")
points(pm[,  .(x,y)], pch=15, cex=0.30, col="#1F5C8B")
points(D[pres==1, .(x,y)], pch=1, cex=0.22, col="#00000060", lwd=0.5)
text(xmin(e)+0.5*dx, ymax(e)+0.075*dy, "Zonas prioritarias de prospección", cex=1.35, font=2)
text(xmin(e)+0.5*dx, ymax(e)+0.035*dy,
     "1 % de mayor potencial y cuartil de menor incertidumbre", cex=0.8, col="grey30")
legend(xmin(e)-0.02*dx, ymin(e)+0.05*dy, bty="n", pch=c(15,15,1), pt.cex=c(1.1,1.1,0.9), cex=0.8,
       col=c("#C1440E","#1F5C8B","black"),
       legend=c(sprintf("Prioridad nacional (%d celdas)", nrow(prio)),
                sprintf("Prioridad mesozoica (%d celdas)", nrow(pm)),
                "Localidad ya conocida"))
mtext(sprintf("Cada celda son 4 km². Superficie priorizada: %s km² (%.2f %% del país)",
              format((nrow(prio)+nrow(pm))*4, big.mark=" "),
              100*(nrow(prio)+nrow(pm))*4/1285216), side=1, line=-2, cex=0.72, col="grey35")
dev.off()
cat("escrita: fig06_priorizacion.png\n")
