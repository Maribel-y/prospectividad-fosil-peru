# ============================================================
# Figuras compuestas para la versión reducida del manuscrito.
# Diez figuras sueltas se agrupan en cinco láminas de varios paneles.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
source("R/00_config.R")
peru <- limite_nacional()
FIG <- file.path(RUTA, "salidas/figuras")
g   <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
D   <- fread(file.path(RUTA,"datos/procesados/predicciones.csv"))
occ <- fread(file.path(RUTA,"datos/procesados/ocurrencias_con_celda.csv"))
e <- ext(peru); dx <- xmax(e)-xmin(e); dy <- ymax(e)-ymin(e)

panel <- function(letra, titulo) {
  mtext(bquote(bold(.(letra))), side=3, line=0.5, adj=0, cex=0.95)
  mtext(titulo, side=3, line=0.5, adj=0.5, cex=0.82, col="grey25")
}
lienzo <- function(margen_der=0.04) {
  plot(ext(xmin(e)-0.04*dx, xmax(e)+margen_der*dx, ymin(e)-0.05*dy, ymax(e)+0.05*dy),
       col=NA, border=NA, axes=FALSE, xlab="", ylab="")
}
mapa_decil <- function(r, paleta="YlOrRd", barra=TRUE, etiquetas=c("bajo","alto")) {
  q  <- unique(quantile(values(r), seq(0,1,0.1), na.rm=TRUE))
  rc <- classify(r, q, include.lowest=TRUE)
  col <- hcl.colors(length(q)-1, paleta, rev=TRUE)
  lienzo(if (barra) 0.20 else 0.04)
  plot(peru, col="grey96", border=NA, add=TRUE)
  plot(rc, add=TRUE, col=col, legend=FALSE)
  plot(peru, border="grey45", add=TRUE, lwd=0.6)
  if (barra) {
    yb <- seq(ymin(e)+0.28*dy, ymin(e)+0.60*dy, length.out=length(col)+1)
    rect(xmax(e)+0.04*dx, yb[-length(yb)], xmax(e)+0.10*dx, yb[-1], col=col, border=NA)
    text(xmax(e)+0.12*dx, yb[c(1,length(yb))], etiquetas, adj=0, cex=0.7)
  }
}

# ---- Lámina 1: inventario y diseño de validación ----------------
png(file.path(FIG,"comp1_datos.png"), width=2000, height=1350, res=210)
par(mfrow=c(1,2), mar=c(0.5,0.5,3,0.5), oma=c(0,0,0,0))
lienzo(); plot(peru, col="grey96", border="grey55", add=TRUE)
mes <- unique(occ[era=="Mesozoico", celda]); todas <- unique(occ$celda)
points(xyFromCell(g, todas), pch=16, cex=0.30, col=adjustcolor("#3A5F8F",0.75))
points(xyFromCell(g, mes),   pch=16, cex=0.45, col="#C1440E")
legend(xmin(e)-0.03*dx, ymin(e)+0.10*dy, bty="n", pch=16, cex=0.65, pt.cex=0.9,
       col=c(adjustcolor("#3A5F8F",0.75),"#C1440E"),
       legend=c(sprintf("Celdas con presencia (n = %d)", length(todas)),
                sprintf("Celdas mesozoicas (n = %d)", length(mes))))
panel("(a)", "Inventario nacional")
r_fold <- rast(g$tierra); values(r_fold) <- NA
Dm <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv")); r_fold[Dm$celda] <- Dm$fold
lienzo(); plot(r_fold, add=TRUE, col=hcl.colors(5,"Set 2"), legend=FALSE)
plot(peru, border="grey35", add=TRUE, lwd=0.6)
points(xyFromCell(g, todas), pch=16, cex=0.14, col="#00000070")
legend(xmin(e)-0.03*dx, ymin(e)+0.12*dy, bty="n", pch=15, cex=0.62, pt.cex=1.1,
       col=hcl.colors(5,"Set 2"), legend=paste("Pliegue",1:5), ncol=2)
panel("(b)", "Bloques espaciales de 100 km")
dev.off(); cat("comp1 lista\n")

# ---- Lámina 2: modelo base y desempeño --------------------------
png(file.path(FIG,"comp2_modelo.png"), width=2100, height=1500, res=210)
layout(matrix(c(1,2,1,3), nrow=2, byrow=TRUE), widths=c(1.05,1))
par(mar=c(0.5,0.5,3,0.5))
mapa_decil(rast(file.path(RUTA,"datos/procesados/favorabilidad_base.tif")))
points(D[pres==1,.(x,y)], pch=1, cex=0.22, col="#00000090", lwd=0.5)
panel("(a)", "Favorabilidad geológica base")

par(mar=c(4.2,4.2,3,1))
roc <- readRDS(file.path(RUTA,"datos/procesados/roc_folds.rds"))
plot(NA, xlim=0:1, ylim=0:1, xlab="Tasa de falsos positivos",
     ylab="Tasa de verdaderos positivos", cex.axis=0.8, cex.lab=0.85)
abline(0,1,lty=3,col="grey60")
for (k in 1:5) lines(roc[[k]]$fpr, roc[[k]]$tpr, col=hcl.colors(5,"Dark 3")[k], lwd=1.7)
legend("bottomright", bty="n", lwd=1.7, cex=0.65, col=hcl.colors(5,"Dark 3"),
       legend=paste("Pliegue",1:5))
panel("(b)", "Curvas ROC por pliegue")

imp <- fread(file.path(RUTA,"salidas/tablas/importancia.csv"))
par(mar=c(4.2,11,3,1))
barplot(rev(imp$Gain[1:8]), horiz=TRUE, names.arg=rev(imp$nombre[1:8]), las=1,
        cex.names=0.68, cex.axis=0.75, col="#3A5F8F", border=NA, xlab="Ganancia relativa")
panel("(c)", "Importancia de covariables")
dev.off(); cat("comp2 lista\n")

# ---- Lámina 3: el efecto del sesgo ------------------------------
r <- rast(file.path(RUTA,"datos/procesados/prospectividad_nacional.tif"))
png(file.path(FIG,"comp3_sesgo.png"), width=2000, height=1400, res=210)
par(mfrow=c(1,2), mar=c(0.5,0.5,3,0.5), oma=c(2.6,0,0,0))
mapa_decil(r[["prospectividad"]], barra=FALSE); panel("(a)", "Proyección observacional")
mapa_decil(r[["prospectividad_sin_sesgo"]], barra=TRUE); panel("(b)", "Sesgo de acceso neutralizado")
mtext(paste("Correlación con la distancia a vías: -0,61 en (a) frente a -0,22 en (b).",
            "El 41,4 % del decil más prospectivo cambia de lugar, y el AUC validado del mapa corregido es 0,061 menor."),
      side=1, outer=TRUE, line=0.8, cex=0.72, col="grey30")
dev.off(); cat("comp3 lista\n")

# ---- Lámina 4: prospectividad e incertidumbre -------------------
png(file.path(FIG,"comp4_productos.png"), width=2000, height=1400, res=210)
par(mfrow=c(1,2), mar=c(0.5,0.5,3,0.5))
mapa_decil(r[["prospectividad_sin_sesgo"]])
points(D[pres==1,.(x,y)], pch=1, cex=0.20, col="#00000080", lwd=0.5)
panel("(a)", "Prospectividad nacional")
mapa_decil(r[["incertidumbre"]], paleta="Blues", etiquetas=c("baja","alta"))
panel("(b)", "Incertidumbre entre pliegues")
dev.off(); cat("comp4 lista\n")

# ---- Lámina 5: estratificación y priorización -------------------
# Panel (b): la priorización que se recomienda en el texto, es decir,
# la robusta a la constante de neutralización (guion 16) para el modelo
# nacional, junto a la priorización del estrato mesozoico (guion 10).
rm_  <- rast(file.path(RUTA,"datos/procesados/prospectividad_mesozoico.tif"))
rob  <- fread(file.path(RUTA,"salidas/tablas/tabla_priorizacion_robusta.csv"))
pm   <- fread(file.path(RUTA,"salidas/tablas/tabla_priorizacion_mesozoico.csv"))
comunes <- length(intersect(rob$celda, pm$celda))
cat("comp5: celdas robustas nacionales", nrow(rob), "| mesozoicas", nrow(pm),
    "| en común", comunes, "\n")
png(file.path(FIG,"comp5_priorizacion.png"), width=2000, height=1400, res=210)
par(mfrow=c(1,2), mar=c(0.5,0.5,3,0.5), oma=c(2.2,0,0,0))
mapa_decil(rm_[["prospectividad_sin_sesgo"]])
points(D[pres_mes==1,.(x,y)], pch=1, cex=0.28, col="#00000090", lwd=0.6)
panel("(a)", "Prospectividad del Mesozoico")
lienzo(); plot(peru, col="grey96", border="grey55", add=TRUE)
points(rob[,.(x,y)], pch=15, cex=0.26, col="#C1440E")
points(pm[,.(x,y)],  pch=15, cex=0.26, col="#1F5C8B")
fmt <- function(n) format(n, big.mark=" ")
legend(xmin(e)-0.03*dx, ymin(e)+0.12*dy, bty="n", pch=15, cex=0.62, pt.cex=1.1,
       col=c("#C1440E","#1F5C8B"),
       legend=c(sprintf("Nacional robusta (%s)", fmt(nrow(rob))),
                sprintf("Mesozoica (%s)", fmt(nrow(pm)))))
panel("(b)", "Zonas prioritarias")
mtext(sprintf("Priorización nacional robusta: %s km² (%s %% del país). Celdas compartidas con la mesozoica: %d.",
              fmt(nrow(rob)*4), format(round(100*nrow(rob)*4/1285216, 2), decimal.mark=","), comunes),
      side=1, outer=TRUE, line=0.6, cex=0.72, col="grey30")
dev.off(); cat("comp5 lista\n")
