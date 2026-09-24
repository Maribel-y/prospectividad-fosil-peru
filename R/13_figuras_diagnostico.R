# ============================================================
# Figuras de diagnóstico: importancia de variables, curvas ROC
# por pliegue espacial y diseño de bloques de validación.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(xgboost); library(ranger)})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
K <- 5; N_FONDO <- 10000
COVAR <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion",
           "dist_vias_km","dist_poblados_km","dist_rios_km","woe_uni")
FORM <- as.formula(paste("~", paste(COVAR, collapse=" + "), "- 1"))
woe_map <- function(tr){ N<-nrow(tr); Dn<-sum(tr$pres)
  w <- tr[, .(n=.N,d=sum(pres)), by=uni]; w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]
  setNames(w$Wmas, w$uni) }
apl <- function(m,u){ v<-m[u]; fifelse(is.na(v),0,as.numeric(v)) }
entrenar <- function(tr){ ent <- rbind(tr[pres==1], tr[pres==0][sample(.N,min(.N,N_FONDO))]); y<-ent$pres
  list(m = xgb.train(params=list(objective="binary:logistic",max_depth=5,eta=0.05,subsample=0.8,
                                 colsample_bytree=0.8,scale_pos_weight=sum(y==0)/sum(y==1),nthread=8),
                     data=xgb.DMatrix(model.matrix(FORM,ent),label=y), nrounds=300, verbose=0),
       ent = ent, y = y) }

# ---- Predicciones por pliegue (para las curvas ROC) -----------
roc_folds <- list()
for (k in 1:K) {
  tr <- D[fold!=k][!(pres==1 & unidad_ignea)]; te <- D[fold==k]
  wm <- woe_map(tr); tr[, woe_uni := apl(wm,uni)]; te[, woe_uni := apl(wm,uni)]
  p <- predict(entrenar(tr)$m, xgb.DMatrix(model.matrix(FORM, te)))
  o <- order(-p); y <- te$pres[o]
  roc_folds[[k]] <- data.table(fpr = cumsum(y==0)/sum(y==0), tpr = cumsum(y==1)/sum(y==1))
}

# ---- Modelo completo para la importancia ----------------------
# Importance of the final national model of script 10 (the mapped one)
imp <- as.data.table(xgb.importance(model = readRDS(file.path(RUTA,"datos/procesados/modelo_final_nacional.rds"))$modelo))[order(-Gain)][1:12]
etiqueta <- c(woe_uni="Unidad geológica (WoE)", dist_vias_km="Distancia a vías",
  dist_poblados_km="Distancia a poblados", dist_rios_km="Distancia a ríos",
  elevacion="Elevación", pendiente="Pendiente", rugosidad="Rugosidad del terreno", exposicion="Exposición de roca",
  log_edad="Edad (log Ma)", edad_na="Edad no disponible", ambMarino="Ambiente marino",
  ambContinental="Ambiente continental", ambTransicional="Ambiente transicional",
  `ambSin dato`="Ambiente sin dato", eraCenozoico="Era cenozoica", eraMesozoico="Era mesozoica",
  `eraPaleozoico_o_mas`="Era paleozoica", `eraSin dato`="Era sin dato")
imp[, nombre := fifelse(is.na(etiqueta[Feature]), Feature, etiqueta[Feature])]
# Se guardan para poder recomponer las figuras sin reentrenar
saveRDS(roc_folds, file.path(RUTA,"datos/procesados/roc_folds.rds"))
fwrite(imp, file.path(RUTA,"salidas/tablas/importancia.csv"))

png(file.path(RUTA,"salidas/figuras/fig08_importancia.png"), width=1700, height=1250, res=210)
par(mar=c(4.2,13,3.2,1.5))
b <- barplot(rev(imp$Gain), horiz=TRUE, names.arg=rev(imp$nombre), las=1, cex.names=0.82,
             col="#3A5F8F", border=NA, xlab="Ganancia relativa")
title(main="Importancia de las covariables", line=1.8, cex.main=1.2)
mtext("Gradient boosting, modelo nacional ajustado sobre todo el país", side=3, line=0.4,
      cex=0.78, col="grey35")
dev.off(); cat("fig08 escrita\n")

# ---- Curvas ROC por pliegue -----------------------------------
png(file.path(RUTA,"salidas/figuras/fig09_roc.png"), width=1500, height=1500, res=210)
par(mar=c(4.4,4.4,3.6,1.2))
plot(NA, xlim=0:1, ylim=0:1, xlab="Tasa de falsos positivos", ylab="Tasa de verdaderos positivos")
abline(0,1,lty=3,col="grey60")
col <- hcl.colors(K,"Dark 3")
for (k in 1:K) lines(roc_folds[[k]]$fpr, roc_folds[[k]]$tpr, col=col[k], lwd=1.9)
title(main="Curvas ROC por pliegue espacial", line=2.0, cex.main=1.2)
mtext("Cada curva corresponde a una región del país retenida como evaluación", side=3, line=0.5,
      cex=0.78, col="grey35")
auc_k <- fread(file.path(RUTA,"salidas/tablas/tabla_comparacion_modelos.csv"))
legend("bottomright", bty="n", lwd=2, col=col, cex=0.82,
       legend=paste("Pliegue", 1:K))
mtext("AUC medio = 0,859 ± 0,052", side=1, line=-1.5, adj=0.95, cex=0.8, col="grey25")
dev.off(); cat("fig09 escrita\n")

# ---- Diseño de bloques espaciales ------------------------------
PROJ <- "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m"
peru <- project(geodata::gadm("PER", level=0, path=file.path(RUTA,"datos/crudos")), PROJ)
r_fold <- rast(g$tierra); values(r_fold) <- NA; r_fold[D$celda] <- D$fold
e <- ext(peru); dx <- xmax(e)-xmin(e); dy <- ymax(e)-ymin(e)
png(file.path(RUTA,"salidas/figuras/fig10_bloques.png"), width=1500, height=1900, res=210)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ext(xmin(e)-0.04*dx, xmax(e)+0.16*dx, ymin(e)-0.10*dy, ymax(e)+0.10*dy),
     col=NA, border=NA, axes=FALSE, xlab="", ylab="")
plot(r_fold, add=TRUE, col=hcl.colors(K,"Set 2"), legend=FALSE)
plot(peru, border="grey35", add=TRUE)
points(D[pres==1, .(x,y)], pch=16, cex=0.18, col="#00000075")
text(xmin(e)+0.5*dx, ymax(e)+0.07*dy, "Bloques espaciales de validación", cex=1.25, font=2)
text(xmin(e)+0.5*dx, ymax(e)+0.032*dy,
     "172 bloques de 100 km repartidos en cinco pliegues", cex=0.78, col="grey30")
legend(xmax(e)+0.02*dx, ymin(e)+0.55*dy, bty="n", pch=15, pt.cex=1.4, cex=0.78,
       col=hcl.colors(K,"Set 2"), legend=paste("Pliegue", 1:K))
dev.off(); cat("fig10 escrita\n")
