# Priorización robusta: celdas que permanecen en el 1 % superior sea cual
# sea la constante empleada para neutralizar el sesgo de accesibilidad.
suppressPackageStartupMessages({library(data.table); library(terra); library(xgboost)})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
ACC <- c("dist_vias_km","dist_poblados_km","dist_rios_km")
COVAR <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion", ACC, "woe_uni")
FORM <- as.formula(paste("~", paste(COVAR, collapse=" + "), "- 1"))
woe_map <- function(tr){ N<-nrow(tr); Dn<-sum(tr$pres)
  w <- tr[, .(n=.N,d=sum(pres)), by=uni]; w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]
  setNames(w$Wmas,w$uni) }
apl <- function(m,u){ v<-m[u]; fifelse(is.na(v),0,as.numeric(v)) }

fin <- readRDS(file.path(RUTA, "datos/procesados/modelo_final_nacional.rds"))   # script 10
Df <- copy(D)[, woe_uni := apl(fin$woe, uni)]
m <- fin$modelo

cte <- list(cero=function(v) 0,
            p25 =function(v) quantile(D[[v]][D$pres==1],0.25,na.rm=TRUE),
            p50 =function(v) median(D[[v]][D$pres==1],na.rm=TRUE),
            p75 =function(v) quantile(D[[v]][D$pres==1],0.75,na.rm=TRUE),
            terr=function(v) median(D[[v]],na.rm=TRUE))
P <- sapply(cte, function(f){ Dc <- copy(Df); for (v in ACC) Dc[[v]] <- f(v)
                              predict(m, xgb.DMatrix(model.matrix(FORM, Dc))) })
enTop <- apply(P, 2, function(p) p >= quantile(p, 0.99))
D[, veces_en_top := rowSums(enTop)]
D[, prosp_media  := rowMeans(P)]

cat("===== PRIORIZACIÓN ROBUSTA A LA CONSTANTE =====\n")
t <- D[, .(celdas = .N), by = veces_en_top][order(-veces_en_top)]
t[, km2 := celdas*4]; print(t)
rob <- D[veces_en_top == 5]
cat("\nCeldas en el 1 % superior bajo las cinco constantes:", nrow(rob),
    sprintf("(%s km², %.2f %% del país)\n", format(nrow(rob)*4, big.mark=" "), 100*nrow(rob)*4/1285216))
cat("Presencias conocidas que caen en ellas:", sum(rob$pres), "\n")
cat("Frente a las", sum(enTop[, "p50"]), "celdas del 1 % con la constante adoptada.\n")

r <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))$tierra
cap <- function(v,n){ x <- rast(r); values(x) <- NA; x[D$celda] <- v; names(x) <- n; x }
writeRaster(c(cap(D$prosp_media,"prospectividad_media"), cap(D$veces_en_top,"robustez")),
            file.path(RUTA,"datos/procesados/prospectividad_robusta.tif"), overwrite=TRUE)
fwrite(rob[, .(celda,x,y,era,amb,uni,prosp_media)],
       file.path(RUTA,"salidas/tablas/tabla_priorizacion_robusta.csv"))
