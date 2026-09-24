# ============================================================
# Riesgos 4 y 5 de la autorrevisión
#  (4) ¿Depende el mapa del valor al que se fijan las covariables
#      de accesibilidad al neutralizar el sesgo?
#  (5) ¿Es estable entre pliegues el peso de evidencia por unidad?
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(xgboost)})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
K <- 5; N_FONDO <- 10000
ACC <- c("dist_vias_km","dist_poblados_km","dist_rios_km")
COVAR <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion", ACC, "woe_uni")
FORM <- as.formula(paste("~", paste(COVAR, collapse=" + "), "- 1"))
woe_map <- function(tr){ N<-nrow(tr); Dn<-sum(tr$pres)
  w <- tr[, .(n=.N,d=sum(pres)), by=uni]; w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]
  setNames(w$Wmas,w$uni) }
apl <- function(m,u){ v<-m[u]; fifelse(is.na(v),0,as.numeric(v)) }

# ---------- (4) Sensibilidad a la constante ----------
fin <- readRDS(file.path(RUTA, "datos/procesados/modelo_final_nacional.rds"))   # script 10
Df <- copy(D)[, woe_uni := apl(fin$woe, uni)]
mfin <- fin$modelo

escenarios <- list(
  "Cero (acceso inmediato)"        = function(v) 0,
  "Percentil 25 de las presencias" = function(v) quantile(D[[v]][D$pres==1], 0.25, na.rm=TRUE),
  "Mediana de las presencias"      = function(v) median(D[[v]][D$pres==1], na.rm=TRUE),
  "Percentil 75 de las presencias" = function(v) quantile(D[[v]][D$pres==1], 0.75, na.rm=TRUE),
  "Mediana del territorio"         = function(v) median(D[[v]], na.rm=TRUE))

pred <- lapply(escenarios, function(f){
  Dc <- copy(Df); for (v in ACC) Dc[[v]] <- f(v)
  predict(mfin, xgb.DMatrix(model.matrix(FORM, Dc)))
})
ref <- pred[["Mediana de las presencias"]]
top <- function(p) which(p >= quantile(p, 0.99))
tab4 <- rbindlist(lapply(names(pred), function(n){
  p <- pred[[n]]
  data.table(escenario = n,
    valor_vias_km = round(escenarios[[n]]("dist_vias_km"), 1),
    cor_con_adoptado = round(cor(p, ref, method="spearman"), 3),
    solapamiento_1pct = round(100*length(intersect(top(p), top(ref)))/length(top(ref)), 1),
    cor_con_dist_vias = round(cor(p, D$dist_vias_km, method="spearman"), 3))
}))
cat("===== (4) SENSIBILIDAD A LA CONSTANTE DE NEUTRALIZACIÓN =====\n"); print(tab4)
p_obs <- predict(mfin, xgb.DMatrix(model.matrix(FORM, Df)))
cat("\nReferencia sin neutralizar: correlación con distancia a vías =",
    round(cor(p_obs, D$dist_vias_km, method="spearman"), 3), "\n")

# ---------- (5) Estabilidad del WoE por unidad ----------
mapas <- lapply(1:K, function(k) woe_map(D[fold!=k][!(pres==1 & unidad_ignea)]))
comunes <- Reduce(intersect, lapply(mapas, names))
M <- sapply(mapas, function(m) m[comunes])
rownames(M) <- comunes

npres <- D[pres==1, .N, by=uni]
info <- data.table(uni = comunes, n_pres = npres$N[match(comunes, npres$uni)])
info[is.na(n_pres), n_pres := 0]
sel <- which(info$n_pres >= 5)

pares <- combn(K, 2)
cor_todas <- mean(apply(pares, 2, function(i) cor(M[,i[1]], M[,i[2]], method="spearman")))
cor_sel   <- mean(apply(pares, 2, function(i) cor(M[sel,i[1]], M[sel,i[2]], method="spearman")))
signo <- apply(M, 1, function(x) all(x > 0) || all(x < 0))

cat("\n===== (5) ESTABILIDAD DEL WoE POR UNIDAD ENTRE PLIEGUES =====\n")
cat("  unidades presentes en los cinco pliegues:", length(comunes),
    "| de ellas con >= 5 presencias:", length(sel), "\n")
cat("  correlación media entre pliegues (todas):", round(cor_todas,3), "\n")
cat("  correlación media entre pliegues (>= 5 presencias):", round(cor_sel,3), "\n")
cat("  unidades con signo consistente en los cinco pliegues:",
    round(100*mean(signo),1), "% (todas) |",
    round(100*mean(signo[sel]),1), "% (>= 5 presencias)\n")

est <- data.table(unidad = comunes[sel], n_pres = info$n_pres[sel],
                  woe_medio = round(rowMeans(M[sel,,drop=FALSE]),2),
                  woe_min = round(apply(M[sel,,drop=FALSE],1,min),2),
                  woe_max = round(apply(M[sel,,drop=FALSE],1,max),2))
est[, amplitud := round(woe_max - woe_min, 2)]
setorder(est, -woe_medio)
cat("\n  Doce unidades más favorables y su variación entre pliegues:\n")
print(head(est, 12))
cat("\n  amplitud mediana entre pliegues:", round(median(est$amplitud),2),
    "| unidades con amplitud > 2:", sum(est$amplitud > 2), "de", nrow(est), "\n")

fwrite(tab4, file.path(RUTA,"salidas/tablas/tabla_sensibilidad_constante.csv"))
fwrite(est,  file.path(RUTA,"salidas/tablas/tabla_estabilidad_woe.csv"))
