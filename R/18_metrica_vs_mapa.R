# ============================================================
# Métrica frente a mapa, con el MISMO modelo y los MISMOS pliegues
#
# El 41,4 % de cambio en el decil prioritario compara dos proyecciones
# del mismo modelo (observacional y neutralizada), mientras que el 0,042
# de AUC procede de la ablación (modelo con y sin covariables de acceso).
# Este guion evalúa las tres variantes sobre idénticos pliegues y fondo:
#   obs : modelo con acceso, proyectado con los valores reales de acceso
#   neu : el mismo modelo, con el acceso fijado a la mediana de las
#         presencias de ENTRENAMIENTO (sin fuga hacia el pliegue de prueba)
#   sin : modelo ajustado sin covariables de acceso (ablación)
# y añade la validación aleatoria del modelo con acceso (H3).
# Se repite con REP semillas para separar el efecto del ruido de muestreo
# del fondo, que es lo que explica las pequeñas diferencias entre las
# cifras de los guiones 08 y 09 (0,859 / 0,857 / 0,855).
# ============================================================
suppressPackageStartupMessages({library(data.table); library(xgboost)})
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
K <- 5; N_FONDO <- 10000; REP <- 10
ACC  <- c("dist_vias_km","dist_poblados_km","dist_rios_km")
BASE <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion","woe_uni")
F_CON <- as.formula(paste("~", paste(c(BASE, ACC), collapse=" + "), "- 1"))
F_SIN <- as.formula(paste("~", paste(BASE, collapse=" + "), "- 1"))

auc <- function(y,p){ r<-rank(p); n1<-sum(y==1); n0<-sum(y==0); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
tss <- function(y,p){ u <- unique(quantile(p, seq(0.5,0.999,0.001), na.rm=TRUE))
  max(sapply(u, function(t) mean(p[y==1]>=t)+mean(p[y==0]<t)-1)) }
pr_auc <- function(y,p){ o<-order(-p); y<-y[o]; tp<-cumsum(y==1); fp<-cumsum(y==0)
  sum(diff(c(0, tp/sum(y==1))) * tp/(tp+fp)) }
woe_map <- function(tr){ N<-nrow(tr); Dn<-sum(tr$pres)
  w <- tr[, .(n=.N,d=sum(pres)), by=uni]; w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]
  setNames(w$Wmas,w$uni) }
apl <- function(m,u){ v<-m[u]; fifelse(is.na(v),0,as.numeric(v)) }
gbm <- function(form, ent){ y <- ent$pres
  xgb.train(params=list(objective="binary:logistic",max_depth=5,eta=0.05,subsample=0.8,
                        colsample_bytree=0.8,scale_pos_weight=sum(y==0)/sum(y==1),nthread=8),
            data=xgb.DMatrix(model.matrix(form,ent),label=y), nrounds=300, verbose=0) }
prep <- function(col, k){
  tr <- D[get(col)!=k][!(pres==1 & unidad_ignea)]; te <- D[get(col)==k]
  wm <- woe_map(tr); tr[, woe_uni := apl(wm,uni)]; te[, woe_uni := apl(wm,uni)]
  list(tr=tr, te=te, ent=rbind(tr[pres==1], tr[pres==0][sample(.N,min(.N,N_FONDO))])) }
met <- function(y,p,pref) setNames(list(auc(y,p), tss(y,p), pr_auc(y,p)), paste0(pref, c("_AUC","_TSS","_PR")))

filas <- list()
for (r in 1:REP) {
  set.seed(r)
  D[, fold_aleatorio := sample(rep_len(1:K, .N))]
  for (k in 1:K) {
    s <- prep("fold", k); te <- s$te; y <- te$pres
    m_con <- gbm(F_CON, s$ent); m_sin <- gbm(F_SIN, s$ent)   # mismo fondo para ambos
    p_obs <- predict(m_con, xgb.DMatrix(model.matrix(F_CON, te)))
    te_n <- copy(te); for (v in ACC) te_n[[v]] <- median(s$tr[[v]][s$tr$pres==1], na.rm=TRUE)
    p_neu <- predict(m_con, xgb.DMatrix(model.matrix(F_CON, te_n)))
    p_sin <- predict(m_sin, xgb.DMatrix(model.matrix(F_SIN, te)))
    top_obs <- p_obs >= quantile(p_obs, 0.9); top_neu <- p_neu >= quantile(p_neu, 0.9)
    sa <- prep("fold_aleatorio", k)
    p_ale <- predict(gbm(F_CON, sa$ent), xgb.DMatrix(model.matrix(F_CON, sa$te)))
    filas[[length(filas)+1]] <- as.data.table(c(list(rep=r, fold=k),
      met(y,p_obs,"obs"), met(y,p_neu,"neu"), met(y,p_sin,"sin"),
      list(ale_AUC = auc(sa$te$pres,p_ale), ale_TSS = tss(sa$te$pres,p_ale),
           cambio_decil = 100*(1 - sum(top_obs & top_neu)/sum(top_obs)),
           cor_vias_obs = cor(p_obs, te$dist_vias_km, method="spearman"),
           cor_vias_neu = cor(p_neu, te$dist_vias_km, method="spearman"))))
  }
  cat("repetición", r, "de", REP, "lista\n")
}
R <- rbindlist(filas)
fwrite(R, file.path(RUTA,"salidas/tablas/tabla_metrica_vs_mapa_reps.csv"))

# media por repetición (sobre pliegues) y luego media y sd entre repeticiones
porrep <- R[, lapply(.SD, mean), by=rep, .SDcols=-c("rep","fold")]
porrep[, `:=`(dif_neu_obs = neu_AUC - obs_AUC, dif_sin_obs = sin_AUC - obs_AUC,
              dif_ale_obs = ale_AUC - obs_AUC, difTSS_ale_obs = ale_TSS - obs_TSS)]
vars <- setdiff(names(porrep), "rep")
res <- data.table(indicador = vars,
                  media = sapply(vars, function(v) mean(porrep[[v]])),
                  sd_entre_rep = sapply(vars, function(v) sd(porrep[[v]])),
                  min = sapply(vars, function(v) min(porrep[[v]])),
                  max = sapply(vars, function(v) max(porrep[[v]])))
cat("\n===== MÉTRICA FRENTE A MAPA (", REP, "repeticiones x", K, "pliegues espaciales ) =====\n")
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x,4) else x)])
# diferencia pareada obs vs neu dentro de cada pliegue y repetición
d <- R$neu_AUC - R$obs_AUC
cat("\nAUC neutralizado - observacional, pareado por pliegue (n =", length(d), "): media",
    round(mean(d),4), "| rango", round(min(d),4), "a", round(max(d),4), "\n")
fwrite(res, file.path(RUTA,"salidas/tablas/tabla_metrica_vs_mapa.csv"))
