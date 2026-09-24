# ============================================================
# Contraste directo de las hipótesis H2 y H3 del anteproyecto
#
# H2: el ML supera al modelo base solo si se controla el sesgo de
#     accesibilidad  -> ablación de las covariables de accesibilidad.
# H3: la validación cruzada aleatoria sobrestima el desempeño frente
#     al bloqueo espacial -> mismos modelos, pliegues aleatorios.
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(terra); library(ranger); library(xgboost); library(maxnet)
})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
K <- 5; N_FONDO <- 10000
ACCESO <- c("dist_vias_km","dist_poblados_km","dist_rios_km")
BASE   <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion","woe_uni")

auc <- function(y,p){ r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1])-n1*(n1+1)/2)/(n1*n0) }
tss <- function(y,p){ u <- unique(quantile(p, seq(0.5,0.999,0.001), na.rm=TRUE))
  max(sapply(u, function(t) mean(p[y==1]>=t)+mean(p[y==0]<t)-1)) }

woe_map <- function(train){ N <- nrow(train); Dn <- sum(train$pres)
  w <- train[, .(n=.N, d=sum(pres)), by=uni]
  w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]; setNames(w$Wmas, w$uni) }

correr <- function(k, covar, columna_fold) {
  tr <- D[get(columna_fold) != k]; te <- D[get(columna_fold) == k]
  tr <- tr[!(pres == 1 & unidad_ignea)]
  if (sum(te$pres) < 3) return(NULL)
  wm <- woe_map(tr)
  tr[, woe_uni := fifelse(is.na(wm[uni]), 0, as.numeric(wm[uni]))]
  te[, woe_uni := fifelse(is.na(wm[uni]), 0, as.numeric(wm[uni]))]
  ent <- rbind(tr[pres==1], tr[pres==0][sample(.N, min(.N, N_FONDO))]); y <- ent$pres
  f <- as.formula(paste("~", paste(covar, collapse=" + "), "- 1"))
  Xm <- function(dt) model.matrix(f, data = dt)

  rf <- ranger(x = as.data.frame(ent[, ..covar]), y = factor(y), probability = TRUE,
               num.trees = 500, min.node.size = 5, num.threads = 8,
               case.weights = ifelse(y==1, sum(y==0)/sum(y==1), 1))
  p_rf <- predict(rf, as.data.frame(te[, ..covar]), num.threads = 8)$predictions[, "1"]
  dtr <- xgb.DMatrix(data = Xm(ent), label = y)
  gb <- xgb.train(params = list(objective="binary:logistic", max_depth=5, eta=0.05,
                                subsample=0.8, colsample_bytree=0.8,
                                scale_pos_weight=sum(y==0)/sum(y==1), nthread=8),
                  data = dtr, nrounds = 300, verbose = 0)
  p_gb <- predict(gb, xgb.DMatrix(Xm(te)))
  data.table(fold = k,
             RF_AUC = auc(te$pres, p_rf), RF_TSS = tss(te$pres, p_rf),
             GB_AUC = auc(te$pres, p_gb), GB_TSS = tss(te$pres, p_gb))
}

resumen <- function(covar, columna_fold, etiqueta) {
  r <- rbindlist(lapply(1:K, function(k) correr(k, covar, columna_fold)))
  data.table(escenario = etiqueta,
             RF_AUC = round(mean(r$RF_AUC),3), RF_TSS = round(mean(r$RF_TSS),3),
             GB_AUC = round(mean(r$GB_AUC),3), GB_TSS = round(mean(r$GB_TSS),3))
}

# H3 necesita pliegues aleatorios (sin estructura espacial)
D[, fold_aleatorio := sample(rep_len(1:K, .N))]

cat("===== H2: ¿aporta el control explícito del sesgo de accesibilidad? =====\n")
h2 <- rbindlist(list(
  resumen(BASE,                "fold", "Sin covariables de accesibilidad"),
  resumen(c(BASE, ACCESO),     "fold", "Con covariables de accesibilidad")))
print(h2)

cat("\n===== H3: validación aleatoria frente a bloqueo espacial =====\n")
h3 <- rbindlist(list(
  resumen(c(BASE, ACCESO), "fold",           "Bloqueo espacial (100 km)"),
  resumen(c(BASE, ACCESO), "fold_aleatorio", "k-fold aleatorio")))
print(h3)
cat("\nSobrestimación del AUC por usar validación aleatoria:",
    round(h3[2, GB_AUC] - h3[1, GB_AUC], 3), "(Gradient Boosting)\n")

fwrite(rbind(h2, h3), file.path(RUTA,"salidas/tablas/tabla_hipotesis.csv"))
