# ============================================================
# Fase 4 (parte B) — Entrenamiento y comparación de modelos
#
# MaxEnt, Random Forest y Gradient Boosting frente al modelo base
# geológico, todos bajo los mismos bloques espaciales de 100 km.
# Dos diseños: nacional (fondo aleatorio + covariables de accesibilidad)
# y mesozoico (target-group background: las demás localidades fósiles
# replican el esfuerzo de colecta).
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(terra); library(ranger); library(xgboost); library(maxnet)
})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
K <- 5; N_FONDO <- 10000

COVAR <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion",
           "dist_vias_km","dist_poblados_km","dist_rios_km","woe_uni")

# ---- Métricas -------------------------------------------------
auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0)
                        (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
tss <- function(y, p) { u <- unique(quantile(p, seq(0.5,0.999,0.001), na.rm=TRUE))
  max(sapply(u, function(t) mean(p[y==1] >= t) + mean(p[y==0] < t) - 1)) }
boyce <- function(y, p, nclas = 10) {
  br <- unique(quantile(p, seq(0,1,length.out=nclas+1), na.rm=TRUE))
  if (length(br) < 3) return(NA_real_)
  br[1] <- -Inf; br[length(br)] <- Inf
  cl <- cut(p, br, labels=FALSE); k <- length(br)-1
  pe <- sapply(1:k, function(i){ E <- sum(cl==i)/length(cl)
                                 if (E==0) NA else (sum(y[cl==i]==1)/sum(y==1))/E })
  if (sum(!is.na(pe)) < 3) return(NA_real_)
  suppressWarnings(cor(1:k, pe, method="spearman", use="complete.obs"))
}

# WoE de unidad estimado solo en entrenamiento (evita fuga de información)
woe_map <- function(train, resp) {
  N <- nrow(train); Dn <- sum(train[[resp]])
  w <- train[, .(n = .N, d = sum(get(resp))), by = uni]
  w[, Wmas := log(((d + 0.5)/Dn)/((n - d + 0.5)/(N - Dn)))]
  setNames(w$Wmas, w$uni)
}
aplicar_woe <- function(mapa, uni) { v <- mapa[uni]; fifelse(is.na(v), 0, as.numeric(v)) }

X <- function(dt) model.matrix(~ era + amb + log_edad + edad_na + elevacion + pendiente +
                                 rugosidad + exposicion + dist_vias_km + dist_poblados_km + dist_rios_km +
                                 woe_uni - 1, data = dt)

# ---- Un pliegue, los cuatro modelos ---------------------------
correr_fold <- function(k, resp, diseno, excluir_sosp) {
  tr <- D[fold != k]; te <- D[fold == k]
  if (excluir_sosp) tr <- tr[!(get(resp) == 1 & unidad_ignea)]
  if (sum(te[[resp]]) < 3 || sum(tr[[resp]]) < 10) return(NULL)

  wm <- woe_map(tr, resp)
  tr[, woe_uni := aplicar_woe(wm, uni)]; te[, woe_uni := aplicar_woe(wm, uni)]

  presencias <- tr[get(resp) == 1]
  fondo <- if (diseno == "tgb") {
    # target-group background: otras localidades fósiles (mismo esfuerzo de colecta)
    rbind(tr[pres == 1 & get(resp) == 0], tr[pres == 0][sample(.N, min(.N, 2000))])
  } else {
    tr[get(resp) == 0][sample(.N, min(.N, N_FONDO))]
  }
  ent <- rbind(presencias, fondo); y <- ent[[resp]]

  pred <- list()
  # Modelo base geológico (referencia)
  mb <- suppressWarnings(glm(as.formula(paste(resp, "~ era + amb + woe_uni")),
                             data = tr, family = binomial))
  pred[["Base geológico"]] <- suppressWarnings(predict(mb, newdata = te, type = "response"))
  # MaxEnt
  pred[["MaxEnt"]] <- tryCatch({
    mx <- maxnet(y, as.data.frame(ent[, ..COVAR]),
                 f = maxnet.formula(y, as.data.frame(ent[, ..COVAR]), classes = "lq"),
                 regmult = 2)
    as.numeric(predict(mx, as.data.frame(te[, ..COVAR]), type = "cloglog"))
  }, error = function(e) { cat("   MaxEnt falló:", conditionMessage(e), "\n"); rep(NA, nrow(te)) })
  # Random Forest
  rf <- ranger(x = as.data.frame(ent[, ..COVAR]), y = factor(y), probability = TRUE,
               num.trees = 500, min.node.size = 5, num.threads = 8,
               case.weights = ifelse(y == 1, sum(y == 0)/sum(y == 1), 1))
  pred[["Random Forest"]] <- predict(rf, as.data.frame(te[, ..COVAR]), num.threads = 8)$predictions[, "1"]
  # Gradient Boosting
  # xgboost 3.x cambió la interfaz de xgboost(); xgb.train sobre DMatrix es estable
  dtr <- xgb.DMatrix(data = X(ent), label = y)
  gb <- xgb.train(params = list(objective = "binary:logistic", max_depth = 5,
                                eta = 0.05, subsample = 0.8, colsample_bytree = 0.8,
                                scale_pos_weight = sum(y == 0)/sum(y == 1), nthread = 8),
                  data = dtr, nrounds = 300, verbose = 0)
  pred[["Gradient Boosting"]] <- predict(gb, xgb.DMatrix(X(te)))

  rbindlist(lapply(names(pred), function(m) {
    p <- pred[[m]]; ok <- !is.na(p)
    if (sum(ok) < 10) return(NULL)
    data.table(modelo = m, fold = k, AUC = auc(te[[resp]][ok], p[ok]),
               TSS = tss(te[[resp]][ok], p[ok]), Boyce = boyce(te[[resp]][ok], p[ok]))
  }))
}

comparar <- function(resp, diseno, etiqueta, excluir_sosp = TRUE) {
  cat("\n=====", etiqueta, "=====\n")
  cat("presencias:", sum(D[[resp]]), "| diseño de fondo:", diseno,
      "| sospechosas excluidas:", excluir_sosp, "\n")
  r <- rbindlist(lapply(1:K, function(k) correr_fold(k, resp, diseno, excluir_sosp)))
  tab <- r[, .(AUC = mean(AUC), AUC_sd = sd(AUC), TSS = mean(TSS),
               Boyce = mean(Boyce, na.rm = TRUE)), by = modelo][order(-AUC)]
  print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
  tab[, escenario := etiqueta]
  tab
}

res <- rbindlist(list(
  comparar("pres",     "aleatorio", "MODELO NACIONAL (todas las localidades)"),
  comparar("pres_mes", "tgb",       "MODELO ESTRATIFICADO AL MESOZOICO"),
  comparar("pres",     "aleatorio", "NACIONAL — sensibilidad: con las 15 sospechosas", FALSE)
))
fwrite(res, file.path(RUTA,"salidas/tablas/tabla_comparacion_modelos.csv"))
cat("\nComparación escrita.\n")
