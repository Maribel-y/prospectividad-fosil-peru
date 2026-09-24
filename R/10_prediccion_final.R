# ============================================================
# Fase 5 — Predicción, incertidumbre y priorización
#
# El mapa de prospectividad se obtiene con el modelo ganador
# (Gradient Boosting). La incertidumbre es la desviación estándar
# entre los cinco modelos entrenados en la validación con bloqueo
# espacial: mide cuánto cambia la predicción según qué región del
# país se haya usado para entrenar.
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(terra); library(xgboost)
})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")
D <- fread(file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
D[, era := factor(era)]; D[, amb := factor(amb)]
g <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
K <- 5; N_FONDO <- 10000

COVAR <- c("era","amb","log_edad","edad_na","elevacion","pendiente","rugosidad","exposicion",
           "dist_vias_km","dist_poblados_km","dist_rios_km","woe_uni")
FORM <- as.formula(paste("~", paste(COVAR, collapse=" + "), "- 1"))

woe_map <- function(train, resp) { N <- nrow(train); Dn <- sum(train[[resp]])
  w <- train[, .(n=.N, d=sum(get(resp))), by=uni]
  w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]; setNames(w$Wmas, w$uni) }
aplicar <- function(m, u) { v <- m[u]; fifelse(is.na(v), 0, as.numeric(v)) }

entrenar_gb <- function(tr, resp, diseno) {
  presencias <- tr[get(resp) == 1]
  fondo <- if (diseno == "tgb")
    rbind(tr[pres == 1 & get(resp) == 0], tr[pres == 0][sample(.N, min(.N, 2000))])
  else tr[get(resp) == 0][sample(.N, min(.N, N_FONDO))]
  ent <- rbind(presencias, fondo); y <- ent[[resp]]
  xgb.train(params = list(objective="binary:logistic", max_depth=5, eta=0.05,
                          subsample=0.8, colsample_bytree=0.8,
                          scale_pos_weight=sum(y==0)/sum(y==1), nthread=8),
            data = xgb.DMatrix(model.matrix(FORM, ent), label = y),
            nrounds = 300, verbose = 0)
}

producir <- function(resp, diseno, sufijo, titulo) {
  cat("\n=====", titulo, "=====\n")
  # Predicciones de los cinco modelos de la validación espacial
  P <- matrix(NA_real_, nrow = nrow(D), ncol = K)
  for (k in 1:K) {
    tr <- D[fold != k][!(get(resp) == 1 & unidad_ignea)]
    wm <- woe_map(tr, resp)
    tr[, woe_uni := aplicar(wm, uni)]
    Dk <- copy(D)[, woe_uni := aplicar(wm, uni)]
    P[, k] <- predict(entrenar_gb(tr, resp, diseno), xgb.DMatrix(model.matrix(FORM, Dk)))
    cat("  pliegue", k, "listo\n")
  }
  # Modelo final con todos los datos
  tr <- D[!(get(resp) == 1 & unidad_ignea)]
  wm <- woe_map(tr, resp); tr[, woe_uni := aplicar(wm, uni)]
  Df <- copy(D)[, woe_uni := aplicar(wm, uni)]
  mfin <- entrenar_gb(tr, resp, diseno)
  # Scripts 15, 16 and 21 reuse this exact model, so every table and the
  # robust priority set describe the same fitted map
  saveRDS(list(modelo = mfin, woe = wm), file.path(RUTA, paste0("datos/procesados/modelo_final_", sufijo, ".rds")))
  pred <- predict(mfin, xgb.DMatrix(model.matrix(FORM, Df)))

  # Predicción con el sesgo de muestreo apagado: las covariables de
  # accesibilidad se fijan al valor típico de las celdas donde SÍ se ha
  # buscado, de modo que el mapa responda a la geología y no a la
  # facilidad de acceso. Es el mapa útil para priorizar campañas.
  Dc <- copy(Df)
  for (v in c("dist_vias_km","dist_poblados_km","dist_rios_km"))
    Dc[[v]] <- median(D[[v]][D[[resp]] == 1], na.rm = TRUE)
  pred_corr <- predict(mfin, xgb.DMatrix(model.matrix(FORM, Dc)))

  incert <- apply(P, 1, sd)
  D[, (paste0("prosp_", sufijo))  := pred]
  D[, (paste0("incert_", sufijo)) := incert]
  D[, (paste0("prospcorr_", sufijo)) := pred_corr]

  cap <- function(v, nm) { r <- rast(g$tierra); values(r) <- NA; r[D$celda] <- v; names(r) <- nm; r }
  writeRaster(c(cap(pred, "prospectividad"), cap(incert, "incertidumbre"),
                cap(pred_corr, "prospectividad_sin_sesgo")),
              file.path(RUTA, paste0("datos/procesados/prospectividad_", sufijo, ".tif")),
              overwrite = TRUE)

  # Priorización: alto potencial y baja incertidumbre
  # Las celdas más prospectivas son también las más inciertas, así que el
  # corte de incertidumbre se toma DENTRO del 1 % superior, no sobre todo el país.
  alto <- which(pred_corr >= quantile(pred_corr, 0.99))
  u_inc <- quantile(incert[alto], 0.25)
  prio <- D[alto][incert[alto] <= u_inc]
  cat("  celdas en el 1 % superior:", length(alto),
      "| priorizadas (cuartil menos incierto):", nrow(prio), "\n")
  cat("  prospectividad: mediana", round(median(pred),4), "| máximo", round(max(pred),3), "\n")
  cat("  incertidumbre : mediana", round(median(incert),4), "| máximo", round(max(incert),3), "\n")
  prio
}

prio_nac <- producir("pres",     "aleatorio", "nacional",  "PROSPECTIVIDAD NACIONAL")
prio_mes <- producir("pres_mes", "tgb",       "mesozoico", "PROSPECTIVIDAD MESOZOICA")

fwrite(D[, .(celda, x, y, pres, pres_mes, era, amb, uni,
             prosp_nacional, prospcorr_nacional, incert_nacional,
             prosp_mesozoico, prospcorr_mesozoico, incert_mesozoico)],
       file.path(RUTA,"datos/procesados/predicciones.csv"))
fwrite(prio_nac[, .(celda, x, y, era, amb, uni, prospcorr_nacional, incert_nacional)],
       file.path(RUTA,"salidas/tablas/tabla_priorizacion_nacional.csv"))
fwrite(prio_mes[, .(celda, x, y, era, amb, uni, prospcorr_mesozoico, incert_mesozoico)],
       file.path(RUTA,"salidas/tablas/tabla_priorizacion_mesozoico.csv"))
cat("\nMapas de prospectividad e incertidumbre escritos.\n")
