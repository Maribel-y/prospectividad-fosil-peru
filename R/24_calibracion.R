# ============================================================
# Calibration of the national GBM under spatial blocking.
# The model is fitted on presences plus a 10 000-cell background with
# class weights, so its raw output is not an absolute probability of
# occurrence. This script quantifies that on the evaluation folds:
#   raw       GBM output as is
#   platt     logistic recalibration fitted on the TRAINING cells of
#             each fold (logit of the raw score as the only predictor)
# and compares both with the no-skill reference (training prevalence),
# through the Brier score and a reliability curve by predicted decile.
# ============================================================
source(file.path(path.expand("~/Desktop/dinos"), "R/comun_modelado.R"))
K <- 5; REP <- 5
logit <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

filas <- list(); curvas <- list()
for (r in 1:REP) { set.seed(r)
  for (k in 1:K) {
    s <- dividir(D$fold == k); te <- s$te; y <- te$pres
    m <- gbm(F_ALL, s$ent)
    p_tr <- pred(m, F_ALL, s$tr); p_te <- pred(m, F_ALL, te)
    cal <- glm(s$tr$pres ~ logit(p_tr), family = binomial())
    p_pl <- as.numeric(plogis(coef(cal)[1] + coef(cal)[2] * logit(p_te)))
    prev <- mean(s$tr$pres)
    brier <- function(p) mean((p - y)^2)
    filas[[length(filas)+1]] <- data.table(rep = r, fold = k, prevalencia = mean(y),
      brier_ref = brier(rep(prev, length(y))), brier_raw = brier(p_te), brier_platt = brier(p_pl),
      media_raw = mean(p_te), media_platt = mean(p_pl), max_platt = max(p_pl))
    dec <- cut(rank(p_te, ties.method = "first"), 10, labels = FALSE)
    curvas[[length(curvas)+1]] <- data.table(rep = r, fold = k, decil = dec, y = y, raw = p_te, platt = p_pl)
  }
  cat("rep", r, "done\n")
}
R <- rbindlist(filas)
R[, `:=`(bss_raw = 1 - brier_raw/brier_ref, bss_platt = 1 - brier_platt/brier_ref)]
res <- R[, lapply(.SD, mean), .SDcols = -c("rep","fold")]
cat("\n===== CALIBRATION (", REP, "reps x", K, "spatial folds ) =====\n"); print(t(round(res, 5)))

C <- rbindlist(curvas)[, .(obs = mean(y), raw = mean(raw), platt = mean(platt), n = .N), by = decil][order(decil)]
cat("\nReliability by predicted decile (pooled):\n"); print(C[, lapply(.SD, function(x) signif(x, 3))])
# Top percentile, where prioritisation actually operates
top <- rbindlist(curvas)[, .(y, platt, raw, q = frank(raw)/.N), by = .(rep, fold)][q >= 0.99]
cat("\nTop 1 % of each fold: observed rate", round(mean(top$y), 4),
    "| mean recalibrated probability", round(mean(top$platt), 4), "\n")

fwrite(R, file.path(RUTA, "salidas/tablas/tabla_calibracion_reps.csv"))
fwrite(rbind(res[, .(indicador = names(res), valor = unlist(res))],
             data.table(indicador = c("top1_obs","top1_platt"), valor = c(mean(top$y), mean(top$platt)))),
       file.path(RUTA, "salidas/tablas/tabla_calibracion.csv"))
fwrite(C, file.path(RUTA, "salidas/tablas/tabla_calibracion_curva.csv"))

png(file.path(RUTA, "salidas/figuras/figS_calibracion.png"), width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
plot(C$raw, C$obs, log = "xy", pch = 19, col = "#4A6FA5", xlab = "Mean raw GBM score (decile)",
     ylab = "Observed presence rate", main = "(a) Raw score")
abline(0, 1, lty = 2, col = "grey50")
plot(C$platt, C$obs, log = "xy", pch = 19, col = "#C1440E", xlab = "Mean recalibrated probability (decile)",
     ylab = "Observed presence rate", main = "(b) Platt-recalibrated")
abline(0, 1, lty = 2, col = "grey50")
dev.off()
