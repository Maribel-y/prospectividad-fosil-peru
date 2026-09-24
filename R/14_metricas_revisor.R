# ============================================================
# Model comparison repeated with REP seeds (background sampling and
# model randomness) over the five spatial folds of 100 km, for both
# designs. Replaces the single-seed comparison and the paired t-test
# over five folds: differences are reported as the mean over
# repetitions, a 95 % interval across repetitions, and the number of
# the REP x K folds in which the difference is positive.
# ============================================================
source(file.path(path.expand("~/Desktop/dinos"), "R/comun_modelado.R"))
suppressPackageStartupMessages({library(ranger); library(maxnet)})
K <- 5; REP <- 10
V_ALL <- c(GEO, TOPO, EXPO, ACC)

ajustar <- function(tr, te, resp, diseno) {
  wm <- { N <- nrow(tr); Dn <- sum(tr[[resp]])
          w <- tr[, .(n = .N, d = sum(get(resp))), by = uni]
          w[, Wmas := log(((d+0.5)/Dn)/((n-d+0.5)/(N-Dn)))]; setNames(w$Wmas, w$uni) }
  tr[, woe_uni := apl(wm, uni)]; te[, woe_uni := apl(wm, uni)]
  fondo <- if (diseno == "tgb") rbind(tr[pres == 1 & get(resp) == 0], tr[pres == 0][sample(.N, min(.N, 2000))])
           else tr[get(resp) == 0][sample(.N, min(.N, N_FONDO))]
  ent <- rbind(tr[get(resp) == 1], fondo); y <- ent[[resp]]; set(ent, j = "pres", value = y)
  Xe <- as.data.frame(ent[, ..V_ALL]); Xt <- as.data.frame(te[, ..V_ALL])
  list(
    "Geological baseline" = { t2 <- copy(tr); set(t2, j = "pres", value = t2[[resp]]); base_glm(t2, te) },
    "MaxEnt" = as.numeric(predict(maxnet(y, Xe, f = maxnet.formula(y, Xe, classes = "lq"), regmult = 2),
                                  Xt, type = "cloglog")),
    "Random forest" = predict(ranger(x = Xe, y = factor(y), probability = TRUE, num.trees = 500,
                                     min.node.size = 5, num.threads = 8,
                                     case.weights = ifelse(y == 1, sum(y == 0)/sum(y == 1), 1)),
                              Xt, num.threads = 8)$predictions[, "1"],
    "Gradient boosting" = pred(gbm(F_ALL, ent), F_ALL, te))
}

filas <- list()
for (dis in list(c("pres", "aleatorio", "National"), c("pres_mes", "tgb", "Mesozoic")))
  for (r in 1:REP) { set.seed(r)
    for (k in 1:K) {
      resp <- dis[1]
      tr <- D[fold != k][!(get(resp) == 1 & unidad_ignea)]; te <- D[fold == k]
      if (sum(te[[resp]]) < 3) next
      p <- ajustar(tr, te, resp, dis[2]); y <- te[[resp]]
      filas[[length(filas)+1]] <- rbindlist(lapply(names(p), function(m) data.table(
        design = dis[3], rep = r, fold = k, model = m, AUC = auc(y, p[[m]]), TSS = tss(y, p[[m]]),
        PR_AUC = pr_auc(y, p[[m]]), Boyce = boyce(y, p[[m]]))))
    }
    cat(dis[3], "rep", r, "done\n")
  }
R <- rbindlist(filas)
fwrite(R, file.path(RUTA, "salidas/tablas/tabla_metricas_por_pliegue.csv"))

# ---- Table of models: mean over repetitions, sd between repetitions
porrep <- R[, lapply(.SD, mean, na.rm = TRUE), by = .(design, model, rep), .SDcols = c("AUC","TSS","PR_AUC","Boyce")]
tab <- porrep[, .(AUC = mean(AUC), AUC_sd_rep = sd(AUC), TSS = mean(TSS), PR_AUC = mean(PR_AUC),
                  Boyce = mean(Boyce)), by = .(design, model)]
fold_sd <- R[, .(AUC_sd_fold = sd(AUC)), by = .(design, model, rep)][, .(AUC_sd_fold = mean(AUC_sd_fold)), by = .(design, model)]
tab <- merge(tab, fold_sd, by = c("design","model"))[order(design, -AUC)]
cat("\n===== MODELS (", REP, "reps x", K, "spatial folds ) =====\n")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(tab, file.path(RUTA, "salidas/tablas/tabla_comparacion_modelos_rep.csv"))

# ---- Paired differences
W <- dcast(R, design + rep + fold ~ model, value.var = c("AUC","PR_AUC"))
dif <- function(des, a, b, met, nombre) {
  x <- W[design == des]; d <- x[[paste0(met, "_", a)]] - x[[paste0(met, "_", b)]]
  dr <- tapply(d, x$rep, mean)
  data.table(design = des, comparison = nombre, mean_diff = mean(dr),
             lo95 = quantile(dr, 0.025), hi95 = quantile(dr, 0.975),
             fold_min = min(d), fold_max = max(d), folds_positive = sprintf("%d/%d", sum(d > 0), length(d)))
}
comp <- rbindlist(list(
  dif("National", "Gradient boosting", "Geological baseline", "AUC",    "AUC: boosting - baseline"),
  dif("National", "Gradient boosting", "MaxEnt",              "AUC",    "AUC: boosting - MaxEnt"),
  dif("National", "Gradient boosting", "Random forest",       "AUC",    "AUC: boosting - random forest"),
  dif("National", "Gradient boosting", "Geological baseline", "PR_AUC", "PR-AUC: boosting - baseline"),
  dif("Mesozoic", "Gradient boosting", "Geological baseline", "AUC",    "AUC: boosting - baseline")))
cat("\n===== PAIRED DIFFERENCES =====\n")
print(comp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
fwrite(comp, file.path(RUTA, "salidas/tablas/tabla_test_pareado.csv"))
