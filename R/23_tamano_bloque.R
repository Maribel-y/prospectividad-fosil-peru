# ============================================================
# Sensitivity to the spatial block size (50, 100, 150, 200 km).
# For each size the blocks are re-drawn and allocated to five folds
# with REP different seeds. Checks whether the three conclusions of the
# paper hold whatever the block size:
#   (1) GBM outperforms the geological baseline
#   (2) random folds overestimate relative to spatial blocking
#   (3) neutralising access lowers the validated AUC and moves the map
# ============================================================
source(file.path(path.expand("~/Desktop/dinos"), "R/comun_modelado.R"))
K <- 5; REP <- 5
LADOS <- c(50, 100, 150, 200) * 1000

asignar <- function(lado, semilla) {
  set.seed(semilla)
  b <- paste0(floor(D$x/lado), "_", floor(D$y/lado)); bl <- unique(b)
  list(fold = sample(rep_len(1:K, length(bl)))[match(b, bl)], n_bloques = length(bl))
}

filas <- list()
for (lado in LADOS) for (r in 1:REP) {
  a <- asignar(lado, r); set.seed(r)
  for (k in 1:K) {
    s <- dividir(a$fold == k); te <- s$te; y <- te$pres
    if (sum(y) < 3) next
    m <- gbm(F_ALL, s$ent)
    p_obs <- pred(m, F_ALL, te); p_neu <- pred(m, F_ALL, neutralizar(te, s$tr))
    filas[[length(filas)+1]] <- data.table(lado_km = lado/1000, bloques = a$n_bloques, rep = r, fold = k,
      n_pres_test = sum(y), base_AUC = auc(y, base_glm(s$tr, te)),
      obs_AUC = auc(y, p_obs), neu_AUC = auc(y, p_neu),
      obs_TSS = tss(y, p_obs), obs_PR = pr_auc(y, p_obs),
      cambio_decil = 100*(1 - sum(p_obs >= quantile(p_obs,.9) & p_neu >= quantile(p_neu,.9)) /
                                sum(p_obs >= quantile(p_obs,.9))))
  }
  cat("block", lado/1000, "km - rep", r, "done\n")
}
R <- rbindlist(filas)

# Random-fold reference (independent of block size)
ale <- rbindlist(lapply(1:REP, function(r) { set.seed(r); f <- sample(rep_len(1:K, nrow(D)))
  rbindlist(lapply(1:K, function(k) { s <- dividir(f == k)
    data.table(rep = r, AUC = auc(s$te$pres, pred(gbm(F_ALL, s$ent), F_ALL, s$te))) })) }))
AUC_ALE <- ale[, mean(AUC), by = rep][, mean(V1)]

porrep <- R[, .(bloques = first(bloques), base = mean(base_AUC), obs = mean(obs_AUC), neu = mean(neu_AUC),
                TSS = mean(obs_TSS), PR = mean(obs_PR), cambio = mean(cambio_decil)), by = .(lado_km, rep)]
res <- porrep[, .(bloques = round(mean(bloques)), AUC_base = mean(base), AUC_gbm = mean(obs),
                  AUC_gbm_sd = sd(obs), AUC_neu = mean(neu), dif_gbm_base = mean(obs - base),
                  dif_neu_obs = mean(neu - obs), optimismo_aleatorio = AUC_ALE - mean(obs),
                  TSS = mean(TSS), PR_AUC = mean(PR), cambio_decil = mean(cambio)), by = lado_km]
cat("\n===== BLOCK SIZE SENSITIVITY (", REP, "reps x", K, "folds; random-fold AUC =", round(AUC_ALE,3), ") =====\n")
print(res[, lapply(.SD, function(x) round(x, 3))])
fwrite(R,   file.path(RUTA, "salidas/tablas/tabla_tamano_bloque_reps.csv"))
fwrite(res, file.path(RUTA, "salidas/tablas/tabla_tamano_bloque.csv"))
