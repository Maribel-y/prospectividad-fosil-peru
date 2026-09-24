# ============================================================
# Systematic ablation: what each component of the framework adds,
# on identical spatial folds and background, repeated with REP seeds.
#   base      geological baseline (WoE + logistic GLM)
#   geo       GBM, geology only
#   geotopo0  GBM, geology + topography
#   geotopo   GBM, geology + topography + rock exposure
#   all_obs   GBM, geology + topography + access, observational projection
#   all_neu   same model, access neutralized at prediction (proposed)
#   all_rand  all_obs evaluated with random folds (no spatial blocking)
# desp_* = % of the top decile of each variant that is not in the top
# decile of the proposed map (priority-set displacement).
# ============================================================
source(file.path(path.expand("~/Desktop/dinos"), "R/comun_modelado.R"))
K <- 5; REP <- 10
filas <- list()
for (r in 1:REP) {
  set.seed(r); D[, fold_aleatorio := sample(rep_len(1:K, .N))]
  for (k in 1:K) {
    s <- dividir(D$fold == k); te <- s$te; y <- te$pres; dv <- te$dist_vias_km
    p_base <- base_glm(s$tr, te)
    p_geo  <- pred(gbm(F_GEO, s$ent), F_GEO, te)
    p_gt0  <- pred(gbm(F_GT0, s$ent), F_GT0, te)
    p_gt   <- pred(gbm(F_GT,  s$ent), F_GT,  te)
    m_all  <- gbm(F_ALL, s$ent)
    p_obs  <- pred(m_all, F_ALL, te)
    p_neu  <- pred(m_all, F_ALL, neutralizar(te, s$tr))
    sa <- dividir(D$fold_aleatorio == k)
    p_rand <- pred(gbm(F_ALL, sa$ent), F_ALL, sa$te)
    filas[[length(filas)+1]] <- as.data.table(c(list(rep=r, fold=k),
      metricas(y,p_base,dv,"base"), metricas(y,p_geo,dv,"geo"), metricas(y,p_gt0,dv,"geotopo0"), metricas(y,p_gt,dv,"geotopo"),
      metricas(y,p_obs,dv,"all_obs"), metricas(y,p_neu,dv,"all_neu"),
      metricas(sa$te$pres,p_rand,sa$te$dist_vias_km,"all_rand"),
      list(desp_base=desplazamiento(p_base,p_neu), desp_geo=desplazamiento(p_geo,p_neu),
           desp_geotopo0=desplazamiento(p_gt0,p_neu), desp_geotopo=desplazamiento(p_gt,p_neu), desp_all_obs=desplazamiento(p_obs,p_neu))))
  }
  cat("repetition", r, "of", REP, "done\n")
}
R <- rbindlist(filas)
fwrite(R, file.path(RUTA,"salidas/tablas/tabla_ablacion_reps.csv"))
porrep <- R[, lapply(.SD, mean, na.rm=TRUE), by=rep, .SDcols=-c("rep","fold")]
vars <- setdiff(names(porrep), "rep")
res <- data.table(indicador=vars, media=sapply(vars, function(v) mean(porrep[[v]])),
                  sd_entre_rep=sapply(vars, function(v) sd(porrep[[v]])))
cat("\n===== ABLATION (", REP, "reps x", K, "spatial folds ) =====\n")
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x,4) else x)], nrows=100)
fwrite(res, file.path(RUTA,"salidas/tablas/tabla_ablacion.csv"))
