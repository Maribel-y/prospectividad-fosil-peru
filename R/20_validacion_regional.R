# ============================================================
# Geographic transferability: leave-one-region-out validation.
# Natural regions from elevation and the Andean divide:
#   Sierra = elevation >= 1000 m; Costa / Selva = below 1000 m, west / east
#   of the divide, taken per 20-km latitude band as the median x of the
#   10 % highest cells. Train on two regions, test on the third, plus the
#   hardest case: train on the densely sampled coast, test on the rest.
# Second block: leave-one-department-out for every department with at
# least MIN_DEP presence cells, with a 10 km exclusion buffer around the
# held-out department so that neighbouring cells do not leak into
# training. Relates transferability to how densely each department has
# been sampled (presences per 1 000 km2).
# ============================================================
source(file.path(path.expand("~/Desktop/dinos"), "R/comun_modelado.R"))
REP <- 10; MIN_DEP <- 15; BUFFER <- 10000
D[, banda := round(y/20000)]
cresta <- D[, .(xc = median(x[elevacion >= quantile(elevacion, 0.9, na.rm=TRUE)], na.rm=TRUE)), by=banda]
D <- merge(D, cresta, by="banda", sort=FALSE)
D[, region := fifelse(elevacion >= 1000, "Sierra", fifelse(x < xc, "Costa", "Selva"))]
print(D[, .(cells=.N, presences=sum(pres), prevalence_pct=round(100*mean(pres),3)), by=region])
fwrite(D[, .(celda, region)], file.path(RUTA,"datos/procesados/region_natural.csv"))

escenarios <- list("Test Costa"  = quote(region == "Costa"),
                   "Test Sierra" = quote(region == "Sierra"),
                   "Test Selva"  = quote(region == "Selva"),
                   "Train Costa only, test Sierra + Selva" = quote(region != "Costa"))
filas <- list()
for (r in 1:REP) for (nm in names(escenarios)) {
  set.seed(r)
  es_test <- D[, eval(escenarios[[nm]])]
  s <- if (grepl("only", nm)) {
    x <- dividir(es_test); x$tr <- x$tr[region == "Costa"]
    x$ent <- rbind(x$tr[pres==1], x$tr[pres==0][sample(.N, min(.N, N_FONDO))]); x
  } else dividir(es_test)
  te <- s$te; y <- te$pres; dv <- te$dist_vias_km
  m_all <- gbm(F_ALL, s$ent)
  p_obs <- pred(m_all, F_ALL, te); p_neu <- pred(m_all, F_ALL, neutralizar(te, s$tr))
  p_gt  <- pred(gbm(F_GT, s$ent), F_GT, te)
  filas[[length(filas)+1]] <- as.data.table(c(list(rep=r, scenario=nm, n_test_pres=sum(y)),
    metricas(y, base_glm(s$tr, te), dv, "base"), metricas(y,p_gt,dv,"geotopo"),
    metricas(y,p_obs,dv,"all_obs"), metricas(y,p_neu,dv,"all_neu"),
    list(desp_obs_neu = desplazamiento(p_obs, p_neu))))
  cat("rep", r, "-", nm, "done\n")
}
R <- rbindlist(filas)
fwrite(R, file.path(RUTA,"salidas/tablas/tabla_validacion_regional_reps.csv"))
res <- R[, lapply(.SD, function(x) round(mean(x, na.rm=TRUE), 4)), by=scenario, .SDcols=-c("rep","scenario")]
cat("\n===== LEAVE-ONE-REGION-OUT (mean of", REP, "reps) =====\n")
print(t(res))
fwrite(res, file.path(RUTA,"salidas/tablas/tabla_validacion_regional.csv"))

# ---- Leave-one-department-out ---------------------------------------
suppressPackageStartupMessages(library(terra))
dep <- project(readRDS(file.path(RUTA, "datos/crudos/gadm/gadm41_PER_1_pk.rds")), PROJ_LAEA <-
               "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m")
pts <- vect(as.matrix(D[, .(x, y)]), crs = PROJ_LAEA)
D[, departamento := extract(dep, pts)$NAME_1]
dens <- D[!is.na(departamento), .(cells = .N, presences = sum(pres)), by = departamento]
dens[, pres_per_1000km2 := 1000*presences/(4*cells)]
deps <- dens[presences >= MIN_DEP][order(-pres_per_1000km2)]$departamento
cat("\nDepartments held out:", paste(deps, collapse = ", "), "\n")
filas <- list()
for (dp in deps) {
  en_dep <- D$departamento %in% dp
  zona <- buffer(dep[dep$NAME_1 == dp], BUFFER)                   # department + 10 km
  cerca <- !is.na(extract(zona, pts)[, 2])
  for (r in 1:REP) { set.seed(r)
    s <- dividir(en_dep); s$tr <- s$tr[!(celda %in% D$celda[cerca & !en_dep])]
    s$ent <- rbind(s$tr[pres==1], s$tr[pres==0][sample(.N, min(.N, N_FONDO))])
    te <- s$te; y <- te$pres; dv <- te$dist_vias_km
    m_all <- gbm(F_ALL, s$ent)
    filas[[length(filas)+1]] <- as.data.table(c(list(rep = r, departamento = dp, n_test_pres = sum(y)),
      metricas(y, base_glm(s$tr, te), dv, "base"),
      metricas(y, pred(m_all, F_ALL, te), dv, "all_obs"),
      metricas(y, pred(m_all, F_ALL, neutralizar(te, s$tr)), dv, "all_neu")))
  }
  cat("department", dp, "done\n")
}
RD <- rbindlist(filas)
resd <- RD[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = departamento, .SDcols = -c("rep","departamento")]
resd <- merge(resd, dens, by = "departamento")[order(-pres_per_1000km2)]
cat("\n===== LEAVE-ONE-DEPARTMENT-OUT (mean of", REP, "reps, buffer", BUFFER/1000, "km) =====\n")
print(resd[, .(departamento, presences, pres_per_1000km2 = round(pres_per_1000km2, 2),
               base_AUC = round(base_AUC, 3), obs_AUC = round(all_obs_AUC, 3), neu_AUC = round(all_neu_AUC, 3),
               neu_PR = round(all_neu_PR, 4))])
cat("Spearman between sampling density and neutralised AUC:",
    round(cor(resd$pres_per_1000km2, resd$all_neu_AUC, method = "spearman"), 3), "\n")
fwrite(RD,   file.path(RUTA, "salidas/tablas/tabla_validacion_departamento_reps.csv"))
fwrite(resd, file.path(RUTA, "salidas/tablas/tabla_validacion_departamento.csv"))
