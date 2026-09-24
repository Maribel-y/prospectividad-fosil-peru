# ============================================================
# Where the model fails: performance of each spatial fold (script 14,
# national gradient boosting, mean of the repetitions) together with
# the departments that contribute most presences to each fold.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")
R <- fread(file.path(RUTA, "salidas/tablas/tabla_metricas_por_pliegue.csv"))
D <- fread(file.path(RUTA, "datos/procesados/tabla_modelado.csv"), select = c("celda","x","y","pres","fold"))
dep <- project(readRDS(file.path(RUTA, "datos/crudos/gadm/gadm41_PER_1_pk.rds")),
               "+proj=laea +lat_0=-10 +lon_0=-75 +datum=WGS84 +units=m")
P <- D[pres == 1]
P[, departamento := extract(dep, vect(as.matrix(P[, .(x, y)]), crs = crs(dep)))$NAME_1]

m <- R[design == "National" & model == "Gradient boosting",
       .(AUC = mean(AUC), PR_AUC = mean(PR_AUC), TSS = mean(TSS)), by = fold]
top <- P[, .N, by = .(fold, departamento)][order(fold, -N)][, .(presences = sum(N),
         main_departments = paste(head(sprintf("%s (%d)", departamento, N), 3), collapse = ", ")), by = fold]
tab <- merge(m, top, by = "fold")[order(AUC)]
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(tab, file.path(RUTA, "salidas/tablas/tabla_rendimiento_por_pliegue.csv"))

# ---- Summary of the leave-one-department-out block of script 20 ----
LD <- fread(file.path(RUTA, "salidas/tablas/tabla_validacion_departamento.csv"))
res <- data.table(
  indicador = c("departments", "median AUC baseline", "median AUC observational", "median AUC neutralised",
                "departments neutralised > baseline", "rho density vs neutralised AUC",
                "rho density vs neutralised PR-AUC"),
  valor = c(nrow(LD), median(LD$base_AUC), median(LD$all_obs_AUC), median(LD$all_neu_AUC),
            sum(LD$all_neu_AUC > LD$base_AUC),
            cor(LD$pres_per_1000km2, LD$all_neu_AUC, method = "spearman"),
            cor(LD$pres_per_1000km2, LD$all_neu_PR, method = "spearman")))
print(res[, .(indicador, valor = round(valor, 3))])
fwrite(res, file.path(RUTA, "salidas/tablas/tabla_validacion_departamento_resumen.csv"))
