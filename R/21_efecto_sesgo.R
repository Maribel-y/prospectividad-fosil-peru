# ============================================================
# Effect of bias neutralization on the national map (Table 7):
# Spearman correlation of each projection with distance to roads,
# settlements and rivers; share of the top decile within 5 km of a
# river; top-decile overlap between projections; and overlap between
# the national and Mesozoic priority sets.
# Optional arguments: processed-data dir and tables dir (defaults: project).
# ============================================================
suppressPackageStartupMessages(library(data.table))
RUTA <- path.expand("~/Desktop/dinos")
a <- commandArgs(TRUE)
PROC <- if (length(a) >= 1) a[1] else file.path(RUTA, "datos/procesados")
TAB  <- if (length(a) >= 2) a[2] else file.path(RUTA, "salidas/tablas")
P <- fread(file.path(PROC, "predicciones.csv"))
M <- fread(file.path(PROC, "tabla_modelado.csv"), select = c("celda","dist_vias_km","dist_poblados_km","dist_rios_km"))
P <- merge(P, M, by = "celda")
rho <- function(p, v) cor(p, v, method = "spearman")
top <- function(p) p >= quantile(p, 0.9)
t_obs <- top(P$prosp_nacional); t_neu <- top(P$prospcorr_nacional)
f <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",")
prio_n <- fread(file.path(TAB, "tabla_priorizacion_nacional.csv"))
prio_m <- fread(file.path(TAB, "tabla_priorizacion_mesozoico.csv"))
tab <- data.table(
  metrica = c("Correlación con distancia a vías (Spearman)", "Correlación con distancia a poblados",
              "Correlación con distancia a ríos", "Decil superior a menos de 5 km de un río (%)",
              "Solapamiento del decil superior entre proyecciones (%)",
              "Solapamiento de celdas priorizadas nacional vs mesozoico"),
  observacional = c(f(rho(P$prosp_nacional, P$dist_vias_km)), f(rho(P$prosp_nacional, P$dist_poblados_km)),
                    f(rho(P$prosp_nacional, P$dist_rios_km)), f(100*mean(P$dist_rios_km[t_obs] < 5), 1), "", ""),
  neutralizado  = c(f(rho(P$prospcorr_nacional, P$dist_vias_km)), f(rho(P$prospcorr_nacional, P$dist_poblados_km)),
                    f(rho(P$prospcorr_nacional, P$dist_rios_km)), f(100*mean(P$dist_rios_km[t_neu] < 5), 1),
                    f(100*sum(t_obs & t_neu)/sum(t_obs), 1),
                    sprintf("%d de %d y %d", length(intersect(prio_n$celda, prio_m$celda)), nrow(prio_n), nrow(prio_m))))
print(tab)
cat("Territory within 5 km of a river (%):", f(100*mean(P$dist_rios_km < 5), 1), "\n")
fwrite(tab, file.path(TAB, "tabla_efecto_sesgo.csv"))
