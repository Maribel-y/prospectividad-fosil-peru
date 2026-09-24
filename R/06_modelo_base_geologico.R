# ============================================================
# Fase 3 — Modelo base geológico (línea base empírica obligatoria)
#
# Pesos de Evidencia (WoE) sobre las covariables geológicas y regresión
# logística evaluada con bloqueo espacial. Ningún modelo de aprendizaje
# automático posterior se acepta si no supera claramente este resultado.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
set.seed(42)
RUTA <- path.expand("~/Desktop/dinos")

g   <- rast(file.path(RUTA,"datos/procesados/grilla_2km.tif"))
geo <- rast(file.path(RUTA,"datos/procesados/covar_geologia_2km.grd"))
dic <- lapply(c("era_geol","ambiente","unidad"), function(n)
        fread(file.path(RUTA,"datos/procesados", paste0("dic_", n, ".csv"))))
names(dic) <- c("era_geol","ambiente","unidad")
pres_celdas <- fread(file.path(RUTA,"datos/procesados/celdas_presencia.csv"))$celda
occ <- fread(file.path(RUTA,"datos/procesados/ocurrencias_con_celda.csv"))
mes_celdas <- unique(occ[era == "Mesozoico", celda])

# ---- Tabla de análisis ---------------------------------------
ter <- which(!is.na(values(g$tierra)))
D <- data.table(celda = ter,
                pres  = as.integer(ter %in% pres_celdas),
                pres_mes = as.integer(ter %in% mes_celdas),
                era_cod = values(geo$era_geol)[ter],
                amb_cod = values(geo$ambiente)[ter],
                uni_cod = values(geo$unidad)[ter],
                edad    = values(geo$edad_media)[ter])
xy <- xyFromCell(g, ter); D[, `:=`(x = xy[,1], y = xy[,2])]

etq <- function(cod, d) { v <- d$etiqueta[match(cod, d$codigo)]; fifelse(is.na(v) | v %in% c("","<Null>"), "Sin dato", v) }
D[, era := etq(era_cod, dic$era_geol)]
D[, amb := etq(amb_cod, dic$ambiente)]
D[, uni := etq(uni_cod, dic$unidad)]
cat("Celdas:", nrow(D), "| presencias:", sum(D$pres), "| mesozoicas:", sum(D$pres_mes), "\n\n")

# ---- Pesos de Evidencia --------------------------------------
woe <- function(dt, var, resp = "pres") {
  N <- nrow(dt); Dn <- sum(dt[[resp]])
  dt[, .(celdas = .N, presencias = sum(get(resp))), by = c(var)][
    , `:=`(
      Wmas  = log(((presencias + 0.5)/Dn) / ((celdas - presencias + 0.5)/(N - Dn))),
      Wmenos= log(((Dn - presencias + 0.5)/Dn) / ((N - celdas - Dn + presencias + 0.5)/(N - Dn))))][
    , contraste := Wmas - Wmenos][
    , s_contraste := sqrt(1/(presencias+0.5) + 1/(celdas-presencias+0.5) +
                          1/(Dn-presencias+0.5) + 1/(N-celdas-Dn+presencias+0.5))][
    , C_estudentizado := round(contraste/s_contraste, 2)][order(-contraste)]
}

cat("===== WoE POR ERA GEOLÓGICA =====\n")
w_era <- woe(D, "era"); print(w_era[, .(era, celdas, presencias, contraste = round(contraste,2), C_estudentizado)])
cat("\n===== WoE POR AMBIENTE SEDIMENTARIO =====\n")
w_amb <- woe(D, "amb"); print(w_amb[, .(amb, celdas, presencias, contraste = round(contraste,2), C_estudentizado)])
cat("\n===== WoE: 12 UNIDADES MÁS FAVORABLES (>= 5 presencias) =====\n")
w_uni <- woe(D, "uni")[presencias >= 5]
print(head(w_uni[, .(uni, celdas, presencias, contraste = round(contraste,2), C_estudentizado)], 12))

fwrite(w_era, file.path(RUTA,"salidas/tablas/woe_era.csv"))
fwrite(w_amb, file.path(RUTA,"salidas/tablas/woe_ambiente.csv"))
fwrite(w_uni, file.path(RUTA,"salidas/tablas/woe_unidad.csv"))

# ---- Bloques espaciales --------------------------------------
# Bloques de 100 km asignados a 5 pliegues: evita que celdas vecinas
# queden a la vez en entrenamiento y evaluación.
LADO <- 100000; K <- 5
D[, bloque := paste0(floor(x/LADO), "_", floor(y/LADO))]
bl <- unique(D$bloque)
D[, fold := match(bloque, bl)]
D[, fold := sample(rep_len(1:K, length(bl)))[fold]]
cat("\nBloques espaciales de", LADO/1000, "km:", length(bl), "| pliegues:", K, "\n")
print(D[, .(celdas = .N, presencias = sum(pres)), by = fold][order(fold)])

# ---- Métricas ------------------------------------------------
auc <- function(y, p) { r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0)
                        (sum(r[y==1]) - n1*(n1+1)/2) / (n1*n0) }
tss <- function(y, p) { u <- unique(quantile(p, seq(0.5,0.999,0.001), na.rm = TRUE))
  max(sapply(u, function(t) { s <- mean(p[y==1] >= t); e <- mean(p[y==0] < t); s + e - 1 })) }
boyce <- function(y, p, nclas = 10) {   # índice de Boyce continuo
  # Un modelo puramente categórico produce pocas probabilidades distintas,
  # así que los cortes se toman sobre los valores únicos disponibles.
  br <- unique(quantile(p, seq(0, 1, length.out = nclas + 1), na.rm = TRUE))
  if (length(br) < 3) return(NA_real_)
  br[1] <- -Inf; br[length(br)] <- Inf
  cl <- cut(p, br, labels = FALSE)
  k  <- length(br) - 1
  pe <- sapply(1:k, function(i) { E <- sum(cl == i)/length(cl)
                                  if (E == 0) NA else (sum(y[cl == i] == 1)/sum(y == 1))/E })
  if (sum(!is.na(pe)) < 3) return(NA_real_)
  suppressWarnings(cor(1:k, pe, method = "spearman", use = "complete.obs"))
}

# WoE de unidad calculado SOLO con el pliegue de entrenamiento (sin fuga)
woe_unidad_train <- function(train, test) {
  w <- woe(train, "uni")[, .(uni, Wmas)]
  v <- w$Wmas[match(test$uni, w$uni)]
  fifelse(is.na(v), 0, v)
}

evaluar <- function(formula, usar_woe = FALSE, resp = "pres") {
  res <- rbindlist(lapply(1:K, function(k) {
    tr <- D[fold != k]; te <- D[fold == k]
    if (sum(te[[resp]]) < 3 || sum(tr[[resp]]) < 10) return(NULL)
    if (usar_woe) { tr <- copy(tr); te <- copy(te)
      tr[, woe_uni := woe_unidad_train(tr, tr)]; te[, woe_uni := woe_unidad_train(D[fold != k], te)] }
    m <- suppressWarnings(glm(formula, data = tr, family = binomial))
    p <- suppressWarnings(predict(m, newdata = te, type = "response"))
    ok <- !is.na(p)
    data.table(fold = k, AUC = auc(te[[resp]][ok], p[ok]),
               TSS = tss(te[[resp]][ok], p[ok]), Boyce = boyce(te[[resp]][ok], p[ok]))
  }))
  res[, .(AUC = mean(AUC), AUC_sd = sd(AUC), TSS = mean(TSS), Boyce = mean(Boyce, na.rm = TRUE))]
}

cat("\n===== MODELO BASE GEOLÓGICO — validación con bloqueo espacial =====\n")
modelos <- list(
  "M0: solo era"                    = list(pres ~ era, FALSE),
  "M1: era + ambiente"              = list(pres ~ era + amb, FALSE),
  "M2: era + ambiente + WoE unidad" = list(pres ~ era + amb + woe_uni, TRUE))
tab <- rbindlist(lapply(names(modelos), function(n) {
  r <- evaluar(modelos[[n]][[1]], modelos[[n]][[2]]); r[, modelo := n]; r }))
setcolorder(tab, "modelo")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])
fwrite(tab, file.path(RUTA,"salidas/tablas/tabla_modelo_base.csv"))

# ---- Modelo final y mapa de favorabilidad base ---------------
D[, woe_uni := woe_unidad_train(D, D)]
mfinal <- suppressWarnings(glm(pres ~ era + amb + woe_uni, data = D, family = binomial))
D[, favorabilidad := predict(mfinal, type = "response")]
r_fav <- rast(g$tierra); values(r_fav) <- NA
r_fav[D$celda] <- D$favorabilidad
names(r_fav) <- "favorabilidad_base"
writeRaster(r_fav, file.path(RUTA,"datos/procesados/favorabilidad_base.tif"), overwrite = TRUE)
fwrite(D[, .(celda, x, y, pres, pres_mes, era, amb, uni, edad, fold, favorabilidad)],
       file.path(RUTA,"datos/procesados/tabla_analisis.csv"))
cat("\nModelo base ajustado y mapa de favorabilidad escrito.\n")
