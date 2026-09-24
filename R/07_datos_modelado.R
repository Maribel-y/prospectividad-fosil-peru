# ============================================================
# Fase 4 (parte A) — Tabla de modelado y control de calidad
# Une geología, topografía y accesibilidad por celda, y marca las
# presencias geológicamente imposibles detectadas en la Fase 3.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
RUTA <- path.expand("~/Desktop/dinos")

D    <- fread(file.path(RUTA,"datos/procesados/tabla_analisis.csv"))
topo <- rast(file.path(RUTA,"datos/procesados/covar_topo_2km.tif"))
acc  <- rast(file.path(RUTA,"datos/procesados/covar_acceso_2km.tif"))
expo <- rast(file.path(RUTA,"datos/procesados/covar_exposicion_2km.tif"))   # guion 02c

for (n in names(topo)) D[[n]] <- values(topo[[n]])[D$celda]
for (n in names(acc))  D[[n]] <- values(acc[[n]])[D$celda]
D[["exposicion"]] <- values(expo)[D$celda]

# ---- Control de calidad: presencias en roca intrusiva ---------
# Un batolito o una tonalita no puede contener fósiles: esas presencias
# delatan error de georreferencia o de asignación de unidad.
IGNEA <- paste0("batolito|tonalita|granodiorita|granito|diorita|monzonita|",
                "gabro|sienita|pórfido|porfido|ortogneis|gneis|migmatita|",
                "intrusivo|plutón|pluton|dique|stock")
D[, unidad_ignea := grepl(IGNEA, tolower(uni))]
sosp <- D[pres == 1 & unidad_ignea]
cat("===== CONTROL DE CALIDAD =====\n")
cat("Presencias en unidades intrusivas:", nrow(sosp), "de", sum(D$pres),
    sprintf("(%.1f %%)\n", 100*nrow(sosp)/sum(D$pres)))
print(head(sosp[, .N, by = uni][order(-N)], 8))

# ---- Preparación de covariables -------------------------------
D[, era := factor(fifelse(is.na(era) | era == "", "Sin dato", era))]
D[, amb := factor(fifelse(is.na(amb) | amb == "", "Sin dato", amb))]
D[, edad_na := as.integer(is.na(edad))]
D[, edad_imp := fifelse(is.na(edad), median(edad, na.rm = TRUE), edad)]
D[, log_edad := log10(edad_imp + 1)]

num <- c("elevacion","pendiente","rugosidad","exposicion","dist_vias_km","dist_poblados_km","dist_rios_km")
for (n in num) D[[n]][is.na(D[[n]])] <- median(D[[n]], na.rm = TRUE)

cat("\nCovariables por celda:", paste(c("era","amb","log_edad","edad_na", num), collapse=", "), "\n")
cat("Celdas completas:", sum(complete.cases(D[, c("era","amb","log_edad", num), with = FALSE])),
    "de", nrow(D), "\n")

fwrite(D, file.path(RUTA,"datos/procesados/tabla_modelado.csv"))
fwrite(sosp[, .(celda, x, y, uni, era, amb)],
       file.path(RUTA,"salidas/tablas/tabla_presencias_sospechosas.csv"))
cat("\nTabla de modelado escrita:", nrow(D), "celdas x", ncol(D), "columnas\n")
