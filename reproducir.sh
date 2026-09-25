#!/bin/sh
# ============================================================
# Reproduces every number, table and figure of the manuscript.
# Usage, from the project root:   sh reproducir.sh [desde]
#   desde  optional step number to resume from (e.g. 07)
# Requirements: R >= 4.6 with renv::restore() done, Python 3, GDAL
# (gdal_translate) and curl. Source data are downloaded by the scripts.
# Every stochastic step fixes its seed inside the script (see the
# seeds table in the supplementary material). Logs go to salidas/*.log.
# ============================================================
set -e
DESDE=${1:-00}
paso() { n=$1; shift
  if [ "$(printf '%s' "$n" | cut -c1-2)" \< "$DESDE" ]; then return; fi
  echo ">> $n"; "$@"; }
R() { Rscript "R/$1" > "salidas/log_$(basename "$1" .R).log" 2>&1; }

paso 00  R 00_censo_viabilidad.R          # PBDB + GBIF, census
paso 00b R 00b_revisar_descartados.R
paso 00c R 00c_recuperar_costeras.R       # coastal recovery, degree-precision filter
paso 01  R 01_grilla_presencias.R         # 2 km grid and response
paso 02  R 02_covariables_topo.R          # SRTM relief
paso 02a R 02a_mosaicos_necesarios.R     # Sentinel-2 tiles that touch Peru
paso 02b sh R/02b_descargar_sentinel2.sh  # Sentinel-2 NDVI p90 -> counts per cell
paso 02c R 02c_covar_exposicion.R         # rock exposure covariate
paso 03  python3 R/descargar_geologia.py
paso 03  python3 R/descargar_geologia_100k.py
paso 03  R 03_covar_geologia.R
paso 04  R 04_covar_accesibilidad.R       # OSM roads/settlements, HydroRIVERS
paso 05  R 05_geologia_combinada.R        # 1:50 000 + 1:100 000 with imputation
paso 05b R 05b_limpiar_geologia.R
paso 06  R 06_modelo_base_geologico.R     # WoE baseline, 100 km blocks (seed 42)
paso 07  R 07_datos_modelado.R            # modelling table
paso 08  R 08_entrenamiento.R             # single-seed comparison (seed 42)
paso 09  R 09_hipotesis.R
paso 10  R 10_prediccion_final.R          # maps, uncertainty, single-projection sets
paso 13  R 13_figuras_diagnostico.R       # ROC by fold, importance
paso 14  R 14_metricas_revisor.R          # Tables 4-5: 10 seeds x 5 folds, both designs
paso 15  R 15_sensibilidad.R              # neutralisation constant, WoE stability
paso 16  R 16_priorizacion_robusta.R      # robust priority set
paso 18  R 18_metrica_vs_mapa.R           # Table 6: metric versus map (seeds 1-10)
paso 19  R 19_ablacion.R                  # ablation (seeds 1-10)
paso 20  R 20_validacion_regional.R       # leave-one-region / department-out
paso 21  R 21_efecto_sesgo.R              # Table 7
paso 22  R 22_localidades_vs_celdas.R     # localities vs presence cells
paso 23  R 23_tamano_bloque.R             # block size 50-200 km (seeds 1-5)
paso 24  R 24_calibracion.R               # Brier score, reliability (seeds 1-5)
paso 25  R 25_dependencia_imputacion.R    # priority sets vs imputed geology
paso 26  R 26_pliegue_debil.R             # weakest fold, department summary
paso 27  R 27_gbif_derived_dataset.R      # datasetKey counts for the GBIF derived dataset
paso 28  R 28_figura_area_estudio.R      # study-area figure
paso 17  R 17b_figuras_compuestas_en.R    # composite figures (English)
echo "Done. Tables in salidas/tablas, figures in salidas/figuras."
