#!/bin/sh
# ============================================================
# Descarga el percentil 90 anual de NDVI de Sentinel-2 (2021) que
# publica ESA WorldCover (compuestos de 10 m, mosaicos de 1 grado;
# datos abiertos CC-BY 4.0, sin cuenta). Solo se lee la banda 1
# (NDVI p90) en su primer resumen interno, de ~20 m (6000 x 6000 px
# por grado): suficiente para resolver afloramientos de decenas de
# metros sin descargar los 10 m nativos.
#
# Cada mosaico se reduce de inmediato a conteos por celda de 2 km
# (guion 02c_conteos_mosaico.R) y el ráster de 20 m se borra, de modo
# que en disco solo quedan los conteos (datos/crudos/sentinel2_conteos).
# Si existe $OUT/necesarios.txt (guion 02a) solo se procesan esos mosaicos.
# Uso: sh R/02b_descargar_sentinel2.sh   (desde la raíz del proyecto)
# ============================================================
export OUT=datos/crudos/sentinel2_conteos
export TMPD="${TMPDIR:-/tmp}/s2_ndvi20"
export BASE=/vsicurl/https://esa-worldcover-s2.s3.eu-central-1.amazonaws.com/ndvi/2021
mkdir -p "$OUT" "$TMPD"

# Mosaicos de 1 grado dentro de los 26 mosaicos de 3 grados de WorldCover
for t in $(cat datos/crudos/worldcover/mosaicos.txt); do
  lat=$(echo "$t" | cut -c2-3); lon=$(echo "$t" | cut -c5-7)
  for dy in 0 1 2; do for dx in 0 1 2; do
    printf "S%02dW%03d\n" $((10#$lat - dy)) $((10#$lon - dx))
  done; done
done | sort -u | { if [ -s "$OUT/necesarios.txt" ]; then grep -F -x -f "$OUT/necesarios.txt"; else cat; fi; } | while read -r m; do
  [ -s "$OUT/conteos_${m}.tif" ] || [ -e "$OUT/sin_mosaico_${m}" ] || echo "$m"
done | xargs -P 6 -n 1 sh -c '
  m=$1; lat=$(echo "$m" | cut -c1-3); f="$TMPD/ndvi20_${m}.tif"
  if gdal_translate -q -b 1 -ovr 0 -co COMPRESS=DEFLATE \
       "$BASE/$lat/ESA_WorldCover_10m_2021_v200_${m}_NDVI.tif" "$f" 2>/dev/null; then
    Rscript R/02c_conteos_mosaico.R "$f" "$OUT/conteos_${m}.tif" >/dev/null 2>&1 \
      && echo "ok $m" || echo "FALLO conteo $m"
    rm -f "$f"
  else
    touch "$OUT/sin_mosaico_${m}"; echo "sin mosaico (mar): $m"
  fi' _
ls "$OUT"/conteos_*.tif | wc -l
