# Prospectividad paleontológica del Perú

Inventario nacional de localidades fósiles y modelado predictivo de prospectividad
con control explícito del sesgo de muestreo, implementado íntegramente en R.

**Estado:** análisis completo (versión v1.1.0).
**Fecha de corte del inventario:** 18 de agosto de 2026.

## Resultado en una línea

Sobre 323 084 celdas de 2 km, 655 localidades independientes (699 celdas de presencia)
y covariables de preservación, exposición (Sentinel-2) y observación, el *gradient
boosting* alcanza un AUC de 0,866 bajo bloqueo espacial (10 repeticiones) frente a 0,759
del modelo base geológico. Neutralizar el sesgo de accesibilidad en la predicción baja el
AUC validado en 0,052 pero cambia el 37,9 % del decil prioritario: la métrica penaliza el
mapa corregido. Se mantiene con bloques de 50 a 200 km.

## Cómo reproducir

Requisitos: R ≥ 4.6 con `renv`, Python 3, GDAL (`gdal_translate`) y `curl`.

```r
renv::restore()      # 70 paquetes fijados en renv.lock
```

```sh
sh reproducir.sh     # toda la cadena, en orden; o  sh reproducir.sh 07  para retomar
```

Los datos **no se versionan**: los scripts los descargan de sus fuentes públicas.
Cada paso aleatorio fija su semilla dentro del script (tabla S1 del manuscrito);
los hiperparámetros están en `salidas/tablas/tabla_hiperparametros.csv`.
El orden y la función de cada script están comentados en `reproducir.sh`.

## Fuentes de datos

| Fuente | Uso | Acceso |
|---|---|---|
| Paleobiology Database | Ocurrencias fósiles | API pública |
| GBIF | Especímenes fósiles con coordenadas | API pública |
| INGEMMET / GEOCATMIN | Cartas geológicas 1:50 000 y 1:100 000 | [API REST](https://geocatmin.ingemmet.gob.pe/arcgis/rest/services) |
| SRTM (vía `geodata`) | Relieve | Descarga automática |
| OpenStreetMap (vía `geodata`) | Vías y centros poblados | Descarga automática |
| HydroRIVERS | Cursos de agua | Descarga automática |
| GADM | Límites administrativos | Descarga automática |

La colección paleontológica del INGEMMET (~20 000 especímenes) **no** está integrada:
requiere solicitud institucional.

## Notas técnicas

Cuatro trampas encontradas durante el desarrollo, documentadas aquí porque
producen resultados que parecen correctos:

1. `maxAllowableOffset` de ArcGIS se interpreta en las unidades de `outSR`. Con
   `outSR=4326` son **grados**: usar `0.001`, no `100`.
2. `terra` 1.9.34: `makeValid()` colapsa la capa completa (157 105 polígonos a 18)
   dejando los atributos descolgados. No usarlo aquí.
3. Al escribir GeoTIFF multicapa con capas categóricas, las celdas sin dato quedan
   con el centinela `-2147483648` y R las lee como una categoría válida, lo que
   infla la cobertura. Se guarda en `.grd` nativo.
4. El campo `AMBIENTE_S` de GEOCATMIN usa el dominio `DGR_AMBIETE_SED`:
   `1` = Continental, `2` = De transición, `3` = Marino.

## Estructura

```
R/            scripts de análisis, en orden de ejecución
datos/        crudos y procesados (no versionados; se regeneran)
salidas/      figuras y tablas del artículo

```

## Limitaciones conocidas

- La carta 1:50 000 cubre el 61 % del país; el resto se completa con la 1:100 000 e imputación.
- El PR-AUC absoluto es bajo (0,054) y la salida no está calibrada: es un puntaje relativo, no una probabilidad.
- La transferencia a regiones sin historia de colecta es débil (AUC 0,59-0,68 al retirar una región entera).
- El límite entre las cartas 1:50 000 y 1:100 000 deja artefactos rectos (p. ej. en Madre de Dios).
- Los `id` de ocurrencias se escribieron sin `bit64`: no sirven para rastrear el registro original.

## Licencia

El código se distribuye con licencia [MIT](LICENSE).
Los datos derivados (tablas y mapas en `salidas/` y `datos/procesados/`) siguen los
términos de sus fuentes: PBDB (CC BY 4.0), GBIF (licencia de cada conjunto de datos),
INGEMMET (obras derivadas con atribución), ESA WorldCover (CC BY 4.0) y
OpenStreetMap (ODbL). Los límites de GADM no se redistribuyen.

## Cita

[Pendiente: referencia del artículo cuando esté publicado.]
