"""Descarga la capa Litología del Mapa Geológico Integrado 1:50 000 (INGEMMET, ago-2025)
vía la API REST de GEOCATMIN. El servicio no soporta paginación: se avanza por rangos de
OBJECTID. Geometría simplificada a 0.001 grados (~110 m; la grilla de análisis es de 2 km).
OJO: maxAllowableOffset va en las unidades de outSR, que aquí son GRADOS, no metros.
Usa curl porque el Python del sistema no tiene bundle de certificados."""
import json, os, subprocess, time

B = ("https://geocatmin.ingemmet.gob.pe/arcgis/rest/services/"
     "SERV_GEOLOGIA_50K_INTEGRADA/MapServer/7/query")
DEST = "/Users/maribelbel/Desktop/dinos/datos/crudos/geologia_50k"
CAMPOS = ("OBJECTID,UNIDAD,LITOLOGIA,E_MAX_MA,E_MIN_MA,AMBIENTE_S,MEDIO_SEDI,"
          "MODO_EMPLA,SISTEMA_MAX,SISTEMA_MIN,SERIE_MAX,TIPO_UNIDAD,CTG_UNIDAD,ETIQUETA")
PASO, MAXID = 1000, 158000

def pedir(desde, hasta, destino, intentos=4):
    cmd = ["curl", "-s", "--max-time", "180", "-G", B,
           "--data-urlencode", f"where=OBJECTID BETWEEN {desde} AND {hasta}",
           "--data-urlencode", f"outFields={CAMPOS}",
           "--data-urlencode", "returnGeometry=true",
           "--data-urlencode", "maxAllowableOffset=0.001",
           "--data-urlencode", "outSR=4326",
           "--data-urlencode", "f=geojson", "-o", destino]
    for i in range(intentos):
        subprocess.run(cmd, check=False)
        try:
            d = json.load(open(destino))
            if "features" in d:
                return len(d["features"])
        except Exception:
            pass
        time.sleep(4 * (i + 1))
    print(f"  FALLO {desde}-{hasta}", flush=True)
    if os.path.exists(destino):
        os.remove(destino)
    return 0

total = 0
for desde in range(1, MAXID, PASO):
    hasta = desde + PASO - 1
    f_out = os.path.join(DEST, f"lote_{desde:06d}.geojson")
    if os.path.exists(f_out):
        try:
            total += len(json.load(open(f_out))["features"]); continue
        except Exception:
            os.remove(f_out)
    n = pedir(desde, hasta, f_out)
    total += n
    if n == 0 and os.path.exists(f_out):
        os.remove(f_out)          # rango vacío: no deja archivo
    if desde % 10000 == 1:
        print(f"  OBJECTID {desde:>6} … acumulado {total:>7} polígonos", flush=True)
    time.sleep(0.3)
print(f"LISTO: {total} polígonos", flush=True)
