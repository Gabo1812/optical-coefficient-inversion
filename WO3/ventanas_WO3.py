import numpy as np
import os
import shutil
import subprocess
import re
import time
import sys
from scipy.interpolate import PchipInterpolator

sys.stdout.reconfigure(encoding="utf-8")
start_time = time.time()

# =========================================================
# CONFIGURACION
# =========================================================
MUESTRA = "S5"
ESPESOR = 131.0

BASE = r"C:\Users\Gabo\Computacional\LabAvanzadoII\WO3"
archivo_exp  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_Dy_T")
archivo_nk   = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_DATANK0_DY")
film0_source = os.path.join(BASE, "FILM0.f")

output_dir  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_ventanas")
tabla_final = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_tabla_final_nk.txt")

# Ventanas: ancho variable, MUCHO solapamiento (paso << ventana)
VENTANA_TRANSICION = 25   # medio-ancho en 350-500 nm
VENTANA_PLANA      = 30   # medio-ancho en 500-1500 nm
PASO_TRANSICION    = 10   # menos denso que antes, ahorra tiempo (cleanup maneja el salto)
PASO_PLANO         = 15   # buen solapamiento en zona plana
LIMITE_ZONA        = 500
centro_inicial     = 550

def ventana_para_centro(centro):
    return VENTANA_TRANSICION if centro < LIMITE_ZONA else VENTANA_PLANA

def paso_para_centro(centro):
    return PASO_TRANSICION if centro < LIMITE_ZONA else PASO_PLANO

# =========================================================
# CARGAR DATOS
# =========================================================
data_exp = np.loadtxt(archivo_exp)
lam_exp  = data_exp[:, 0]
T_exp    = data_exp[:, 1]

data_nk = np.loadtxt(archivo_nk)
lam_nk  = data_nk[:, 0]
n0_all  = data_nk[:, 1]
k0_all  = data_nk[:, 2]

if not np.allclose(lam_exp, lam_nk):
    raise ValueError("lambda de Dy_T y DATANK0_DY no coinciden")

lam_min, lam_max = lam_exp.min(), lam_exp.max()

# Centros con paso variable
def generar_centros(lam_min, lam_max, centro_inicial):
    centros = set()
    c = centro_inicial
    while c <= lam_max:
        centros.add(c)
        c += paso_para_centro(c)
    c = centro_inicial
    while c >= lam_min:
        centros.add(c)
        c -= paso_para_centro(c)
    return sorted(centros)

centros_todos = generar_centros(lam_min, lam_max, centro_inicial)

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# =========================================================
# CREAR VENTANAS
# =========================================================
print("\nCREANDO VENTANAS (con alto solapamiento)")

ventanas         = []
nobs_por_ventana = {}
centro_por_ventana = {}
medio_ancho_por_ventana = {}

for i, centro in enumerate(centros_todos, 1):
    folder = os.path.join(output_dir, f"win_{i:03d}_c{centro}")
    os.makedirs(folder, exist_ok=True)

    medio_ancho = ventana_para_centro(centro)
    mask  = np.abs(lam_exp - centro) <= medio_ancho
    lam_w = lam_exp[mask]
    T_w   = T_exp[mask]

    if len(lam_w) == 0:
        continue

    n_w = PchipInterpolator(lam_nk, n0_all)(lam_w)
    k_w = PchipInterpolator(lam_nk, k0_all)(lam_w)
    k_w = np.clip(k_w, 0.0, None)

    with open(os.path.join(folder, "Dy_T"), "w") as f:
        for l, t in zip(lam_w, T_w):
            f.write(f"{l:.6f} {t:.6f}\n")
    with open(os.path.join(folder, "DATANK0_DY"), "w") as f:
        for l, n, k in zip(lam_w, n_w, k_w):
            f.write(f"{l:.6f} {n:.6f} {k:.6f}\n")

    nobs = len(lam_w)
    ventanas.append(folder)
    nobs_por_ventana[folder] = nobs
    centro_por_ventana[folder] = centro
    medio_ancho_por_ventana[folder] = medio_ancho

print(f"  Total de ventanas: {len(ventanas)}")

# =========================================================
# COMPILAR POR NOBS UNICO
# =========================================================
nobs_unicos = sorted(set(nobs_por_ventana.values()))
print("\nNOBS unicos:", nobs_unicos)

def patch_film0(source, dest, nobs, thickness):
    with open(source) as f:
        content = f.read()
    content = re.sub(
        r"(PARAMETER\s*\(\s*NOBS\s*=\s*)\d+(\s*,\s*NOBS2\s*=\s*2\s*\*\s*NOBS\s*\))",
        rf"\g<1>{nobs}\g<2>", content, flags=re.IGNORECASE)
    content = re.sub(
        r"(HF\s*=\s*)[0-9]+\.?[0-9]*D0",
        rf"\g<1>{thickness:.1f}D0", content)
    content = re.sub(
        r"CALL\s+SODA_LIME_GLASS\s*\(",
        "CALL FUSED_QUARTZ(", content, flags=re.IGNORECASE)

    # Limites POR PUNTO anclados a la semilla de PUMA (no un MINK/MAXK
    # uniforme por ventana). Cada punto queda confinado a un entorno
    # estrecho alrededor de SU PROPIO valor semilla:
    #   n: [semilla-0.03, semilla+0.03]
    #   k: [0.5*k_semilla, 2.0*k_semilla]  si k_semilla > 1e-8
    #      [0, 1e-6]                       si k_semilla ~ 0 (fija k≈0)
    pattern = (
        r"DO\s+I=1,NOBS\s*\r?\n"
        r"\s*L\(I\)=MINN\s*\r?\n"
        r"\s*U\(I\)=MAXN\s*\r?\n"
        r"\s*ENDDO\s*\r?\n"
        r"\s*DO\s+I=NOBS\+1,NOBS2\s*\r?\n"
        r"\s*L\(I\)=MINK\s*\r?\n"
        r"\s*U\(I\)=MAXK\s*\r?\n"
        r"\s*ENDDO"
    )
    replacement = (
        "       DO I=1,NOBS\r\n"
        "          L(I)=NF(I)-0.03D0\r\n"
        "          U(I)=NF(I)+0.03D0\r\n"
        "       ENDDO\r\n"
        "       DO I=NOBS+1,NOBS2\r\n"
        "          J=I-NOBS\r\n"
        "          IF(KF(J).GT.1.0D-8) THEN\r\n"
        "             L(I)=0.5D0*KF(J)\r\n"
        "             U(I)=2.0D0*KF(J)\r\n"
        "          ELSE\r\n"
        "             L(I)=0.0D0\r\n"
        "             U(I)=1.0D-6\r\n"
        "          ENDIF\r\n"
        "       ENDDO"
    )
    content, nsub = re.subn(pattern, replacement, content, count=1)
    if nsub == 0:
        raise RuntimeError("No se encontro el bloque de limites L(I)/U(I) a reemplazar")

    with open(dest, "w") as f:
        f.write(content)

ejecutables = {}
print("\nCOMPILANDO FILM0")

for nobs in nobs_unicos:
    film_patch = os.path.join(output_dir, f"FILM0_{nobs}.f")
    exe_path   = os.path.join(output_dir, f"film0_{nobs}.exe")
    patch_film0(film0_source, film_patch, nobs, ESPESOR)
    subprocess.run(["gfortran", "-O2", "-w", "-o", exe_path, film_patch], check=True)
    ejecutables[nobs] = exe_path
    print(f"  NOBS={nobs} -> OK")

# =========================================================
# EJECUTAR FILM0 EN CADA VENTANA
# =========================================================
print("\nEJECUTANDO FILM0")

suma_n = np.zeros(len(lam_exp))
suma_k = np.zeros(len(lam_exp))
suma_peso = np.zeros(len(lam_exp))

idx_lookup = {lam: i for i, lam in enumerate(lam_exp)}

for folder in ventanas:
    nobs      = nobs_por_ventana[folder]
    centro    = centro_por_ventana[folder]
    medio_ancho = medio_ancho_por_ventana[folder]
    exe_local = os.path.join(folder, "film0_exec.exe")
    shutil.copy(ejecutables[nobs], exe_local)

    t0 = time.time()
    result = subprocess.run([exe_local], cwd=folder,
                            capture_output=True, text=True)
    dt = time.time() - t0

    nombre = os.path.basename(folder)
    datank = os.path.join(folder, "DATANK1_DY")

    if not os.path.exists(datank):
        print(f"  {nombre} | FALLO sin resultado")
        continue

    largest_func, iterations = None, None
    for line in result.stdout.splitlines():
        if "LARGEST FUNCTION" in line:
            try: largest_func = float(line.split(":")[-1].strip())
            except: pass
        if "ITERATIONS" in line:
            try: iterations = int(line.split(":")[-1].strip())
            except: pass

    diag = f" | t={dt:.1f}s"
    if largest_func is not None: diag += f" | error={largest_func:.2e}"
    if iterations is not None: diag += f" | iter={iterations}"
    print(f"  {nombre} | OK{diag}")

    data = np.loadtxt(datank)
    lam_w, n_w, k_w = data[:, 0], data[:, 1], data[:, 2]

    dist = np.abs(lam_w - centro)
    peso = np.clip(1.0 - dist / medio_ancho, 0.0, None)

    for l, n, k, w in zip(lam_w, n_w, k_w, peso):
        idx = idx_lookup.get(round(l, 6))
        if idx is None:
            idx = int(np.argmin(np.abs(lam_exp - l)))
        suma_n[idx]    += w * n
        suma_k[idx]    += w * k
        suma_peso[idx] += w

# =========================================================
# CALCULAR PROMEDIO PONDERADO FINAL
# =========================================================
print("\nCALCULANDO PROMEDIO PONDERADO POR PUNTO")

mask_validos = suma_peso > 0
lam_final = lam_exp[mask_validos]
n_final   = suma_n[mask_validos] / suma_peso[mask_validos]
k_final   = suma_k[mask_validos] / suma_peso[mask_validos]

from scipy.signal import savgol_filter

WIN_SMOOTH_N = 11
WIN_SMOOTH_K = 31
POLY_ORDER   = 2

if len(n_final) > WIN_SMOOTH_N:
    n_final = savgol_filter(n_final, WIN_SMOOTH_N, POLY_ORDER)
if len(k_final) > WIN_SMOOTH_K:
    k_final = savgol_filter(k_final, WIN_SMOOTH_K, POLY_ORDER)

k_final = np.clip(k_final, 0.0, None)

print(f"\nSuavizado aplicado: n (ventana={WIN_SMOOTH_N}), k (ventana={WIN_SMOOTH_K})")

with open(tabla_final, "w") as f:
    for l, n, k in zip(lam_final, n_final, k_final):
        f.write(f"{l:.6f}   {n:.6f}   {k:.6f}\n")

print(f"\nTABLA FINAL: {tabla_final}  ({len(lam_final)} puntos, malla experimental completa)")

# =========================================================
# LIMPIAR CARPETAS DE VENTANAS (ya no se necesitan)
# =========================================================
shutil.rmtree(output_dir, ignore_errors=True)
print(f"Carpetas temporales de ventanas eliminadas ({output_dir})")

elapsed = time.time() - start_time
print(f"\nPipeline terminado en {elapsed:.2f} segundos")