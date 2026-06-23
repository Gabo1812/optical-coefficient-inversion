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
MUESTRA = "S2"
ESPESOR = 95.0

BASE = r"C:\Users\Gabo\Computacional\LabAvanzadoII\WO3"
archivo_exp  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_Dy_T")
archivo_nk   = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_DATANK0_DY_constante")
film0_source = os.path.join(BASE, "FILM0.f")

output_dir  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_ventanas_constante")
tabla_final = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_tabla_final_nk_constante.txt")

# Mismas ventanas que la corrida con semilla PUMA, para comparacion directa
VENTANA_TRANSICION = 25
VENTANA_PLANA      = 30
PASO_TRANSICION    = 10
PASO_PLANO         = 15
LIMITE_ZONA        = 500
centro_inicial     = 550

N_CONSTANTE = 2.0   # valor "al tanteo" - estimacion generica, sin informacion previa
K_CONSTANTE = 0.0   # se asume material transparente, sin saber nada del borde UV

# Limites ANCHOS y GENERICOS (no anclados a ninguna semilla, ya que la
# semilla es constante y no aporta informacion real sobre n,k locales)
N_MIN_GLOBAL = 1.5
N_MAX_GLOBAL = 2.5
K_MIN_GLOBAL = 0.0
K_MAX_GLOBAL = 0.01   # suficientemente amplio para capturar absorcion real

def ventana_para_centro(centro):
    return VENTANA_TRANSICION if centro < LIMITE_ZONA else VENTANA_PLANA

def paso_para_centro(centro):
    return PASO_TRANSICION if centro < LIMITE_ZONA else PASO_PLANO

# =========================================================
# GENERAR SEMILLA CONSTANTE (n y k fijos en toda la malla)
# =========================================================
print("\nGENERANDO SEMILLA CONSTANTE")
data_exp_ref = np.loadtxt(archivo_exp)
lam_ref = data_exp_ref[:, 0]

with open(archivo_nk, "w") as f:
    for l in lam_ref:
        f.write(f"{l:.6f} {N_CONSTANTE:.6f} {K_CONSTANTE:.6f}\n")

print(f"  Semilla constante: n={N_CONSTANTE}, k={K_CONSTANTE}, "
      f"{len(lam_ref)} puntos [{lam_ref.min():.0f}-{lam_ref.max():.0f}] nm")

# =========================================================
# CARGAR DATOS
# =========================================================
data_exp = np.loadtxt(archivo_exp)
lam_exp  = data_exp[:, 0]
T_exp    = data_exp[:, 1]

data_nk = np.loadtxt(archivo_nk)   # el que acabamos de generar arriba
lam_nk  = data_nk[:, 0]
n0_all  = data_nk[:, 1]   # constante
k0_all  = data_nk[:, 2]   # constante (0)

if not np.allclose(lam_exp, lam_nk):
    raise ValueError("lambda de Dy_T y semilla no coinciden")

lam_min, lam_max = lam_exp.min(), lam_exp.max()

def generar_centros(lam_min, lam_max, centro_inicial):
    centros = set()
    c = centro_inicial
    while c <= lam_max:
        centros.add(c); c += paso_para_centro(c)
    c = centro_inicial
    while c >= lam_min:
        centros.add(c); c -= paso_para_centro(c)
    return sorted(centros)

centros_todos = generar_centros(lam_min, lam_max, centro_inicial)

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# =========================================================
# PATCH FILM0 — limites ANCHOS y GENERICOS (no por punto/semilla)
# =========================================================
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
        f"          L(I)={N_MIN_GLOBAL}D0\r\n"
        f"          U(I)={N_MAX_GLOBAL}D0\r\n"
        "       ENDDO\r\n"
        "       DO I=NOBS+1,NOBS2\r\n"
        f"          L(I)={K_MIN_GLOBAL}D0\r\n"
        f"          U(I)={K_MAX_GLOBAL}D0\r\n"
        "       ENDDO"
    )
    content, nsub = re.subn(pattern, replacement, content, count=1)
    if nsub == 0:
        raise RuntimeError("No se encontro el bloque de limites a reemplazar")
    with open(dest, "w") as f:
        f.write(content)

# =========================================================
# CREAR VENTANAS
# =========================================================
print("\nCREANDO VENTANAS (semilla constante n=2.0, k=0)")
ventanas, nobs_por_ventana = [], {}
centro_por_ventana, medio_ancho_por_ventana = {}, {}

for i, centro in enumerate(centros_todos, 1):
    folder = os.path.join(output_dir, f"win_{i:03d}_c{centro}")
    os.makedirs(folder, exist_ok=True)
    medio_ancho = ventana_para_centro(centro)
    mask  = np.abs(lam_exp - centro) <= medio_ancho
    lam_w = lam_exp[mask]; T_w = T_exp[mask]
    if len(lam_w) == 0: continue

    n_w = np.full_like(lam_w, n0_all[0])   # constante
    k_w = np.full_like(lam_w, k0_all[0])   # constante (0)

    with open(os.path.join(folder, "Dy_T"), "w") as f:
        for l, t in zip(lam_w, T_w): f.write(f"{l:.6f} {t:.6f}\n")
    with open(os.path.join(folder, "DATANK0_DY"), "w") as f:
        for l, n, k in zip(lam_w, n_w, k_w): f.write(f"{l:.6f} {n:.6f} {k:.6f}\n")

    nobs = len(lam_w)
    ventanas.append(folder)
    nobs_por_ventana[folder] = nobs
    centro_por_ventana[folder] = centro
    medio_ancho_por_ventana[folder] = medio_ancho

print(f"  Total de ventanas: {len(ventanas)}")

nobs_unicos = sorted(set(nobs_por_ventana.values()))
print("NOBS unicos:", nobs_unicos)

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
# EJECUTAR
# =========================================================
print("\nEJECUTANDO FILM0")
suma_n = np.zeros(len(lam_exp)); suma_k = np.zeros(len(lam_exp)); suma_peso = np.zeros(len(lam_exp))
idx_lookup = {round(l, 6): i for i, l in enumerate(lam_exp)}

for folder in ventanas:
    nobs = nobs_por_ventana[folder]; centro = centro_por_ventana[folder]
    medio_ancho = medio_ancho_por_ventana[folder]
    exe_local = os.path.join(folder, "film0_exec.exe")
    shutil.copy(ejecutables[nobs], exe_local)

    t0 = time.time()
    result = subprocess.run([exe_local], cwd=folder, capture_output=True, text=True)
    dt = time.time() - t0

    nombre = os.path.basename(folder)
    datank = os.path.join(folder, "DATANK1_DY")
    if not os.path.exists(datank):
        print(f"  {nombre} | FALLO sin resultado | t={dt:.1f}s")
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
        suma_n[idx] += w * n; suma_k[idx] += w * k; suma_peso[idx] += w

mask_validos = suma_peso > 0
lam_final = lam_exp[mask_validos]
n_final   = suma_n[mask_validos] / suma_peso[mask_validos]
k_final   = suma_k[mask_validos] / suma_peso[mask_validos]
k_final   = np.clip(k_final, 0.0, None)

with open(tabla_final, "w") as f:
    for l, n, k in zip(lam_final, n_final, k_final):
        f.write(f"{l:.6f}   {n:.6f}   {k:.6f}\n")

print(f"\nTABLA FINAL: {tabla_final}  ({len(lam_final)} puntos)")
print(f"Pipeline terminado en {time.time()-start_time:.2f} segundos")