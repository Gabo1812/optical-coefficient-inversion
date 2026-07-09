import numpy as np
import os
import shutil
import subprocess
import re
import time
import sys
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

start_time = time.time()

# =========================================================
# CONFIGURACION — parametros que el usuario puede cambiar
# =========================================================
MUESTRA = "S2"      # nombre de la muestra (S2 o S5)
ESPESOR = 95.0      # espesor de la pelicula en nm

VENTANA = 25        # medio-ancho de cada ventana espectral en nm
PASO    = 10        # separacion entre centros de ventana en nm
                    # (PASO < VENTANA genera solapamiento entre ventanas)

BASE = r"C:\Users\Gabo\Computacional\LabAvanzadoII\WO3"
archivo_exp  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_Dy_T")
archivo_nk   = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_DATANK0_DY")
film0_source = os.path.join(BASE, "FILM0.f")
tabla_final  = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_tabla_final_nk.txt")
output_dir   = os.path.join(BASE, "WO3_fortran", f"{MUESTRA}_ventanas_tmp")

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
    raise ValueError("Las longitudes de onda de Dy_T y DATANK0_DY no coinciden")

lam_min, lam_max = lam_exp.min(), lam_exp.max()
centros_todos = list(range(int(lam_min), int(lam_max) + 1, PASO))

if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)

# =========================================================
# CREAR VENTANAS
# =========================================================
print(f"\nCREANDO VENTANAS | ventana=+-{VENTANA} nm | paso={PASO} nm")

ventanas           = []
nobs_por_ventana   = {}
centro_por_ventana = {}

for i, centro in enumerate(centros_todos, 1):
    folder = os.path.join(output_dir, f"win_{i:03d}_c{centro}")
    os.makedirs(folder, exist_ok=True)

    mask  = np.abs(lam_exp - centro) <= VENTANA
    lam_w = lam_exp[mask]
    T_w   = T_exp[mask]

    if len(lam_w) == 0:
        continue

    n_w = PchipInterpolator(lam_nk, n0_all)(lam_w)
    k_w = np.clip(PchipInterpolator(lam_nk, k0_all)(lam_w), 0.0, None)

    with open(os.path.join(folder, "Dy_T"), "w") as f:
        for l, t in zip(lam_w, T_w):
            f.write(f"{l:.6f} {t:.6f}\n")
    with open(os.path.join(folder, "DATANK0_DY"), "w") as f:
        for l, n, k in zip(lam_w, n_w, k_w):
            f.write(f"{l:.6f} {n:.6f} {k:.6f}\n")

    ventanas.append(folder)
    nobs_por_ventana[folder]   = len(lam_w)
    centro_por_ventana[folder] = centro

print(f"  Total de ventanas: {len(ventanas)}")

# =========================================================
# PARCHEAR Y COMPILAR FILM0 POR NOBS UNICO
# =========================================================
def patch_film0(source, dest, nobs, thickness):
    """
    Modifica FILM0.f para usar:
    - El NOBS correcto para esta ventana
    - El espesor de la pelicula indicado
    - Sustrato de cuarzo amorfo (Malitson) en lugar de vidrio soda-lime
    - Limites de busqueda por punto anclados a la semilla de PUMA:
        n: [semilla - 0.03, semilla + 0.03]
        k: [0.5*semilla, 2.0*semilla] si k > 1e-8, sino [0, 1e-6]
    """
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
        raise RuntimeError("No se encontro el bloque de limites en FILM0.f")

    with open(dest, "w") as f:
        f.write(content)

nobs_unicos = sorted(set(nobs_por_ventana.values()))
print(f"\nCOMPILANDO FILM0 | NOBS unicos: {nobs_unicos}")

ejecutables = {}
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

suma_n    = np.zeros(len(lam_exp))
suma_k    = np.zeros(len(lam_exp))
suma_peso = np.zeros(len(lam_exp))
idx_lookup = {lam: i for i, lam in enumerate(lam_exp)}

for folder in ventanas:
    nobs    = nobs_por_ventana[folder]
    centro  = centro_por_ventana[folder]
    exe_local = os.path.join(folder, "film0_exec.exe")
    shutil.copy(ejecutables[nobs], exe_local)

    t0 = time.time()
    result = subprocess.run([exe_local], cwd=folder,
                            capture_output=True, text=True)
    dt = time.time() - t0

    datank = os.path.join(folder, "DATANK1_DY")
    if not os.path.exists(datank):
        print(f"  {os.path.basename(folder)} | FALLO | t={dt:.1f}s")
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
    if iterations  is not None: diag += f" | iter={iterations}"
    print(f"  {os.path.basename(folder)} | OK{diag}")

    data = np.loadtxt(datank)
    lam_w, n_w, k_w = data[:, 0], data[:, 1], data[:, 2]

    # Peso triangular: maximo en el centro, cero en el borde
    peso = np.clip(1.0 - np.abs(lam_w - centro) / VENTANA, 0.0, None)

    for l, n, k, w in zip(lam_w, n_w, k_w, peso):
        idx = idx_lookup.get(round(l, 6))
        if idx is None:
            idx = int(np.argmin(np.abs(lam_exp - l)))
        suma_n[idx]    += w * n
        suma_k[idx]    += w * k
        suma_peso[idx] += w

# =========================================================
# PROMEDIO PONDERADO Y SUAVIZADO FINAL
# =========================================================
print("\nCALCULANDO RESULTADO FINAL")

mask_validos = suma_peso > 0
lam_final = lam_exp[mask_validos]
n_final   = suma_n[mask_validos] / suma_peso[mask_validos]
k_final   = suma_k[mask_validos] / suma_peso[mask_validos]

n_final = savgol_filter(n_final, window_length=11, polyorder=2)
k_final = savgol_filter(k_final, window_length=31, polyorder=2)
k_final = np.clip(k_final, 0.0, None)

with open(tabla_final, "w") as f:
    for l, n, k in zip(lam_final, n_final, k_final):
        f.write(f"{l:.6f}   {n:.6f}   {k:.6f}\n")

shutil.rmtree(output_dir, ignore_errors=True)

print(f"\nTABLA FINAL: {tabla_final}  ({len(lam_final)} puntos)")
print(f"Pipeline terminado en {time.time()-start_time:.2f} segundos")