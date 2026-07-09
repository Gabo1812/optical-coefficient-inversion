import numpy as np
import os
import shutil
import subprocess
import re
import time

start_time = time.time()

# =========================================================
# CONFIGURACIÓN
# =========================================================
espesor_film = 167.0
archivo_exp  = "Dy_T"
archivo_nk   = "DATANK0_DY"

film0_source = "ITO_fortran/FILM0.f"

output_dir  = "ITO_ventanas"
tabla_final = "tabla_final_nk.txt"

# Solapamiento: paso = ventana/2
ventana = 40
centros = list(range(450, 1001, 20))  # centros cada 20 nm

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
    raise ValueError("λ de Dy_T y DATANK0_DY no coinciden")

os.makedirs(output_dir, exist_ok=True)

# =========================================================
# CREAR VENTANAS
# =========================================================

print("\nCREANDO VENTANAS")

ventanas         = []
nobs_por_ventana = {}
diagnosticos     = {}  # {folder: {"error": float, "iter": int}}

for i, centro in enumerate(centros, 1):

    folder = os.path.join(output_dir, f"win_{i:03d}")
    os.makedirs(folder, exist_ok=True)

    mask  = np.abs(lam_exp - centro) <= ventana
    lam_w = lam_exp[mask]
    T_w   = T_exp[mask]
    n0_w  = n0_all[mask]
    k0_w  = k0_all[mask]

    if len(lam_w) == 0:
        continue

    # ---- VENTANA ASIMÉTRICA EN BORDES:
    # Si la ventana tiene menos puntos que una ventana central completa,
    # extender el lado que SÍ tiene datos para compensar.
    nobs_central = int(2 * ventana / np.min(np.diff(lam_exp)) + 1)

    if len(lam_w) < nobs_central:
        lado_izq = centro - lam_exp[0]   # cuánto hay a la izquierda
        lado_der = lam_exp[-1] - centro  # cuánto hay a la derecha

        if lado_izq < ventana:
            # borde izquierdo — extender hacia la derecha
            extra = nobs_central - len(lam_w)
            nueva_lam_max = lam_w[-1] + extra * np.min(np.diff(lam_exp))
            nueva_lam_max = min(nueva_lam_max, lam_exp[-1])
            mask = (lam_exp >= lam_exp[0]) & (lam_exp <= nueva_lam_max)
        else:
            # borde derecho — extender hacia la izquierda
            extra = nobs_central - len(lam_w)
            nueva_lam_min = lam_w[0] - extra * np.min(np.diff(lam_exp))
            nueva_lam_min = max(nueva_lam_min, lam_exp[0])
            mask = (lam_exp >= nueva_lam_min) & (lam_exp <= lam_exp[-1])

        lam_w = lam_exp[mask]
        T_w   = T_exp[mask]
        n0_w  = n0_all[mask]
        k0_w  = k0_all[mask]

    with open(os.path.join(folder, "Dy_T"), "w") as f:
        for l, t in zip(lam_w, T_w):
            f.write(f"{l:.6f} {t:.6f}\n")

    with open(os.path.join(folder, "DATANK0_DY"), "w") as f:
        for l, n, k in zip(lam_w, n0_w, k0_w):
            f.write(f"{l:.6f} {n:.6f} {k:.6f}\n")

    nobs = len(lam_w)

    # Formato consistente con el parser de read_nobs()
    with open(os.path.join(folder, "meta.txt"), "w", encoding="utf-8") as f:
        f.write(f"Ventana centrada en λ={centro} nm ±{ventana} nm\n")
        f.write(f"NOBS = {nobs}\n")
        f.write(f"Rango: {lam_w[0]:.1f} – {lam_w[-1]:.1f} nm\n")
        f.write(f"Recuerda cambiar PARAMETER (NOBS = NOBS) en FILM0.f\n")

    ventanas.append(folder)
    nobs_por_ventana[folder] = nobs

    print(f"  win_{i:03d} | centro = {centro} nm | NOBS = {nobs}")

# =========================================================
# DETECTAR NOBS ÚNICOS
# =========================================================

nobs_unicos = sorted(set(nobs_por_ventana.values()))
print("\nNOBS únicos encontrados:", nobs_unicos)

# =========================================================
# PARCHEAR Y COMPILAR FILM0 POR NOBS ÚNICO
# =========================================================

def patch_film0(source, dest, nobs, thickness):
    with open(source) as f:
        content = f.read()

    # Parchear NOBS (maneja NOBS2=2*NOBS en la misma línea)
    content = re.sub(
        r"(PARAMETER\s*\(\s*NOBS\s*=\s*)\d+(\s*,\s*NOBS2\s*=\s*2\s*\*\s*NOBS\s*\))",
        rf"\g<1>{nobs}\g<2>",
        content,
        flags=re.IGNORECASE
    )

    # Parchear espesor HF (maneja comentario inline con !)
    content = re.sub(
        r"(HF\s*=\s*)[0-9]+\.?[0-9]*D0",
        rf"\g<1>{thickness:.1f}D0",
        content
    )

    with open(dest, "w") as f:
        f.write(content)


ejecutables = {}
print(f"\nEspesor del film: {espesor_film} nm")
print("\nCOMPILANDO FILM0")

for nobs in nobs_unicos:

    film_patch = os.path.join(output_dir, f"FILM0_{nobs}.f")
    exe_path   = os.path.join(output_dir, f"film0_{nobs}.exe")

    patch_film0(film0_source, film_patch, nobs, espesor_film)

    subprocess.run(
        ["gfortran", "-O2", "-w", "-o", exe_path, film_patch],
        check=True
    )

    ejecutables[nobs] = exe_path
    print(f"  NOBS = {nobs} → OK")

# =========================================================
# EJECUTAR FILM0 EN CADA VENTANA
# =========================================================

print("\nEJECUTANDO FILM0")

for folder in ventanas:

    nobs      = nobs_por_ventana[folder]
    exe_local = os.path.join(folder, "film0_exec.exe")

    shutil.copy(ejecutables[nobs], exe_local)

    result = subprocess.run(
        [exe_local],
        cwd=folder,
        capture_output=True,
        text=True
    )

    nombre = os.path.basename(folder)
    datank = os.path.join(folder, "DATANK1_DY")
    estado = "✓" if os.path.exists(datank) else "✗ sin resultado"

    # Extraer métricas del stdout del Fortran
    largest_func  = None
    iterations    = None
    for line in result.stdout.splitlines():
        if "LARGEST FUNCTION" in line:
            try: largest_func = float(line.split(":")[-1].strip())
            except: pass
        if "ITERATIONS" in line:
            try: iterations = int(line.split(":")[-1].strip())
            except: pass

    diag = ""
    if largest_func is not None:
        diag += f" | error={largest_func:.2e}"
    if iterations is not None:
        diag += f" | iter={iterations}"

    print(f"  {nombre} | {estado}{diag}")
    diagnosticos[folder] = {"error": largest_func, "iter": iterations}

# =========================================================
# ADVERTENCIA: VENTANAS CON MAL AJUSTE
# =========================================================

errores_validos = {f: d["error"] for f, d in diagnosticos.items() if d["error"] is not None}
if errores_validos:
    error_medio = np.mean(list(errores_validos.values()))
    umbral      = error_medio * 10   # ventanas con error 10x mayor al promedio
    malas       = [f for f, e in errores_validos.items() if e > umbral]
    if malas:
        print(f"\n⚠ VENTANAS CON MAL AJUSTE (error > {umbral:.2e}):")
        for f in malas:
            print(f"  {os.path.basename(f)} | error = {errores_validos[f]:.2e}")
    else:
        print(f"\n✓ Todas las ventanas convergieron bien (error medio = {error_medio:.2e})")

# =========================================================
# EXTRAER PUNTO CENTRAL DE CADA VENTANA
# =========================================================

print("\nEXTRAYENDO RESULTADOS")

resultados = []

for folder in ventanas:

    datank = os.path.join(folder, "DATANK1_DY")
    if not os.path.exists(datank):
        print(f"  ✗ Sin resultado: {os.path.basename(folder)}")
        continue

    # Leer centro desde meta.txt
    lam0 = None
    with open(os.path.join(folder, "meta.txt"), encoding="utf-8") as f:
        for line in f:
            m = re.search(r"λ=(\d+)", line)
            if m:
                lam0 = float(m.group(1))
                break

    if lam0 is None:
        print(f"  ✗ No se encontró λ0 en meta.txt: {os.path.basename(folder)}")
        continue

    data = np.loadtxt(datank)
    lam  = data[:, 0]
    n    = data[:, 1]
    k    = data[:, 2]

    idx = np.argmin(np.abs(lam - lam0))
    resultados.append((lam[idx], n[idx], k[idx]))
    print(f"  {os.path.basename(folder)} | λ={lam0} nm → n={n[idx]:.4f}, k={k[idx]:.6f}")

# =========================================================
# GUARDAR TABLA FINAL
# =========================================================

resultados = sorted(resultados, key=lambda x: x[0])

with open(tabla_final, "w") as f:
    for lam, n, k in resultados:
        f.write(f"{lam:.6f}   {n:.6f}   {k:.6f}\n")

print(f"\nTABLA FINAL: {tabla_final}  ({len(resultados)} puntos)")

# =========================================================
# LIMPIAR VENTANAS TEMPORALES
# =========================================================

print("\nBORRANDO VENTANAS TEMPORALES")
shutil.rmtree(output_dir)

# =========================================================
# TIEMPO TOTAL
# =========================================================

elapsed = time.time() - start_time
print(f"\nPipeline terminado en {elapsed:.2f} segundos")