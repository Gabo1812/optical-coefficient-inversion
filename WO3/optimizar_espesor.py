import numpy as np
import os
import shutil
import subprocess
import re
import time
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

sys.stdout.reconfigure(encoding="utf-8")
start_time = time.time()

# =========================================================
# CONFIGURACION POR MUESTRA
# =========================================================
MUESTRAS = {
    "S2": {"espesor_puma": 95.0},
    "S5": {"espesor_puma": 131.0},
}

CENTRO       = 550.0  # nm — punto de referencia estable (k≈0, sin ambiguedad)
MEDIO_ANCHO  = 20.0   # nm

BASE = r"C:\Users\Gabo\Computacional\LabAvanzadoII\WO3"
film0_source = os.path.join(BASE, "FILM0.f")

# =========================================================
# PATCH FILM0 (limites por punto anclados a semilla, igual que el pipeline principal)
# =========================================================
def patch_film0(source, dest, nobs, thickness):
    with open(source) as f:
        content = f.read()
    content = re.sub(
        r"(PARAMETER\s*\(\s*NOBS\s*=\s*)\d+(\s*,\s*NOBS2\s*=\s*2\s*\*\s*NOBS\s*\))",
        rf"\g<1>{nobs}\g<2>", content, flags=re.IGNORECASE)
    content = re.sub(
        r"(HF\s*=\s*)[0-9]+\.?[0-9]*D0",
        rf"\g<1>{thickness:.4f}D0", content)
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
        raise RuntimeError("No se encontro el bloque de limites a reemplazar")
    with open(dest, "w") as f:
        f.write(content)

# =========================================================
# PIPELINE COMPLETO PARA UNA MUESTRA
# =========================================================
def procesar_muestra(muestra, espesor_puma):
    print(f"\n{'='*60}\nProcesando {muestra} (espesor PUMA = {espesor_puma:.0f} nm)\n{'='*60}")

    espesores_candidatos = np.arange(espesor_puma - 25, espesor_puma + 10.5, 1)

    archivo_exp  = os.path.join(BASE, "WO3_fortran", f"{muestra}_Dy_T")
    archivo_nk   = os.path.join(BASE, "WO3_fortran", f"{muestra}_DATANK0_DY")

    output_dir = os.path.join(BASE, "WO3_fortran", f"{muestra}_barrido_espesor")
    out_fig     = os.path.join(BASE, "WO3_fortran", f"{muestra}_barrido_espesor.pdf")
    out_fig_png = os.path.join(BASE, "WO3_fortran", f"{muestra}_barrido_espesor.png")
    out_tabla  = os.path.join(BASE, "WO3_fortran", f"{muestra}_barrido_espesor.txt")

    # ── Cargar datos y construir la ventana en 550 nm ──────────────────────
    data_exp = np.loadtxt(archivo_exp)
    lam_exp  = data_exp[:, 0]; T_exp = data_exp[:, 1]

    data_nk = np.loadtxt(archivo_nk)
    lam_nk  = data_nk[:, 0]; n0_all = data_nk[:, 1]; k0_all = data_nk[:, 2]

    mask  = np.abs(lam_exp - CENTRO) <= MEDIO_ANCHO
    lam_w = lam_exp[mask]; T_w = T_exp[mask]
    n_w   = PchipInterpolator(lam_nk, n0_all)(lam_w)
    k_w   = np.clip(PchipInterpolator(lam_nk, k0_all)(lam_w), 0.0, None)
    NOBS  = len(lam_w)

    print(f"Ventana en {CENTRO} nm +-{MEDIO_ANCHO} nm: NOBS={NOBS}")
    print(f"Barriendo espesor: {espesores_candidatos.min():.0f} a "
          f"{espesores_candidatos.max():.0f} nm, paso 1 nm\n")

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # ── Barrido de espesores ────────────────────────────────────────────────
    resultados = []   # (espesor, error)

    for espesor in espesores_candidatos:
        folder = os.path.join(output_dir, f"esp_{espesor:.0f}")
        os.makedirs(folder, exist_ok=True)

        with open(os.path.join(folder, "Dy_T"), "w") as f:
            for l, t in zip(lam_w, T_w): f.write(f"{l:.6f} {t:.6f}\n")
        with open(os.path.join(folder, "DATANK0_DY"), "w") as f:
            for l, n, k in zip(lam_w, n_w, k_w): f.write(f"{l:.6f} {n:.6f} {k:.6f}\n")

        film_patch = os.path.join(folder, "FILM0.f")
        exe_path   = os.path.join(folder, "film0.exe")
        patch_film0(film0_source, film_patch, NOBS, espesor)
        subprocess.run(["gfortran", "-O2", "-w", "-o", exe_path, film_patch], check=True)

        result = subprocess.run([exe_path], cwd=folder, capture_output=True, text=True)

        largest_func = None
        for line in result.stdout.splitlines():
            if "LARGEST FUNCTION" in line:
                try: largest_func = float(line.split(":")[-1].strip())
                except: pass

        if largest_func is not None:
            resultados.append((espesor, largest_func))
            print(f"  espesor={espesor:6.1f} nm | error={largest_func:.3e}")
        else:
            print(f"  espesor={espesor:6.1f} nm | SIN RESULTADO")

    # ── Guardar tabla ────────────────────────────────────────────────────────
    resultados = np.array(resultados)
    espesores_arr, errores_arr = resultados[:, 0], resultados[:, 1]

    with open(out_tabla, "w") as f:
        for e, err in resultados:
            f.write(f"{e:.2f}   {err:.6e}\n")

    idx_min = np.argmin(errores_arr)
    espesor_optimo = espesores_arr[idx_min]
    error_optimo   = errores_arr[idx_min]

    print(f"\nEspesor optimo encontrado: {espesor_optimo:.1f} nm "
          f"(error={error_optimo:.3e})")
    print(f"Espesor de referencia PUMA: {espesor_puma:.1f} nm")
    print(f"Diferencia: {abs(espesor_optimo - espesor_puma):.1f} nm")

    # ── Figura (sin titulo, con la muestra como etiqueta en la esquina) ─────
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(espesores_arr, errores_arr, "o-", color="steelblue",
            linewidth=3.0, markersize=8)
    ax.axvline(espesor_puma, color="darkorange", linestyle="--", linewidth=2.5,
               label=f"Espesor PUMA ({espesor_puma:.0f} nm)")
    ax.axvline(espesor_optimo, color="crimson", linestyle=":", linewidth=2.5,
               label=f"Minimo SPG ({espesor_optimo:.0f} nm)")
    ax.set_yscale("log")
    ax.set_xlabel("Espesor (nm)", fontsize=19)
    ax.set_ylabel("Error de ajuste (LARGEST FUNCTION)", fontsize=19)
    ax.tick_params(labelsize=17)
    ax.legend(fontsize=15)
    ax.grid(True, alpha=0.3, which="both")

    # Recorte del eje x visible: se calcula y guarda el barrido completo en la
    # tabla, pero en la grafica se recortan solo 5 nm del extremo izquierdo
    # (la parte plana que menos aporta), dejando el resto del rango calculado
    # (-25 a +10 nm respecto a PUMA) tal cual.
    ax.set_xlim(espesor_puma - 20, espesor_puma + 10)

    ax.text(0.04, 0.95, muestra, transform=ax.transAxes,
            fontsize=17, fontweight="bold", color="black",
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="black", alpha=0.9))

    plt.tight_layout()
    plt.savefig(out_fig, bbox_inches="tight", dpi=300)
    plt.savefig(out_fig_png, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"\nGuardado: {out_tabla}")
    print(f"Guardado: {out_fig}")
    print(f"Guardado: {out_fig_png}")

    shutil.rmtree(output_dir, ignore_errors=True)

# =========================================================
# CORRER PARA TODAS LAS MUESTRAS
# =========================================================
for muestra, params in MUESTRAS.items():
    procesar_muestra(muestra, params["espesor_puma"])

elapsed = time.time() - start_time
print(f"\nPipeline terminado en {elapsed:.2f} segundos (ambas muestras)")