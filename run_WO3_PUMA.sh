#!/bin/bash
#SBATCH --partition=serial
#SBATCH --account=2026-i-fs0751
#SBATCH --qos=curso-cpu-long
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="WO3_PUMA"
#SBATCH --mail-user=jose.alvarezcastrillo@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/dev/null
#SBATCH --error=slurm-%x-%j.err

# =============================================================
# USO:
#   sbatch run_WO3_PUMA_1.sh S2 1   ← Barrido amplio  (INIT=0) — YA HECHO P1-P4
#   sbatch run_WO3_PUMA_1.sh S2 2   ← Refinamiento    (INIT=9) — PRÓXIMO para S2/S5
#   sbatch run_WO3_PUMA_1.sh S2 3   ← Ajuste fino     (INIT=9) — actualizar QUAD_3/DFIX_3/INFFIX_3
#   sbatch run_WO3_PUMA_1.sh S8 1   ← NUEVO barrido amplio S8  (INIT=0) — P1/P2 no convergieron
#
# Parámetros actualizados (segunda llamada):
#   LMIN=300 (antes 350), NOBS=150 (antes 100), INFLE=340-460 paso 20
#   QUAD_2 con valores reales de P4 (S2) y P6 (S5)
#
# Para Run 2 (INIT=9): el -sol.txt del mejor run anterior debe estar en SRC.
# Lanzar S2+S5 Run 2 simultáneamente (2 slots). S8 Run 1 cuando libere uno.
# =============================================================

SAMPLE_KEY="${1}"
RUN_NUM="${2}"

if [[ -z "$SAMPLE_KEY" || -z "$RUN_NUM" ]]; then
    echo "ERROR: uso: sbatch run_WO3_PUMA.sh <S2|S5|S8> <1|2|3>"
    exit 1
fi

# ----------------------------------------------------------
# Parámetros por muestra
# ----------------------------------------------------------
case "$SAMPLE_KEY" in
    S2)
        FNAME="WO3_S2_10sccm"
        # Run 1 — barrido amplio (YA HECHO: P1–P4, mejor P4)
        DMIN_1="0065"; DMAX_1="0150"; DSTEP_1=5
        INFMIN_1=0340; INFMAX_1=0460; INFSTEP_1=20
        N0INI=1.90; N0FIN=2.50; N0STEP=0.10
        NFINI=1.80; NFFIN=2.20; NFSTEP=0.10
        K0INI=0.00; K0FIN=0.30; K0STEP=0.02
        # Run 2 — refinamiento desde P4-sol.txt (INIT=9)
        QUAD_2="2.651887e-05"             # QE real de P4
        DMIN_2="0088"; DMAX_2="0102"; DSTEP_2=1
        INFMIN_2=0340; INFMAX_2=0460; INFSTEP_2=20
        # Run 3 — ajuste fino (actualizar QUAD_3/DFIX_3/INFFIX_3 con valores de Run 2)
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="<D_RUN2>"; INFFIX_3="<INF_RUN2>"
        ;;
    S5)
        FNAME="WO3_S5_7sccm"
        # Run 1 — barrido amplio (YA HECHO: P1–P6, mejor P6)
        DMIN_1="0100"; DMAX_1="0160"; DSTEP_1=5
        INFMIN_1=0340; INFMAX_1=0460; INFSTEP_1=20
        N0INI=1.90; N0FIN=2.50; N0STEP=0.10
        NFINI=1.80; NFFIN=2.20; NFSTEP=0.10
        K0INI=0.00; K0FIN=0.30; K0STEP=0.02
        # Run 2 — refinamiento desde P6-sol.txt (INIT=9)
        QUAD_2="5.030000e-05"             # QE real de P6
        DMIN_2="0126"; DMAX_2="0136"; DSTEP_2=1
        INFMIN_2=0340; INFMAX_2=0460; INFSTEP_2=20
        # Run 3 — ajuste fino (actualizar con valores de Run 2)
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="<D_RUN2>"; INFFIX_3="<INF_RUN2>"
        ;;
    S8)
        FNAME="WO3_S8_5sccm"
        # Run 1 — NUEVO barrido amplio (P1–P2 no convergieron, INIT=0)
        DMIN_1="0090"; DMAX_1="0120"; DSTEP_1=2
        INFMIN_1=0340; INFMAX_1=0460; INFSTEP_1=20
        N0INI=2.00; N0FIN=2.70; N0STEP=0.10
        NFINI=1.80; NFFIN=2.20; NFSTEP=0.10
        K0INI=0.05; K0FIN=0.40; K0STEP=0.05
        # Run 2 — actualizar tras Run 1
        QUAD_2="<QUAD_RUN1>"
        DMIN_2="<DMIN>"; DMAX_2="<DMAX>"; DSTEP_2=1
        INFMIN_2=0340; INFMAX_2=0460; INFSTEP_2=20
        # Run 3 — ajuste fino (actualizar con valores de Run 2)
        QUAD_3="<QUAD_RUN2>"
        DFIX_3="<D_RUN2>"; INFFIX_3="<INF_RUN2>"
        ;;
    *)
        echo "ERROR: muestra desconocida '$SAMPLE_KEY'. Usar S2, S5 o S8."
        exit 1
        ;;
esac

# ----------------------------------------------------------
# Parámetros PUMA comunes
# ----------------------------------------------------------
NLAYERS=4; SLAYER=2; SUBSTRATE=60; DATATYPE=T; NOBS=150
LMIN=0300; LMAX=1500

# ----------------------------------------------------------
SRC="/home/jose.alvarezcastrillo/labavanzado/Lab_Avanzado"
TMPDIR=$(mktemp -d)
echo "Muestra: $FNAME  |  Run: $RUN_NUM"
echo "Corriendo PUMA en: $TMPDIR"

cp $SRC/puma_mod           $TMPDIR/
cp $SRC/${FNAME}-dat.txt   $TMPDIR/
cp $SRC/${FNAME}-sol.txt   $TMPDIR/ 2>/dev/null

cd $TMPDIR
rm -f $SRC/${FNAME}-inf.txt

# ----------------------------------------------------------
# RUN 1 — Barrido amplio (INIT=0)
# n y k: rangos amplios definidos arriba
# ----------------------------------------------------------
if [[ "$RUN_NUM" == "1" ]]; then
    echo "=== RUN 1 — Barrido amplio (INIT=0, LMIN=$LMIN, NOBS=$NOBS) ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 3000 1e+100 0 \
        $DMIN_1 $DMAX_1 $DSTEP_1 \
        $INFMIN_1 $INFMAX_1 $INFSTEP_1 \
        $N0INI $N0FIN $N0STEP \
        $NFINI $NFFIN $NFSTEP \
        $K0INI $K0FIN $K0STEP

# ----------------------------------------------------------
# RUN 2 — Refinamiento (INIT=9)
# Requiere -sol.txt del Run 1 en SRC
# n y k: NO se pasan, PUMA los toma del -sol.txt
# ----------------------------------------------------------
elif [[ "$RUN_NUM" == "2" ]]; then
    echo "=== RUN 2 — Refinamiento (INIT=9, LMIN=$LMIN, NOBS=$NOBS) ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 5000 $QUAD_2 9 \
        $DMIN_2 $DMAX_2 $DSTEP_2 \
        $INFMIN_2 $INFMAX_2 $INFSTEP_2

# ----------------------------------------------------------
# RUN 3 — Ajuste fino (INIT=9, espesor e inflexión fijos)
# Requiere -sol.txt del Run 2 en SRC
# Actualizar QUAD_3, DFIX_3, INFFIX_3 con los valores del Run 2
# ----------------------------------------------------------
elif [[ "$RUN_NUM" == "3" ]]; then
    echo "=== RUN 3 — Ajuste fino ==="
    ./puma_mod $FNAME \
        $NLAYERS $SLAYER $SUBSTRATE $DATATYPE $NOBS \
        $LMIN $LMAX 50000 $QUAD_3 9 \
        $DFIX_3 $DFIX_3 1 \
        $INFFIX_3 $INFFIX_3 1

else
    echo "ERROR: run desconocido '$RUN_NUM'. Usar 1, 2 o 3."
    exit 1
fi

echo "=== PUMA finalizado ==="

# ----------------------------------------------------------
# Copiar resultados a SRC y limpiar
# ----------------------------------------------------------
cp ${FNAME}-sol.txt ${FNAME}-inf.txt $SRC/ 2>/dev/null
echo "Resultados copiados a $SRC"

rm -rf $TMPDIR
echo "Carpeta temporal eliminada."