#!/bin/bash
#SBATCH --partition=serial
#SBATCH --account=2026-i-fs0751
#SBATCH --qos=curso-cpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="PUMA_2"
#SBATCH --mail-user=jose.alvarezcastrillo@ucr.ac.cr
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/dev/null
#SBATCH --error=slurm-%x-%j.err

# Directorio de trabajo (donde está puma_mod, los datos y caerán los resultados)
SRC="/home/jose.alvarezcastrillo/labavanzado/Lab_Avanzado"

# Carpeta temporal para correr PUMA (se borra sola al terminar)
TMPDIR=$(mktemp -d)
echo "Corriendo PUMA en carpeta temporal: $TMPDIR"

# ----------------------------------------------------------
# Copiar ejecutable, datos y solución anterior a TMPDIR
# El -sol.txt es necesario para iniptype=9
# ----------------------------------------------------------
cp $SRC/puma_mod            $TMPDIR/
cp $SRC/ITO_dataset-dat.txt $TMPDIR/
cp $SRC/*-sol.txt           $TMPDIR/ 2>/dev/null

cd $TMPDIR

# ----------------------------------------------------------
# Borrar solo el -inf.txt anterior
# El -sol.txt NO se borra: PUMA lo necesita como punto de partida
# ----------------------------------------------------------
rm -f $SRC/*-inf.txt
echo "-inf.txt anterior eliminado."

# ----------------------------------------------------------
# Ejecutar PUMA
# ----------------------------------------------------------

# RUN 1: barrido amplio (iniptype=0)
#./puma_mod ITO_dataset 4 2 70 T 100 0535 1000 3000 1e+100 0 0166 0168 1 0850 0900 25 1.9 2.1 0.02 1.9 2.0 0.02 0.00 0.01 0.002

# RUN 2: refinamiento (iniptype=9 — requiere -sol.txt del run anterior)
#./puma_mod ITO_dataset 4 2 70 T 100 0535 1000 5000 3.885226e-05 9 0164 0166 0.1 0825 0875 25

# RUN 3: ajuste fino (iniptype=9)
./puma_mod ITO_dataset 4 2 70 T 100 0535 1000 500000 3.789647e-05 9 0164 0166 0.1 0825 0825 100

# ----------------------------------------------------------
# Copiar resultados a SRC y limpiar carpeta temporal
# -sol.txt sobreescribe el anterior con la nueva solución
# -inf.txt es el reporte final
# ----------------------------------------------------------
cp *-sol.txt *-inf.txt $SRC/ 2>/dev/null
echo "Resultados copiados a $SRC"

rm -rf $TMPDIR
echo "Carpeta temporal eliminada."
