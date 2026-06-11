#!/bin/bash
#SBATCH --partition=serial
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="FORTRAN"
#SBATCH --output=PUMA_first.out.%j
#SBATCH --mail-user=jose.alvarezcastrillo@ucr.ac.cr
#SBATCH --mail-type=END,FAIL

# Ir al directorio de trabajo
cd /home/jose.alvarezcastrillo/labavanzado/Lab_Avanzado/FORTRAN

# Compilar el programa Fortran usando ruta completa
/opt/ohpc/pub/compiler/gcc/12.2.0/bin/gfortran -O2 -o FILMO FILM0.f

# Ejecutar el programa
./FILMO
