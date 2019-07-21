#PBS -N quick_file_convert -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00

python3.4 file_conversion.py


