#PBS -q lr_batch
#PBS -A ac_etci
#PBS -N HT1
#PBS -l nodes=1:ppn=12:lr2
#PBS -l walltime=24:00:00
#PBS -o /global/scratch/acoffer/g4Data/ShellScriptOutputs/matHT_out24_1
#PBS -e /global/scratch/acoffer/g4Data/ShellScriptOutputs/matHT_err24_1
#PBS -M abcoffer@berkeley.edu
#PBS -m bae

/bin/bash
module load matlab/R2011b

cd /global/home/users/acoffer/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes

unset DISPLAY

matlab -nodisplay -nosplash -nodesktop -r HTparallel 

ssh n0000.scs00 "cd /global/home/users/acoffer/ETCI/simulation/acoff_matlab/ShellScripts/ ; HTMatlabScript1.sh"

