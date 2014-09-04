#PBS -q lr_batch
#PBS -A ac_etci
#PBS -N pB1
#PBS -l nodes=1:ppn=1:lr3
#PBS -l walltime=2:00:00
#PBS -o /global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/ShellScripts/matHT_out2hr_1
#PBS -e /global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/ShellScripts/matHT_err2hr_1
#PBS -M abcoffer@berkeley.edu
#PBS -m bae

/bin/bash
module load matlab/R2011b

#cd /global/home/users/acoffer/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes
#cd /global/home/users/acoffer/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes
cd /global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam

unset DISPLAY

matlab -nodisplay -nosplash -nodesktop -r RunCollection8.m

#ssh n0000.scs00 "cd /global/home/users/acoffer/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes/ ; HTMatlabScript1.sh"

