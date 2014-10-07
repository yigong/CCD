import numpy as np

class slurmGenerator():
    ''' a class to generate slurm script '''
    
    def __init__(self, fileNameIdx):
        fileName = 'eSample_%s.slm'%(int(fileNameIdx))
        self.fileObj = open(fileName, 'w')
        self.fileNameIdx = fileNameIdx
        
    def writeAll(self):
        
        
        lines = ['#!/bin/bash -I']
        lines.append('#Job name:')
        lines.append('#SBATCH --job-name=eSample_%s'%(self.fileNameIdx))
        
        lines.append('#Account:')
        lines.append('#SBATCH -A, --account=ac_etci')
        lines.append('#SBATCH -e, --error=output_link/%s.err'%(self.fileNameIdx))
        lines.append('#SBATCH -o, --output=output_link/%s.out'%(self.fileNameIdx))
        
        lines.append('#Partition')
        lines.append('#SBATCH --partition=lr2')
        lines.append('#QoS:')
        lines.append('#SBATCH --qos=lr_normal')
        lines.append('#Processors:')
        lines.append('#SBATCH --nodes=1')
        lines.append('#SBATCH --ntasks-per-node=12')
        
        lines.append('#Wall clock limit:')
        lines.append('#SBATCH --time=12:00:00')
        
        lines.append('source /global/scratch/ygzhang/g4_work/misc/env.sh')
        lines.append('cd /global/scratch/ygzhang/g4_work/eSample')
        lines.append('./E_in_Nitro eSample_%s.mac >output_link/eSample_%s.out &\n'\
                     %(self.fileNameIdx, self.fileNameIdx))
        lines.append('''FAIL=0
for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
done

echo $FAIL

if [ "$FAIL" == "0" ];
then
echo "YAY!"
else
echo "FAIL! ($FAIL)"
fi''')
        self.fileObj.write('\n'.join(lines))
        self.fileObj.close()