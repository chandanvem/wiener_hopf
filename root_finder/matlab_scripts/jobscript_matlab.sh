#!/bin/bash
#PBS -N DISP_RELATION
#PBS -j oe
#PBS -l nodes=1:ppn=24
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR

NODEFILE=$PBS_NODEFILE
PPN=$(cat $NODEFILE | wc -l)
echo $NODEFILE
echo $PPN

module load apps/matlab/R2021b-new

matlab -nodisplay -nodesktop -r root_finder > vlog.${PBS_JOBID}
#matlab -nodisplay -nodesktop -r slicer_interpolation_main_specific_dir > vlog.${PBS_JOBID}

#matlab -nodisplay -nodesktop -r slicer_interpolation_main_gp_4 > vlog.${PBS_JOBID}

