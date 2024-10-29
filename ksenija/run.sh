#!/bin/bash
#PBS -N parpack
cd $PBS_O_WORKDIR/
source /opt/intel/oneapi/setvars.sh
EXECDIR=/opt/intel/oneapi/mpi/latest/bin

sizes=(20 50)

for size in "${sizes[@]}"; do
   for ((run_num=1; run_num<4; run_num++))
   do
       input_file="${size}input_${run_num}"
       echo $size $run_num > "$input_file"
       $EXECDIR/mpirun -np 5 ./myparpack "$input_file" > disorder${run_num}.out 2> disorder${run_num}.err
   done
done
