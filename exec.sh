#!/bin/bash
module load mpi
athinput="/LyraShared/bcm2vn/athena/outputs/athinput.ion_lya"
problem_id=`cat $athinput | head -7 | tail -1`
basename=`grep -oP '(?<==\ ).*?(?=\ #)' <<< "$problem_id"`
mkdir "/LyraShared/bcm2vn/athena/outputs/$basename"
cd "./$basename"
filename="/LyraShared/bcm2vn/athena/outputs/$basename/$basename.txt"
cat /LyraShared/bcm2vn/athena/outputs/athinput.ion_lya >> $filename
mpirun -np 16 /LyraShared/bcm2vn/athena/bin/athena -i /LyraShared/bcm2vn/athena/outputs/athinput.ion_lya >> $filename
