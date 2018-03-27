#PBS -l walltime=1024:00:00,nodes=1:ppn=12

cd ~/ERCOFTAC/in-static-out-massFlow/

mpirun -np 6 -machinefile $PBS_NODEFILE twoPhaseMixingCentralDyMFoam -parallel | tee -a log.piso


