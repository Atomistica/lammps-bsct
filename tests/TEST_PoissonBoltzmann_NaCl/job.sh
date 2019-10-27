#! /bin/bash -x 
#MSUB -l nodes=4:ppn=20
#MSUB -l walltime=08:00:00
#MSUB -l pmem=4000mb
#MSUB -m abe
#MSUB -M moepse@yahoo.de 
#/## // MSUB -v MPIRUN_OPTIONS="--bind-to core --map-by socket:PE=5 -report-bindings"
#MSUB -v EXECUTABLE=/work/ws/nemo/fr_jh1130-bsct-0/local_modules/lammps/07Aug19/src/lmp_nemo
#MSUB -v EXEC_RUN_OPTIONS="-in"
#MSUB -v EXEC_SCRIPT=/work/ws/nemo/fr_jh1130-bsct-0/local_modules/lammps/07Aug19/src/USER-BSCT/tests/TEST_PoissonBoltzmann_NaCl/production.in 
#MSUB -N My500V_test

module use /work/ws/nemo/fr_jh1130-bsct-0/modulefiles/
module load lammps/07Aug19

TASK_COUNT=${MOAB_PROCCOUNT}
echo "${EXECUTABLE} running on ${MOAB_PROCCOUNT} cores with ${TASK_COUNT} MPI-tasks"
startexe="mpirun -np ${TASK_COUNT} ${EXECUTABLE} ${EXEC_RUN_OPTIONS} ${EXEC_SCRIPT}"
echo $startexe
exec $startexe