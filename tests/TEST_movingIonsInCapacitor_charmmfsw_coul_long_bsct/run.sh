#!/bin/bash
# sample: run three comparative tests on NEMO

module purge
module use /work/ws/nemo/fr_jh1130-bsct-0/modulefiles/
module load lammps/07Aug19

NTASKS=20
LMP_EXE=lmp_nemo_20191015_bsct_charmmfsw

mpirun -n ${NTASKS} ${LMP_EXE} -in lammps_Q_zero.in | tee lammps_Q_zero.out
#mpirun -n ${NTASKS} ${LMP_EXE} -in lammps_Q_opp.in | tee lammps_Q_opp.out
#mpirun -n ${NTASKS} ${LMP_EXE} -in lammps_Q_same.in | tee lammps_Q_same.out

python eval.py lammps_Q_zero.out dump.Q_zero.custom > eval_Q_zero.out
#python eval.py lammps_Q_opp.out  dump.Q_opp.custom  > eval_Q_opp.out
#python eval.py lammps_Q_same.out dump.Q_same.custom > eval_Q_same.out

gnuplot eval_Q_zero.plt
#gnuplot eval_Q_opp.plt
#gnuplot eval_Q_same.plt
