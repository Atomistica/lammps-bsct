# NaCl (Na+, Cl-) in vacuum between plate capacitor

* lammps.in is standard test input
* lammps_no_BSCT.in is reference test input with BSCT optimization switched off
* lammps_Q_*.in are three test variants for different initial charges on electrodes
* eval.py is generic evaluation script for Python 2, run 'python eval.py --help' to see options
* eval_Q_*.plt are accroding gnuplot scripts for generating plots from eval.py output
* run.sh is sample on how to run and postprocess these three tests on NEMO