#! /bin/sh

if [ -z "$1" ]; then
	echo "Syntax: run_test.sh <LAMMPS-executable>"
	exit 999
fi

for a0 in 3.0 4.0 5.0 7.0 10.0; do

  cat lammps.in | sed "s/%%a0/${a0}/g" | $1 > OUT.$a0
  mv log.lammps log.lammps.$a0
  mv dump.custom dump.custom.$a0

done
