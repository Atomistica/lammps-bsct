#! /bin/sh

if [ -z "$1" ]; then
	echo "Syntax: run_test.sh <LAMMPS-executable>"
	exit 999
fi

for a0 in 2.0 3.0 4.0 5.0; do

  cat lammps.in | sed "s/%%a0/${a0}/g" > lammps.in.tmp
  $1 -in lammps.in.tmp > OUT.$a0
  mv log.lammps log.lammps.$a0
  mv dump.custom dump.custom.$a0

done