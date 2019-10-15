#! /bin/bash

if [ -z "$1" ]; then
  echo "Syntax: run_tests.sh <executable> [-n<number-of-processes>]"
  exit 999
fi

if [ ! -e "$1" ]; then
  echo "Executable $1 does not exist."
  exit 1
fi

CMD=$(readlink -f $1)

# specivy env var MPICMD to use other mpi command than default "mpirun"
if [ -z "${MPICMD}" ]; then
  MPICMD="mpirun"
fi
if [ -z "${MPINPOPT}" ]; then
  MPICMD="-np"
fi

np=""
OPTIND=2
while getopts ":n:" opt; do
  case $opt in
    n)
      np=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      exit 1
      ;;
  esac
done

nok=0
nfailed=0

if [ -n "$np" ]; then
  CMD="${MPICMD} ${MPINPOPT} $np $CMD"
fi
echo "${CMD}"
for i in TEST_*; do

  echo $i

  cd $i

  if [ -e run_test.sh ]; then
    sh run_test.sh $CMD > OUT
  else
  if [ -e lammps.in ]; then
    $CMD -in lammps.in > OUT
  else
    echo "Error: Don't know how to run test in directory $i"
    exit 999
  fi
  fi
  if [ $? -eq 0 ]; then
    python eval.py
    if [ $? -eq 0 ]; then
      echo ".ok."
      let nok=$nok+1
    else
      echo "     .failed."
      let nfailed=$nfailed+1
    fi
  else
    echo "     .failed."
    let nfailed=$nfailed+1
  fi

  cd ..

done

let ntot=$nok+$nfailed

echo "Ran $ntot tests; $nfailed failures, $nok successes."
if [ $nfailed -gt 0 ]; then
  echo "WARNING: Some tests failed."
fi
