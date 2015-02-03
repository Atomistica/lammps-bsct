#! /usr/bin/env python

from StringIO import StringIO

import numpy as np

import ase.io as io

###

def loadstate(logfn, dumpfn):
    f = open(logfn)
    l = f.readline()
    while len(l.split()) == 0 or l.split()[0] != 'Step':
        l = f.readline()
    l = f.readline().split()
    ect, ecoul, elong, lx = map(float, [l[4], l[6], l[5], l[8]])
    f.close()

    types, charges = np.loadtxt(dumpfn, skiprows=9, usecols=[1, 5], unpack=True)

    return lx, ect, ecoul, elong, types, charges

###

M = 1.747558 # Madelung constant of Na-Cl
nat = 64     # Number of atoms

X = 1.084
U = 5.0
V = 5.0
p = 2

###

for a0 in [3.0, 4.0, 5.0, 7.0, 10.0]:
    print '=== log.lammps.{}, dump.custom.{} ==='.format(a0, a0)

    lx, ect, ecoul, elong, types, charges = loadstate('log.lammps.{}'.format(a0), 'dump.custom.{}'.format(a0))
    charges1 = charges[types==1]
    charges2 = charges[types==2]

    # Charges should have equal magnitude but opposite sign

    assert np.all(np.abs(charges1+charges2) < 1e-3)

    # Check Coulomb energy

    r0 = lx/2 # Lattice constant
    charge = np.mean(charges1-charges2)/2
    ecoul_check = -charge**2/r0 * M * nat

    print 'Coulomb energy error:', (ecoul+elong-ecoul_check)/ecoul_check

    assert abs(ecoul+elong-ecoul_check) < 1e-3

    # Check energy from charge-transfer model

    ect_check  =  -X*charge
    ect_check += 0.5*U*charge**2
    ect_check += 0.5*V*charge**p
    ect_check *= nat

    print 'CT energy error:', (ect-ect_check)/ect_check

    assert abs(ect-ect_check) < 1e-4

    # Check charges
    # Total energy:    -q**2*M/r0 - X*q + 0.5*(U+V)*q**2
    # Derivative:      -2*q*M/r0 - X + (U+V)*q

    # Equilibrium at
    charge_check = X/(U+V-2*M/r0)

    print 'Charge error:', abs(charge-charge_check)/charge_check
    assert abs(charge-charge_check) < 0.01