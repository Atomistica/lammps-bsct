#! /usr/bin/env python

from StringIO import StringIO

import numpy as np

###

def loaddump(fn, col):
	data = []
	f = open(fn)
	l = f.readline()
	while l:
		for i in range(9):
			l = f.readline()
		for i in range(2):
			data += [float(l.split()[col])]
			l = f.readline()
	return np.array(data)

###

f = open('OUT')
l = f.readline()
while l.split()[0] != 'Step':
    l = f.readline()
s = StringIO()
l = f.readline()
while l.split()[0] != 'Loop':
    s.write(l)
    l = f.readline()
step, ect, ecoul, etot = np.loadtxt(StringIO(s.getvalue()), usecols=[0, 4, 6, 3], unpack=True)
step = np.array(step, dtype=int)
types = np.array(loaddump('dump.custom', 1), dtype=int)
z = loaddump('dump.custom', 4)
dist = z[types==2]-z[types==1]
charges = loaddump('dump.custom', 5)
charges1 = charges[types==1]
charges2 = charges[types==2]

###

# Charges should have equal magnitude but opposite sign

assert np.all(np.abs(charges1+charges2) < 1e-9)

###

X = 1.0
U = 0.1
V = 0.1
p = 1

# Check Coulomb energy

ecoul_check = charges1*charges2/dist

assert np.all(np.abs(ecoul-ecoul_check) < 1e-3)

# Check energy from charge-transfer model

ect_check  = X*(charges1-charges2)
ect_check += 0.5*U*charges1**2 + 0.5*U*charges2**2
ect_check += 0.5*V*np.abs(charges1)**p + 0.5*V*np.abs(charges2)**p

assert np.all(np.abs(ect-ect_check) < 1e-3)

# Check charges
# Total energy:    -q**2/dist + 2*X*q + U*q**2 + V*|q|
# Derivative:      -2*q/dist + 2*X + 2*U*q + V*sign(q)
# Equilibrium at:  q = (2*X+V)/(2.0/dist-2*U)

charges_check = (2*X+V)/(2.0/dist-2*U)

assert np.all(np.abs(charges1-charges_check) < 1e-3)
