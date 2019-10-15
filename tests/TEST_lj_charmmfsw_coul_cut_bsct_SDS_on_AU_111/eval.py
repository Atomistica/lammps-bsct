#! /usr/bin/env python

from StringIO import StringIO

import numpy as np

###

def loaddump(fn, col):
    data = []
    f = open(fn)
    l = f.readline()
    while l:
        while ( len( l.split() ) < 2 ) or ( ( l.split()[0] != 'ITEM:' ) and ( l.split()[1] != 'ATOMS' ) ):
            l = f.readline()

        l = f.readline()
        while l and (len( l.split() ) > 0) and (l.split()[0] != 'ITEM:'):
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
# thermo style:
# Step Temp Press TotEng KinEng PotEng f_ct E_bond E_angle E_dihed E_impro E_pair E_vdwl E_coul E_long E_tail
step, temp, press, etot, ekin, epot, ect, ebond, eangle, edihed, eimpro, epair, evdwl, ecoul, elong, etail = np.loadtxt(
  StringIO(s.getvalue()), usecols=np.arange(0,16,1,dtype=int), unpack=True)
step = np.array(step, dtype=int)
# dump style:
# id mol type x y z q
mol   = np.array(loaddump('dump.custom', 3), dtype=int)
types = np.array(loaddump('dump.custom', 2), dtype=int)
z = loaddump('dump.custom', 5)
#dist = z[types==2]-z[types==1]
charges = loaddump('dump.custom', 6)
charges1 = charges[types==8]  # AU
charges2 = charges[types==12] # SOD

# if verbose:
#        print 'Coulomb energy error:', (ecoul+elong-ecoul_check)/ecoul_check

###

# Charges should have equal magnitude but opposite sign

# assert np.all(np.abs(charges1+charges2) < 1e-9)

###

# BSCT style:
# type 8 1.0 0.0 1.0 1.0 type 12 -1.0 0.0 1.0 1.0

X = 1.0
V = 0.0
U = 1.0
p = 1.0

# Check Coulomb energy

# dummy, adapt for point charge above surface

# ecoul_check = charges1*charges2/dist

# assert np.all(np.abs(ecoul-ecoul_check) < 1e-3)

# Check energy from charge-transfer model

# ect_check  = X*(charges1-charges2)
# ect_check += 0.5*U*charges1**2 + 0.5*U*charges2**2
# ect_check += 0.5*V*np.abs(charges1)**p + 0.5*V*np.abs(charges2)**p

# assert np.all(np.abs(ect-ect_check) < 1e-3)

# Check charges
# Total energy:    -q**2/dist + 2*X*q + U*q**2 + V*|q|
# Derivative:      -2*q/dist + 2*X + 2*U*q + V*sign(q)
# Equilibrium at:  q = (2*X+V)/(2.0/dist-2*U)

# charges_check = (2*X+V)/(2.0/dist-2*U)

# assert np.all(np.abs(charges1-charges_check) < 1e-3)
