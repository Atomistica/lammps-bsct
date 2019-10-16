#! /usr/bin/env python

from StringIO import StringIO
import sys
import numpy as np

# options, parameters
verbose = True
debug = False

substrate_type  = 1
ion_type        = 2

dump_file = 'dump.custom'
log_file  = 'lammps.out'

# BSCT style:
# type 1 1.0 0.0 1.0 1.0 type 2 -1.0 0.0 1.0 1.0

X = 1.0
V = 0.0
U = 1.0
p = 1.0

###

def loaddump(fn, col):
    global verbose
    counter = 0
    data = []
    f = open(fn)
    l = f.readline()
    while l:
        while ( len( l.split() ) < 2 ) or not ( l.split()[0] == 'ITEM:' and l.split()[1] == 'ATOMS' ):
            if debug: sys.stderr.write( "Line {:06d}: '{}', ~ {} fields".format(counter,l.strip(),len(l.split())))
            l = f.readline()
            counter += 1

        if debug: sys.stderr.write( "Line {:06d}: '{}', ~ {} fields".format(counter,l.strip(),len(l.split())))
        if debug: sys.stderr.write( "Now reading atom data")
        l = f.readline()
        counter +=1

        d = []
        while l and (len( l.split() ) > 0) and (l.split()[0] != 'ITEM:'):
            if debug: sys.stderr.write( "Line {:06}: '{}', ~ {} fields".format(counter,l.strip(),len(l.split())))
            d += [float(l.split()[col])]
            l = f.readline()
            counter += 1 
        data += [d]

    if debug: sys.stderr.write( "Read {:d} lines.".format(counter))
    return np.array(data)

###

# read log fil (thermo output)
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

# read dump file
# dump style: id mol type x y z q
mol   = np.array(loaddump('dump.custom', 3), dtype=int)
if verbose: sys.stderr.write( "mol.shape{}".format(mol.shape))
types = np.array(loaddump('dump.custom', 2), dtype=int)
if verbose: sys.stderr.write( "types.shape{}".format(types.shape))
z = loaddump('dump.custom', 5)
if verbose: sys.stderr.write( "z.shape{}".format(z.shape))
#dist = z[types==2]-z[types==1]
charges = loaddump('dump.custom', 6)
if verbose: sys.stderr.write( "charges.shape{}".format(charges.shape))

charges_substrate = np.array(
  [ charges[i,types[i,:]==substrate_type] for i in range(charges.shape[0]) ] )  # AU
if verbose: sys.stderr.write( "charges_substrate.shape({})".format(charges_substrate.shape))

charges_ion = np.array(
  [ charges[i,types[i,:]==ion_type] for i in range(charges.shape[0]) ] )  # SOD
if verbose: sys.stderr.write( "charges_ion.shape({})".format(charges_ion.shape))

charge_substrate = np.sum(charges_substrate, axis=1)
if verbose: sys.stderr.write( "charge_substrate.shape({})".format(charge_substrate.shape))
charge_ion = np.sum(charges_ion, axis=1)
if verbose: sys.stderr.write( "charge_ion.shape({})".format(charge_ion.shape))

charge = np.vstack((charge_substrate,charge_ion)).T
if verbose: sys.stderr.write( "charge.shape({})".format(charge.shape))
print "charge_substrate charge_ion"
for i in range(charge.shape[0]):
    print "{:15.14e} {:15.14e}".format(*charge[i,:])

# if verbose:
#        print 'Coulomb energy error:', (ecoul+elong-ecoul_check)/ecoul_check

###

# Charges should have equal magnitude but opposite sign

# assert np.all(np.abs(charges1+charges2) < 1e-9)

###


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
