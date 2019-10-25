#!/usr/bin/env python
"""Evaluate BSCT tests"""

from StringIO import StringIO
import sys
import numpy as np

# options, parameters
default_verbose = True
default_debug   = False

default_cathode_type  = 1
default_anode_type    = 2
default_cation_type   = 3
default_anion_type    = 4

default_dump_file = 'dump.custom'
default_log_file  = 'log.lammps'

# BSCT style:

hubbard_U_in_eV     = 6.8 # eV
eV_to_kcal          = 3.82929e-23
eV_to_kcal_per_mole = 23.0605

joule_to_kcal        = 0.239e-3    # conversion from J to kcal (kcal / J)
delta_phi            = 0.001       # potential difference between electrodes (V)
faraday_constant     = 96485.33212 # C / mol

default_X = joule_to_kcal * faraday_constant * delta_phi / 2.0
default_U = hubbard_U_in_eV * eV_to_kcal_per_mole
default_V = 0.0
default_p = 1.0

### loaders

def loaddump(fn, col):
    global verbose, debug
    counter = 0
    data = []
    f = open(fn)
    l = f.readline()
    while l:
        #print(counter)
        while ( len( l.split() ) < 2 ) or not ( l.split()[0] == 'ITEM:' and l.split()[1] == 'ATOMS' ):
            if debug: sys.stderr.write( "Line {:06d}: '{}', ~ {} fields\n".format(counter,l.strip(),len(l.split())))
            l = f.readline()
            counter += 1

        if debug: sys.stderr.write( "Line {:06d}: '{}', ~ {} fields\n".format(counter,l.strip(),len(l.split())))
        if debug: sys.stderr.write( "Now reading atom data")
        l = f.readline()
        counter +=1

        d = []
        while l and (len( l.split() ) > 0) and (l.split()[0] != 'ITEM:'):
            if debug: sys.stderr.write( "Line {:06}: '{}', ~ {} fields\n".format(counter,l.strip(),len(l.split())))
            d += [float(l.split()[col])]
            l = f.readline()
            counter += 1 
        data += [d]

    if debug: sys.stderr.write( "Read {:d} lines.\n".format(counter))
    return np.array(data)

def loadlog(fn):
    # read log fil (thermo output)
    f = open(fn)
    l = f.readline()
    counter = 0
    while len(l.split()) == 0 or l.split()[0] != 'Step':
        l = f.readline()
    s = StringIO()
    l = f.readline()
    while l and (len(l.split()) == 0 or l.split()[0] != 'Loop'):
        # ignore empty lines and warnings
        #if len(l.split()) > 0 and l.split()[0] != 'WARNING:':
        #use only lines of correct length
        if len(l.split())==16:
          s.write(l)
        l = f.readline()
    return s

### evaluation

def eval(    
    cathode_type = default_cathode_type,
    anode_type   = default_anode_type,
    cation_type  = default_cation_type,
    anion_type   = default_anion_type,
    
    dump_file    = default_dump_file,
    log_file     = default_log_file,
        
    X = default_X,
    U = default_U,
    V = default_V,
    p = default_p,
    
    **kwargs ):    
    global verbose, debug

    if debug:
        sys.stderr.write( "Got {:d} unknown keyword arguments:\n".format(len(kwargs)) ) 
        for k,v in kwargs.items():
            sys.stderr.write( "  {} : {}\n".format(k,v) ) 

    s = loadlog(log_file)
    try:
      # thermo style:
      # Step Temp Press TotEng KinEng PotEng f_ct E_bond E_angle E_dihed E_impro E_pair E_vdwl E_coul E_long E_tail
      step, temp, press, etot, ekin, epot, ect, ebond, eangle, edihed, eimpro, epair, evdwl, ecoul, elong, etail = np.loadtxt(
        StringIO(s.getvalue()), usecols=np.arange(0,16,1,dtype=int), unpack=True)
    except:
      # reference thermo style (no BSCT):
      # Step Temp Press TotEng KinEng PotEng E_bond E_angle E_dihed E_impro E_pair E_vdwl E_coul E_long E_tail
      step, temp, press, etot, ekin, epot, ebond, eangle, edihed, eimpro, epair, evdwl, ecoul, elong, etail = np.loadtxt(
        StringIO(s.getvalue()), usecols=np.arange(0,15,1,dtype=int), unpack=True)

    step = np.array(step, dtype=int)
    
    # read dump file
    # dump style: id mol type x y z q
    mol   = np.array(loaddump(dump_file, 3), dtype=int)
    if verbose: sys.stderr.write( "mol.shape{}\n".format(mol.shape))
    types = np.array(loaddump(dump_file, 2), dtype=int)
    if verbose: sys.stderr.write( "types.shape{}\n".format(types.shape))
    z = loaddump(dump_file, 5)
    if verbose: sys.stderr.write( "z.shape{}\n".format(z.shape))
    
    charges = loaddump(dump_file, 6)
    if verbose: sys.stderr.write( "charges.shape{}\n".format(charges.shape))
    
    charges_cathode = np.array(
      [ charges[i,types[i,:]==cathode_type] for i in range(charges.shape[0]) ] )  # AU
    if verbose: sys.stderr.write( "charges_cathode.shape({})\n".format(charges_cathode.shape))
    
    z_cathode = np.array(
      [ z[i,types[i,:]==cathode_type] for i in range(z.shape[0]) ] )  # SOD
    if verbose: sys.stderr.write( "z_cathode.shape({})\n".format(z_cathode.shape))
    
    charges_anode = np.array(
      [ charges[i,types[i,:]==anode_type] for i in range(charges.shape[0]) ] )  # AU
    if verbose: sys.stderr.write( "charges_anode.shape({})".format(charges_anode.shape))
    
    z_anode = np.array(
      [ z[i,types[i,:]==anode_type] for i in range(z.shape[0]) ] )  # SOD
    if verbose: sys.stderr.write( "z_anode.shape({})\n".format(z_anode.shape))
    
    charges_cation = np.array(
      [ charges[i,types[i,:]==cation_type] for i in range(charges.shape[0]) ] )  # SOD
    if verbose: sys.stderr.write( "charges_cation.shape({})\n".format(charges_cation.shape))
    
    z_cation = np.array(
      [ z[i,types[i,:]==cation_type] for i in range(z.shape[0]) ] )  # SOD
    if verbose: sys.stderr.write( "z_cation.shape({})\n".format(z_cation.shape))
    
    charges_anion = np.array(
      [ charges[i,types[i,:]==anion_type] for i in range(charges.shape[0]) ] )  # CL
    if verbose: sys.stderr.write( "charges_anion.shape({})\n".format(charges_anion.shape))
    
    z_anion = np.array(
      [ z[i,types[i,:]==anion_type] for i in range(z.shape[0]) ] )  # CL
    if verbose: sys.stderr.write( "z_anion.shape({})\n".format(z_anion.shape))
    
    # charges upper and lower electrode
    charge_cathode = np.sum(charges_cathode, axis=1)
    if verbose: sys.stderr.write( "charge_cathode.shape({})\n".format(charge_cathode.shape))
    
    charge_anode = np.sum(charges_anode, axis=1)
    if verbose: sys.stderr.write( "charge_anode.shape({})\n".format(charge_anode.shape))
    
    # position lower and upper interface
    z_lower_interface = np.max(z_cathode, axis=1)
    if verbose: sys.stderr.write( "z_lower_interface.shape({})\n".format(z_lower_interface.shape))
    
    z_upper_interface = np.min(z_anode, axis=1)
    if verbose: sys.stderr.write( "z_upper_interface.shape({})\n".format(z_upper_interface.shape))
    
    # in case of more than > 1 ion
    charge_cation = np.sum(charges_cation, axis=1)
    if verbose: sys.stderr.write( "charge_cation.shape({})\n".format(charge_cation.shape))
    
    zcom_cation = np.mean(z_cation, axis=1) 
    if verbose: sys.stderr.write( "zcom_cation.shape({})\n".format(zcom_cation.shape))
    
    charge_anion = np.sum(charges_anion, axis=1)
    if verbose: sys.stderr.write( "charge_anion.shape({})\n".format(charge_anion.shape))
    
    zcom_anion = np.mean(z_anion, axis=1) 
    if verbose: sys.stderr.write( "zcom_anion.shape({})\n".format(zcom_anion.shape))
    
    dcom_cation = np.abs(zcom_cation - z_lower_interface)
    dcom_anion  = np.abs(zcom_anion - z_upper_interface)
    
    output = np.vstack((
      dcom_cation, dcom_anion,
      charge_cation, charge_anion,
      charge_cathode,charge_anode )).T
    if verbose: sys.stderr.write( "output.shape({})\n".format(output.shape))
    print "zcom_cation zcom_anion charge_cation charge_anion charge_cathode charge_anode"
    
    for i in range(output.shape[0]):
        print ("{:15.14e} "*output.shape[1]).format(*output[i,:])
    
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

### command line parser

def main():
    global verbose, debug
    import argparse

    # in order to have both:
    # * preformatted help text and ...
    # * automatic display of defaults
    class ArgumentDefaultsAndRawDescriptionHelpFormatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
      pass

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

    parser.add_argument('log_file', nargs='?',
        help='LAMMPS stdout log', type=str, default=default_log_file)
    parser.add_argument('dump_file', nargs='?',
        help='LAMMPS custom dump', type=str, default=default_dump_file)


    parser.add_argument('--cathode-type', metavar='INT',
        help='cathode atom type', type=int, default=default_cathode_type)
    parser.add_argument('--anode-type', metavar='INT',
        help='anode atom type', type=int, default=default_anode_type)
    parser.add_argument('--cation-type', metavar='INT',
        help='cation atom type', type=int, default=default_cation_type)
    parser.add_argument('--anion-type', metavar='INT',
        help='anion atom type', type=int, default=default_anion_type)

    # BSCT style:
    parser.add_argument('-X','--electrochemical-potential', dest='X', metavar='FLOAT',
        help="one elctrode's electrochemical potential, system treated symmetrical", 
        type=float, default=default_X)
    parser.add_argument('-U', '--hubbard-u','--chemical-hardness', dest='U', metavar='FLOAT',
        help='Hubbard U or chemical hardness', type=float, default=default_U)
    parser.add_argument('-V', metavar='FLOAT',
        help='band structure contribution', type=float, default=default_V)
    parser.add_argument('-p', metavar='FLOAT',
        help='exponent', type=float, default=default_p)

    parser.add_argument('--debug', default=False, required=False,
                           action='store_true', dest="debug", help='debug flag')
    parser.add_argument('--verbose', default=False, required=False,
                           action='store_true', dest="verbose", help='verbose flag')

    args = parser.parse_args()

    verbose = args.verbose
    debug = args.debug
    
    kwargs = dict((k,v) for k,v in vars(args).items() if k!="message_type")

    if debug:
        sys.stderr.write( "Parsed {:d} command line options and parameters:\n".format(len(kwargs)) ) 
        for k,v in kwargs.items():
            sys.stderr.write( "  {} : {}\n".format(k,v) ) 
    
    eval(**kwargs)

if __name__ == '__main__':
    main()
