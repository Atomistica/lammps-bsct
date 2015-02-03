LAMMPS-BSCT
===========

The files in this directory are a user-contributed package for LAMMPS.

This package was created by Tommi Järvi and Lars Pastewka. Please contact Lars
Pastewka (lars.pastewka@kit.edu) for questions and suggestions.

Lars Pastewka 
Karlsruhe Institute of Technology (KIT) 
Institute for Applied Materials (IAM) 
Kaiserstraße 12, 76131 Karlsruhe 
e-mail: lars.pastewka@kit.edu

PACKAGE DESCRIPTION
-------------------

This is a LAMMPS (http://lammps.sandia.gov/) fix style implementing the band
structure based charge transfer model described in Pastewka et al., Phys. Rev. B
83, 165418 (2011) (http://dx.doi.org/10.1103/PhysRevB.83.165418). It can also be used for simpler charge
transfer schemes.

It is necessary to use the provided Coulomb styles. Syntax:

> pair_style cout/cut/bsct cutoff

Or:

> pair_style coul/long/bsct cutoff  
> pair_coeff * *  
> kspace_style pppm/bsct prec  

Parameters are equivalent to coul/cut, coul/long and pppm. */bsct styles provide
the electrostatic potential required during the charge equilibration process.

The bsct fix takes care of actually equilibrating the charges. Syntax:

> fix ID group-ID bsct type atom-type X V U p keyword values

There needs to be a "type" specification for every atomic type that is part of
group-ID. The charge equilibration finds the charges q_i for each atom-type specified in the fix that minimizes the total energy functional

E = E_coul + sum_i ( X*q_i + 0.5*U*q_i**2 + V*q_i**p ),

where E_coul is the Coulomb interaction energy of the full system. The parametrs
X, V, U and p must be given for each atom type.

Example:

> fix 1 all bsct type 1 1.0 2.0 2.0 2.0 type 2 -1.0 2.0 2.0 2.0

In addition, the fix takes optional keywords.

keyword = nevery or prec or max_iter or qtot or history or beta or beta_mul or beta_max or mode or log
* nevery N: equilibrate every N steps (default N=1)
* prec P: stop optimization if maximum variation in chemical potential drops
    below P (default P=0.05)
* max_iter N: stop optimization after N iterations if not converged
    (default N=1000)
* qtot Q: total charge of subsystem (default Q=0.0)
* history N: history for Anderson mixer (default N=3)
* beta B: Anderson mixing factor (default B=0.1)
* beta_mul M: multiply Anderson mixing factor by M if energy goes down. Reset 
    to B if energy goes up (default M=1.0)
* beta_max B: don't increase Anderson mixing factor to a value larger than B
    (default B=0.1)
* mode M: Anderson mode (default M=2)
* log L: log verbosity (default L=0)

INSTALLATION
------------

In your LAMMPS src directory type:

> git clone https://github.com/Atomistica/lammps-bsct.git USER-BSCT  
> make yes-user-bsct

Note that the kspace package needs to be included and that the GNU
Scientific Library (GSL) is needed to compile this package.

Please run the tests located in the "tests" subdirectory after compiling. The
script "run_tests.sh" needs to be executed with the LAMMPS executable as a
parameter.

OTHER NOTES
-----------

* 2011-12-19  
  Initial release of package.

* 2015-02-03  
  Ported to recent (30 Jan 2015) version of LAMMPS.
