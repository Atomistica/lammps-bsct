/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/long/bsct,PairCoulLongBSCT)

#else

#ifndef LMP_PAIR_COUL_LONG_BSCT_H
#define LMP_PAIR_COUL_LONG_BSCT_H

#include "pair_coul_long.h"

namespace LAMMPS_NS {

class PairCoulLongBSCT : virtual public PairCoulLong {
 public:
  PairCoulLongBSCT(class LAMMPS *);
  virtual ~PairCoulLongBSCT() {}

  void reset_g_ewald();
  void phi(double &, double *);
};

}

#endif
#endif
