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

PairStyle(lj/charmmfsw/coul/charmmfsh/bsct,PairLJCharmmfswCoulCharmmfshBSCT)

#else

#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_CHARMMFSH_BSCT_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_CHARMMFSH_BSCT_H

#include "pair_coul_cut.h"

namespace LAMMPS_NS {

class PairLJCharmmfswCoulCharmmfshBSCT : virtual public PairCoulCut {
 public:
  PairLJCharmmfswCoulCharmmfshBSCT(class LAMMPS *);
  virtual ~PairLJCharmmfswCoulCharmmfshBSCT() {};

  void compute_potential(double &ecoultot, double *phi);
};

}

#endif
#endif
