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

PairStyle(lj/charmmfsw/coul/long/bsct,PairLJCharmmfswCoulLongBSCT)

#else

#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_BSCT_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_BSCT_H

//0#include "pair_coul_long.h"
#include "pair_lj_charmmfsw_coul_long.h" //1

namespace LAMMPS_NS {

//0class PairLJCharmmfswCoulLongBSCT : virtual public PairCoulLong {
class PairLJCharmmfswCoulLongBSCT : virtual public PairLJCharmmfswCoulLong { //1
 public:
  PairLJCharmmfswCoulLongBSCT(class LAMMPS *);
  virtual ~PairLJCharmmfswCoulLongBSCT() {}

  void reset_g_ewald();
  void compute_potential(double &, double *);
protected: //1
double **scale; //1
};

}

#endif
#endif
