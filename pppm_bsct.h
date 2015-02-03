/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Tommi Järvi (Fh-IWM), Lars Pastewka (KIT)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/bsct,PPPMBSCT)

#else

#ifndef LMP_PPPMBSCT_H
#define LMP_PPPMBSCT_H

#include "lmptype.h"
#include "mpi.h"

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMBSCT : public PPPM {
 public:
  PPPMBSCT(class LAMMPS *, int, char **);
  virtual ~PPPMBSCT();

  virtual void phi(FFT_SCALAR *phi);           // Get electrostatic potential (phi)

 protected:
  virtual void phi_slabcorr(FFT_SCALAR *phi);  // Get slab correction to phi
};

}

#endif
#endif
