/* ----------------------------------------------------------------------
//   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
//   http://lammps.sandia.gov, Sandia National Laboratories
//   Steve Plimpton, sjplimp@sandia.gov
//
//   Copyright (2003) Sandia Corporation.  Under the terms of Contract
//   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
//   certain rights in this software.  This software is distributed under 
//   the GNU General Public License.
//
//   See the README file in the top-level LAMMPS directory.
//------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Tommi Järvi (Fh-IWM), Lars Pastewka (KIT)
------------------------------------------------------------------------- */

#ifndef ANDERSON_MIXER_H
#define ANDERSON_MIXER_H

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "error.h"
#include "mpi.h"

namespace LAMMPS_NS {

/*
  Anderson mixer for solving charges
*/
class AndersonMixer {
 public:  
  AndersonMixer();
  AndersonMixer(int nhist);
  ~AndersonMixer();
  void mix(int iter, int n, gsl_vector* xi, const gsl_vector* yi, double beta, MPI_Comm &world, Error *&error);

 private:
  gsl_matrix* F_hist;  // residual history
  gsl_matrix* x_hist;  // vector history
  int M;               // length of history
  int n_hist;          // size of history arrays
};

}

#endif