///* ----------------------------------------------------------------------
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

#ifdef FIX_CLASS

FixStyle(bsct,FixBSCT)

#else

#ifndef LMP_FIX_BSCT_H
#define LMP_FIX_BSCT_H

#include "stdio.h"
#include "fix.h"
#include "update.h"
#include "pointers.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "fix_bsct_parameters_t.h"
#include "pair_coul_long_bsct.h"
#include "pair_coul_cut_bsct.h"
#include "pppm_bsct.h"


namespace LAMMPS_NS {


  /*
    Anderson mixer for solving charges
  */
  class Anderson_Mixer {

  public:  
    Anderson_Mixer();
    Anderson_Mixer(int nhist);
    ~Anderson_Mixer();
    void mix(int iter, int n, gsl_vector* xi, const gsl_vector* yi, double beta, MPI_Comm &world, Error *&error);

  private:
    gsl_matrix* F_hist;  // residual history
    gsl_matrix* x_hist;  // vector history
    int M;               // length of history
    int n_hist;          // size of history arrays
  };

  /*
    The BSCT fix class
  */
class FixBSCT : public Fix {
 public:
    FixBSCT(class LAMMPS *, int, char **);
    ~FixBSCT();
    int setmask();
    void setup_pre_force(int);
    void pre_force(int);
    double memory_usage();
    void setup(int);
    void init();
    void grow_arrays(int);
    void copy_arrays(int, int);
    void set_arrays(int);
    int pack_exchange(int, double *);
    int unpack_exchange(int, double *);
    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);
    double compute_scalar();

  private:
    // per-atom arrays
    double *phi;                // electrostatic potential of each atom

    // Internal Coulomb classes
    PPPMBSCT *pppm;
    PairCoulLongBSCT *pcl;
    PairCoulCutBSCT *pcc;

    // General settings
    int log;                    // 0 -> no logging
                                // 1 -> log to log.lammps
                                // 2 -> more log, incl. Newton's method
                                // 3 -> debug
                                // Note: log is set to zero on all but
                                //       one node so that only one
                                //       node outputs.
    int plog;                   // plog=log, and not set to zero on
				// other nodes.
    double qtot;                // total charge of atoms within this charge transfer fix
    double dGimax;              // convergence if |dG/dq_i| < dGimax
    int anderson_mode;          // 1 - Newton's iteration, 2 - original, 3 - first 2, then 1

    // Anderson mixing parameters
    int max_iter;               // max number of iterations
    int and_history;            // history parameter for Anderson mixer
    double beta;                // mixing parameter
    double beta_max;            // max for beta
    double beta_mul;            // increase factor

    // Other variables
    int nparams;                   // length of parameters vector
    FixBSCT_parameters_t *params;  // model parameters [0,N] with parameters in [1,N] corresp. to atom type

    double lambda;              // current value of Lagrange multiplier
    int ncq;                    // number of local charges within fix
				// on this timestep: determined every
				// time step and used for error
				// checking.
    int iter;                   // current iteration
    int ncall;                  // how manyeth call?

    // functions
    void update_charges_anderson();
    void adjust_total_charge();
    void gsl_absmax_sqnorm(gsl_vector*, double*, double *);
    double nextq(double q, double phi, double lambda, double X, double U, double V, double p);
    void FixBSCT_fdf(const gsl_vector *x, double *f, gsl_vector *g);
    void get_lambda();
    void update_kspace();
};

}

#endif
#endif
