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

#include "anderson_mixer.h"

using namespace LAMMPS_NS;

AndersonMixer::AndersonMixer() {
  M = 3;
  n_hist = 0;
}

AndersonMixer::AndersonMixer(int nhist) {
  M = nhist;
  n_hist = 0;
}

AndersonMixer::~AndersonMixer() {
  if(n_hist > 0) {
    gsl_matrix_free(F_hist);
    gsl_matrix_free(x_hist);
  }
}


/* ----------------------------------------------------------------------
   Anderson mixer

   iter              input: current iteration (first time should be 1)
   n                 input: dimension of vectors
   xi                input: previous iteration, output: new vector
   yi                input: new trial charge vector
   beta              input: mixing parameter
   ---------------------------------------------------------------------- */

void AndersonMixer::mix(int iter, int n, gsl_vector* xi, const gsl_vector* yi, double beta, MPI_Comm &world, Error *&error) {

  // --- variables

  int i, j;                   // loops
  int nm;                     // order of matrix equation
  int me;                     // my MPI rank
  int status;                 // GSL status
  double *buf;                // MPI buffer for matrix equation
  double *buf2;               // MPI buffer for matrix equation
  int bufsize;                // length of MPI buffer for matrix equation
  int nbuf;                   // looping over buffer elements
  double *coeff;              // coefficients for new vector (size M)
  int singular;               // is matrix A singular? (0 = not)
  int temp;


  // --- check

  if(M<0) {
    error->all(FLERR, "AndersonMixer::mix: Requested mixing history < 0!");
  }


  // --- allocate F_hist, x_hist if needed

  if(n_hist != n) {
    if(iter != 1) {
      error->all(FLERR, "AndersonMixer::mix: Bug: Would have to resize history arrays after first iteration (might be ok in parallel)!");
    }
    if(n_hist > 0) {
      gsl_matrix_free(F_hist);
      gsl_matrix_free(x_hist);
    }
    if(n > 0) {
      F_hist = gsl_matrix_alloc(M+1, n);
      x_hist = gsl_matrix_alloc(M+1, n);
      n_hist = n;
    }
    else {
      F_hist = gsl_matrix_alloc(M+1, 1);
      x_hist = gsl_matrix_alloc(M+1, 1);
      n_hist = 1;
    }
  }


  // --- setup

  // number of history that should be stored for us
  if(iter-1 < M) nm = iter-1;
  else nm = M;

  // input and residual
  for(i=0;i<n;i++) {
    gsl_matrix_set(F_hist, 0, i, gsl_vector_get(yi,i) - gsl_vector_get(xi,i));  // current residual
    gsl_matrix_set(x_hist, 0, i, gsl_vector_get(xi,i));  // previous input vector
  }

  // IF NO HISTORY, RETURN WITH SIMPLE MIXING
  if(nm==0) {

    // shift history
    temp = nm + 1;
    if(temp > M) temp = M;
    for(j=temp;j>0;j--) {
      for(i=0;i<n;i++) {
        gsl_matrix_set(F_hist, j, i, gsl_matrix_get(F_hist, j-1, i));
        gsl_matrix_set(x_hist, j, i, gsl_matrix_get(x_hist, j-1, i));
      }
    }

    for(i=0;i<n;i++) {
      gsl_vector_set(xi, i, (1.0-beta)*gsl_vector_get(xi, i) + beta*gsl_vector_get(yi, i));
    }
    return;
  }
  else if(nm < 0) {
    error->all(FLERR, "AndersonMixer::mix(): Bug!");
  }

  // allocate MPI buffers and coefficients
  bufsize = (nm)*(nm) + (nm);
  buf = new double[bufsize];
  buf2 = new double[bufsize];
  coeff = new double[nm];

  // setup matrix equation
  if(n > 0) {
    gsl_vector_view c1, c2, c3;
    gsl_vector *t1 = gsl_vector_alloc(n);
    gsl_vector *t2 = gsl_vector_alloc(n);
    nbuf = 0;
    c1 = gsl_matrix_row(F_hist, 0);  // col. 0
    for(i=1;i<=nm;i++) {
      c2 = gsl_matrix_row(F_hist, i);  // col. i
      gsl_vector_memcpy(t1, &(c1.vector));
      gsl_vector_sub(t1, &(c2.vector));  // t1 = col. 0 - col. i
      gsl_blas_ddot(t1, &(c1.vector), &buf[nbuf]);  // -> b[i] = (c0 - ci) .* c0
      nbuf++;
      for(j=1;j<=nm;j++) {
        c3 = gsl_matrix_row(F_hist, j);  // col. j
        gsl_vector_memcpy(t2, &(c1.vector));
        gsl_vector_sub(t2, &(c3.vector));  // t2 = col. 0 - col. j
        gsl_blas_ddot(t1, t2, &buf[nbuf]);  // -> A[i][j] = (c0 - ci) .* (c0 - cj)
        nbuf++;
      }
    }
    if(nbuf!=bufsize) {
      error->all(FLERR, "Mixer bug!");
    }
    gsl_vector_free(t1);
    gsl_vector_free(t2);
  }
  else {
    for(i=0;i<bufsize;i++) {
      buf[i] = 0.0;
    }
  }

  // MPI: Sum buffer on root node buf -> buf2
  MPI_Reduce(buf, buf2, bufsize, MPI_DOUBLE, MPI_SUM, 0, world);

  // --- solve matrix equation to get coefficients for constructing new vector
  //             !----------------------------
  //             ! solve A*z=b -> z-->b
  //             ! (eq.(4.3) in Eyert)
  //             !----------------------------

  // MPI root node
  MPI_Comm_rank(world,&me);
  if(me==0) {
    gsl_matrix *A = gsl_matrix_alloc(nm, nm);
    gsl_vector *b = gsl_vector_alloc(nm);
    // construct A and b
    nbuf = 0;
    for(i=0;i<nm;i++) {
      gsl_vector_set(b, i, buf2[nbuf]);
      nbuf++;
      for(j=0;j<nm;j++) {
        gsl_matrix_set(A, j, i, buf2[nbuf]);  // j,i instead of i,j because C array convention?
        nbuf++;
      }
    }
    // solve A*x = b, put x to b
    {
      gsl_permutation *p = gsl_permutation_alloc(nm);
      gsl_vector *x = gsl_vector_alloc(nm);
      int signum;

      // LU decomposition
      status = gsl_linalg_LU_decomp(A, p, &signum);
      if(status) {
        printf("GSL error: %s\n", gsl_strerror(status));
        printf("    error code: %d\n", status);
        error->all(FLERR, "AndersonMixer::mix(): Problem with GSL.");
      }

      // check if singular => if, use simple mixing
      singular = 0;
      for(i=0;i<nm;i++) {
        if(gsl_matrix_get(A, i, i) == 0.0) {
          singular = 1;
          break;
        }
      }

      // solve matrix equation
      if(singular==0) {
        status = gsl_linalg_LU_solve(A, p, b, x);
        if(status) {
          printf("GSL error: %s\n", gsl_strerror(status));
          printf("    error code: %d\n", status);
          error->all(FLERR, "AndersonMixer::mix(): Problem with GSL.");
        }
        for(i=0;i<nm;i++) {
          coeff[i] = gsl_vector_get(x, i);
        }
      }
      gsl_vector_free(x);
      gsl_permutation_free(p);
    }
    gsl_matrix_free(A);
    gsl_vector_free(b);
  }  // end of MPI root node

  // MPI: Broadcast info if A was singular
  MPI_Bcast(&singular, 1, MPI_INT, 0, world);

  // Calculate new charges
  if(singular==0) {
    double xb[n], Fb[n];

    // MPI: Broadcast optimal linear combination
    MPI_Bcast(coeff, nm, MPI_DOUBLE, 0, world);

    for(i=0;i<n;i++) {
      xb[i] = gsl_vector_get(xi,i);
      Fb[i] = gsl_matrix_get(F_hist, 0, i);
    }
    for(j=1;j<=nm;j++) {
      for(i=0;i<n;i++) {
        xb[i] += coeff[j-1]*(gsl_matrix_get(x_hist, j, i) - gsl_vector_get(xi,i));
        Fb[i] += coeff[j-1]*(gsl_matrix_get(F_hist, j, i) - gsl_matrix_get(F_hist, 0, i));
      }
    }
    // new vector
    for(i=0;i<n;i++) {
      gsl_vector_set(xi, i, xb[i] + beta*Fb[i]);
    }
  }
  else {
    // use simple mixing
    for(i=0;i<n;i++) {
      gsl_vector_set(xi, i, (1.0-beta)*gsl_matrix_get(x_hist, 0, i) + beta*gsl_vector_get(yi,i));
    }
  }


  // --- shift history (0->1 etc.)

  if(iter < M) nm = iter;
  else nm = M;

  for(j=nm;j>0;j--) {
    for(i=0;i<n;i++) {
      gsl_matrix_set(F_hist, j, i, gsl_matrix_get(F_hist, j-1, i));
      gsl_matrix_set(x_hist, j, i, gsl_matrix_get(x_hist, j-1, i));
    }
  }


  // --- free memory

  delete [] buf;
  delete [] buf2;
  delete [] coeff;
}
