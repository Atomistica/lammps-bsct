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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_coul_cut.h"
#include "pair_coul_cut_bsct.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutBSCT::PairCoulCutBSCT(LAMMPS *lmp) : PairCoulCut(lmp) {}

/* ---------------------------------------------------------------------- */

void PairCoulCutBSCT::phi(double &ecoultot, double *phi)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,r2inv,rinv,forcecoul,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double phiprefactor, phii, phij;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
/*
      if (j < nall) factor_coul = 1.0;
      else {
        factor_coul = special_coul[j/nall];
        j %= nall;
      }
*/
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);
        forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
        fpair = factor_coul*forcecoul * r2inv;

        // Tommi
        phiprefactor = factor_coul * qqrd2e * scale[itype][jtype] * rinv;
        phii = phiprefactor*q[j];
        phij = phiprefactor*qtmp;

        ecoul = phiprefactor * qtmp*q[j];

        // accumulate phi
        phi[i] += phii;
        if (newton_pair || j < nlocal) phi[j] += phij;

        // accumulate total electrostatic energy
        //   - adapted from Pair::ev_tally()
        if (newton_pair) {
          ecoultot += ecoul;
        }
        else {
          if (i < nlocal) {
            ecoultot += ecoul*0.5;
          }
          if (j < nlocal) {
            ecoultot += ecoul*0.5;
          }
        }
      }

    }
  }

}
