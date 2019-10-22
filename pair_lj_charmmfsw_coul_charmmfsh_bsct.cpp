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
#include "pair_lj_charmmfsw_coul_charmmfsh_bsct.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmfswCoulCharmmfshBSCT::PairLJCharmmfswCoulCharmmfshBSCT(LAMMPS *lmp) : PairLJCharmmfswCoulCharmmfsh(lmp) {}

/* ---------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfshBSCT::compute_potential(double &ecoultot, double *phi)
{
  // Parts removed from original PairLJCharmmfswCoulCharmmfsh::compute commented
  // by //--, parts removed from PairCoulCut::compute commented with //__
  // modified parts suffixed by //**, inserted parts by //++
  // inserted multi-line parts marked by // bsct { .. }
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul; //** evdwl,evdwl12,evdwl6,fpair
  double r,rinv,rsq,r2inv,factor_coul; //** r3invr6inv,forcelj,forcecoul,factor_lj
  //-- double switch1;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double phiprefactor, phii, phij; //++

  ecoul = 0.0; //** evdwl = ecoul = 0.0;
  //-- ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  // int nall = nlocal + atom->nghost; // where from, why necessary?
  double *special_coul = force->special_coul;
  //-- double *special_lj = force->special_lj;
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
      //-- factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      //__ jtype = type[j];

      //__ if (rsq < cutsq[itype][jtype]) {
      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        rinv = sqrt(r2inv); //++

        //set scale[itype][jtype] to "1"
        //__ forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;
        //__ forcecoul = qqrd2e * qtmp*q[j]*rinv;
        //-- forcecoul = qqrd2e * qtmp*q[j]*
        //--  (sqrt(r2inv) - r*cut_coulinv*cut_coulinv); // force shifting
        //-- fpair = factor_coul*forcecoul * r2inv;

        // Do we have to adapt following potential computation to be consistent
        // with CHARMM standard force shifting, i.e. integrate force computation
        // above? Not yet done in the following:
        phiprefactor = factor_coul * qqrd2e * rinv;
        phii = phiprefactor*q[j];
        phij = phiprefactor*qtmp; // qtmp is q[i]

        //__ ecoul = phiprefactor * qtmp*q[j];
        ecoul = qqrd2e * qtmp*q[j]*
          (sqrt(r2inv) + cut_coulinv*cut_coulinv*r - 2.0*cut_coulinv);
        ecoul *= factor_coul; // (1/r + 1/r_c^2 - 2/r_c) = (r-r_c)^2/(r r_c^2)

        // bsct {
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
        // } bsct
      }

    }
  }

}
