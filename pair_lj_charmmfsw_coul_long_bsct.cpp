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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_charmmfsw_coul_long_bsct.h"
//0#include "pair_coul_long.h"
#include "pair_lj_charmmfsw_coul_long.h" //1
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

//0PairLJCharmmfswCoulLongBSCT::PairLJCharmmfswCoulLongBSCT(LAMMPS *lmp) : PairCoulLong(lmp) {}
PairLJCharmmfswCoulLongBSCT::PairLJCharmmfswCoulLongBSCT(LAMMPS *lmp) : PairLJCharmmfswCoulLong(lmp) {} //1

/* ---------------------------------------------------------------------- */

void PairLJCharmmfswCoulLongBSCT::reset_g_ewald()
{
 if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;
  if (ncoultablebits) init_tables(cut_coul,NULL);
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmfswCoulLongBSCT::compute_potential(double &ecoultot, double *phi)
{
  // Parts removed from original PairLJCharmmfswCoulLong::compute
  // commented by //--, modified parts suffixed with //**,
  // inserted parts marked by //++ or // bsct { ... // } bsct
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul; //** evdwl,evdwl12,evdwl6,fpair
  double fraction,table;
  double r,rinv,r2inv,rsq,factor_coul; //** r3inv,r6inv,forcelj,forcecoul,factor_lj
  double grij,expm2,prefactor,t,erfc;
  //-- double switch1;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double phii, phij;  //++
  double phiprefactor;

  ecoul = 0.0; //** evdwl = ecoul = 0.0;
  //-- ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
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

      if (rsq < cut_coulsq) { //** if (rsq < cut_bothsq)
        r2inv = 1.0/rsq;
        //-- if (rsq < cut_coulsq)
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * qtmp*q[j]/r;
          //-- forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          //-- if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

          // bsct {
          phiprefactor = qqrd2e/r;
          phii = phiprefactor*q[j]*erfc;
          phij = phiprefactor*qtmp*erfc;
          if (factor_coul < 1.0) {
            // Tommi: erfc not included here in pair_coul_long either
            phii -= (1.0-factor_coul)*phiprefactor*q[j];
            phij -= (1.0-factor_coul)*phiprefactor*qtmp;
          }
          // } bsct

        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];

          //-- forcecoul = qtmp*q[j] * table;
          //-- if (factor_coul < 1.0) {
          //--   table = ctable[itable] + fraction*dctable[itable];
          //--   prefactor = qtmp*q[j] * table;
          //--   forcecoul -= (1.0-factor_coul)*prefactor;
          //-- }

          // bsct {
          phiprefactor = table;
          phii = phiprefactor*q[j];
          phij = phiprefactor*qtmp;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            //-- set scale[itype][jtype] = 1
            //-- prefactor = scale[itype][jtype] * qtmp*q[j] * table;
            prefactor = qtmp*q[j] * table;
            phiprefactor = table;
            phii -= (1.0-factor_coul)*phiprefactor*q[j];
            phij -= (1.0-factor_coul)*phiprefactor*qtmp;
          }
          // } bsct
        }

        // accumulate phi
        phi[i] += phii; //++
        if (newton_pair || j < nlocal) phi[j] += phij; //++

        // calculate electrostatic energy
        if (!ncoultablebits || rsq <= tabinnersq) {
          ecoul = prefactor*erfc;
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        } else {
          table = etable[itable] + fraction*detable[itable];
          ecoul = qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ptable[itable] + fraction*dptable[itable];
            prefactor = qtmp*q[j] * table;
            ecoul -= (1.0-factor_coul)*prefactor;
          }
        }

        // bsct {
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
