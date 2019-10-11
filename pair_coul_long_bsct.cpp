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
#include "pair_coul_long_bsct.h"
#include "pair_coul_long.h"
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

PairCoulLongBSCT::PairCoulLongBSCT(LAMMPS *lmp) : PairCoulLong(lmp) {}

/* ---------------------------------------------------------------------- */

void PairCoulLongBSCT::reset_g_ewald()
{
 if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;
  if (ncoultablebits) init_tables(cut_coul,NULL);
}

/* ---------------------------------------------------------------------- */

void PairCoulLongBSCT::phi(double &ecoultot, double *phi)
{
  int i,j,ii,jj,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;
  double phii, phij;  // Tommi

  ecoul = 0.0;
  /*
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  */

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * scale[itype][jtype] * qtmp*q[j]/r;
          //forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          //if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

          // Tommi
          {
            double phiprefactor;
            phiprefactor = qqrd2e/r;
            phii = phiprefactor*q[j]*erfc;
            phij = phiprefactor*qtmp*erfc;
            if (factor_coul < 1.0) {
              // Tommi: erfc not included here in pair_coul_long either
              phii -= (1.0-factor_coul)*phiprefactor*q[j];
              phij -= (1.0-factor_coul)*phiprefactor*qtmp;
            }
          }
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction*dftable[itable];
          /*
          forcecoul = scale[itype][jtype] * qtmp*q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction*dctable[itable];
            prefactor = scale[itype][jtype] * qtmp*q[j] * table;
            forcecoul -= (1.0-factor_coul)*prefactor;
          }
          */

          // Tommi
          {
            double phiprefactor;
            phiprefactor = table;
            phii = phiprefactor*q[j];
            phij = phiprefactor*qtmp;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = scale[itype][jtype] * qtmp*q[j] * table;
              phiprefactor = table;
              phii -= (1.0-factor_coul)*phiprefactor*q[j];
              phij -= (1.0-factor_coul)*phiprefactor*qtmp;
            }
          }
        }

        // Tommi

        // accumulate phi
        phi[i] += phii;
        if (newton_pair || j < nlocal) phi[j] += phij;

        // calculate electrostatic energy
        if (!ncoultablebits || rsq <= tabinnersq)
          ecoul = prefactor*erfc;
        else {
          table = etable[itable] + fraction*detable[itable];
          ecoul = scale[itype][jtype] * qtmp*q[j] * table;
        }
        if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;

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

  //if (vflag_fdotr) virial_fdotr_compute();
}
