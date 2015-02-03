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
   Contributing authors: Tommi Järvi (Fh-IWM), Lars Pastewka (KIT)
------------------------------------------------------------------------- */

/*
1. Check nbuf
2. There is an if differentiation_flag != 1 in poisson_peratom. Check if
   we can only compute potential if differentiation_flag != 1.
*/

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pppm_bsct.h"
#include "atom.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

PPPMBSCT::PPPMBSCT(LAMMPS *lmp, int narg, char **arg) : PPPM(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMBSCT::~PPPMBSCT()
{
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver
------------------------------------------------------------------------- */

void PPPMBSCT::poisson_peratom()
{
  PPPM::poisson_peratom();

  int i,j,k;

  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        phi_brick[k][j][i] = u_brick[k][j][i];
      }
}

/* ----------------------------------------------------------------------
   Return electrostatic potential (Compute has to have been run
   prior to this as the actual computation is done there. This only
   extracts the potential from the internal arrays.)
------------------------------------------------------------------------- */

void PPPMBSCT::phi(FFT_SCALAR *phi)
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;

  int nlocal = atom->nlocal;

  const double qscale = qqrd2e * scale;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          phi[i] += qscale*x0*u_brick[mz][my][mx];
        }
      }
    }
  }

  for (i = 0; i < nlocal; i++) {
    phi[i] -= 2*qscale*(g_ewald*q[i]/MY_PIS + MY_PI2*qsum / (g_ewald*g_ewald*volume));
  }

  // slab correction
  if(slabflag) phi_slabcorr(phi);
}


/* ----------------------------------------------------------------------
    Slab correction for phi. Not tested too well.
 ------------------------------------------------------------------------- */

void PPPMBSCT::phi_slabcorr(double *phi)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections

  double ffact = -4.0*MY_PI*dipole_all/volume;
  for (int i = 0; i < nlocal; i++) {
    phi[i] += qqrd2e*(-ffact)*atom->x[i][2];
  }
}

/* ----------------------------------------------------------------------
   allocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMBSCT::allocate_peratom()
{
  PPPM::allocate_peratom();

  memory->create3d_offset(phi_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                          nxlo_out,nxhi_out,"pppm:phi_brick");
}

/* ----------------------------------------------------------------------
   deallocate per-atom memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMBSCT::deallocate_peratom()
{
  PPPM::deallocate_peratom();

  memory->destroy3d_offset(phi_brick,nzlo_out,nylo_out,nxlo_out);
}
