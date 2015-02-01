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

#include "stdlib.h"
#include "string.h"
#include "fix_bsct.h"
#include "atom.h"
#include "input.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math.h"
#include "comm.h"
#include "pair.h"
#include "universe.h"
#include "mpi.h"
#include "anderson_mixer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contructor and parameter setting
------------------------------------------------------------------------- */

FixBSCT::FixBSCT(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // --- set Fix flags etc.
  //     - for flags, see fix.h
  //     - for callback, see atom.h

  restart_peratom = 0;                // restart not supported yet
  create_attribute = 1;               // per-atom data exists
  atom->add_callback(0);
  peratom_flag = 1;                   // storing per-atom data
  size_peratom_cols = 0;              // simple per-atom vector
  peratom_freq = 1;                   // per-atom data available every step

  scalar_flag = 1;                    // compute_scalar() exists and calculates BSCT energy
  extscalar = 1;                      // scalar is extensive (not intensive)
  global_freq = 1;                    // scalar available every timestep

  // allocate internal per-atom arrays
  // (also setup() allocates these, but needed here too)
  phi = NULL;
  grow_arrays(atom->nmax);


  // --- set default parameters

  nevery = 1;                 // update charges every nevery steps (inherited from Fix)
  qtot = 0.0;                 // Total charge of atoms within fix
  anderson_mode = 2;          // see fix_bsct.h
  and_history = 3;            // see fix_bsct.h
  log = 0;                    // see fix_bsct.h
  plog = 0;                   // see fix_bsct.h
  nparams = 0;                // length of parameters vector
  dGimax = 0.05;              // Anderson & C-P: convergence if |dG/dq_i| < dGimax
  max_iter = 1000;            // And: max number of iterations
  beta = 0.1;                 // And: mixing parameter
  beta_max = 0.1;             // And: max for beta
  beta_mul = 1.0;             // And: increase factor

  bool slab = false;          // slab calculation requested?
  std::string slab_volfactor;      // slab_volfactor for pppm
  std::string pppmprec = "1.0e-4"; // Internal PPPM precision
  std::string coulcut = "10.0";    // Internal Coulomb cutoff

  // --- allocate and reset parameter array

  // allocate
  nparams = atom->ntypes+1;
  if(nparams > 0) {
    params = new FixBSCT_parameters_t[nparams];
  }
  else {
    params = NULL;
    nparams = 0;
  }

  // reset
  for(int i=0;i<nparams;i++) {
    params[i].set = false;
    params[i].X = 0.0;
    params[i].U = 0.0;
    params[i].V = 0.0;
    params[i].p = 0.0;
  }


  // --- determine parameters from input

  {
    int iarg = 3;  // where we are on the input line
    int np = 0;    // number of atom type params read in

    while(iarg < narg) {

      if(strcmp(arg[iarg],"nevery")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        nevery = atoi(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"log")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        log = atoi(arg[iarg + 1]);
        plog = log;
        iarg += 2;

        // log only for one processor
        if(comm->me!=0) log = 0;
      }

      else if(strcmp(arg[iarg],"qtot")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        qtot = atof(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"mode")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        anderson_mode = atoi(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"history")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        and_history = atoi(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"beta")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        beta = atof(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"beta_max")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        beta_max = atof(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"beta_mul")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        beta_mul = atof(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"prec")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        dGimax = atof(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"max_iter")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        max_iter = atoi(arg[iarg + 1]);
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"type")==0) {
        int nt;

        if(narg < iarg + 6) error->all(FLERR, "Illegal fix bsct command");
        iarg++;

        np++;

        // atom type
        nt = atoi(arg[iarg]);
        iarg++;

        if(nt < 0 || nt > atom->ntypes) {
          error->all(FLERR, "Illegal fix bsct command: Invalid atom type.");
        }
        if(params[nt].set) {
          error->all(FLERR, "Illegal fix bsct command: Trying to add two sets of values to an atom type.");
        }

        params[nt].set = true;
        params[nt].X = atof(arg[iarg]);
        iarg++;
        params[nt].U = atof(arg[iarg]);
        iarg++;
        params[nt].V = atof(arg[iarg]);
        iarg++;
        params[nt].p = atof(arg[iarg]);
        iarg++;
      }

      else if(strcmp(arg[iarg],"slab")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        slab_volfactor = arg[iarg + 1];
        slab = true;
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"pppmprec")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        pppmprec = arg[iarg + 1];
        iarg += 2;
      }

      else if(strcmp(arg[iarg],"coulcut")==0) {
        if(narg < iarg + 2) error->all(FLERR, "Illegal fix bsct command");
        coulcut = arg[iarg + 1];
        iarg += 2;
      }

      else {
        error->all(FLERR, "Illegal fix bsct command");
      }
    }

    // checks
    if(np < 1) {
      error->all(FLERR, "Illegal fix bsct command: No charge transfer parameters given.");
    }

    // check some defaults
    if(beta_max < beta) beta_max = beta;
  }

  // get Coulomb classes from forces object

  if (force->kspace != NULL) {
    pppm = dynamic_cast<PPPMBSCT*>(force->kspace);
    if (!pppm) {
      error->all(FLERR,"Please use kspace_style pppm/bsct.");
    }
    pcl = dynamic_cast<PairCoulLongBSCT*>(force->pair);
    if (!pcl) {
      error->all(FLERR,"Please use pair_style coul/long/bsct.");
    }
  }
  else {
    pcc = dynamic_cast<PairCoulCutBSCT*>(force->pair);
    if (!pcl && !pcc) {
      error->all(FLERR,"Unsupported pair style. Please use pair_style "
                 "coul/cut/bsct or kspace_style pppm/bsct and pair_style "
                 "coul/long/bsct.");
    }
  }


  // --- report

  if(comm->me==0) {
    fprintf(logfile, "Band structure charge transfer fix parameters:\n");
    fprintf(logfile, "  Update charges every %d steps\n", nevery);
    fprintf(logfile, "  Logging level (0-3)  = %d\n", log);
    fprintf(logfile, "  Total charge = %f\n", qtot);
    fprintf(logfile, "  Convergence criterion for |dG/dq_i| = %f\n", dGimax);
    if(pppm != NULL) fprintf(logfile, "  PPPM accurafy = %s\n", pppmprec.c_str());
    if(slab) fprintf(logfile, "  PPPM slab correction width = %s\n", slab_volfactor.c_str());
    fprintf(logfile, "  Short-range Coulomb cutoff = %s\n", coulcut.c_str());
    fprintf(logfile, "  Convergence criterion for |dG/dq_i| = %f\n", dGimax);
    fprintf(logfile, "  Anderson mixing options:\n");
    fprintf(logfile, "    Mode = %d\n", anderson_mode);
    fprintf(logfile, "    Maximum iterations = %d\n", max_iter);
    fprintf(logfile, "    Mixing parameters:\n");
    fprintf(logfile, "       beta = %f\n", beta);
    fprintf(logfile, "       beta_max = %f\n", beta_max);
    fprintf(logfile, "       beta_mul = %f\n", beta_mul);
    fprintf(logfile, "       history = %d\n", and_history);

    fprintf(logfile, "  Atom type parameters:\n");
    fprintf(logfile,"       type               X               U               V               p\n");
    for(int i=0;i<nparams;i++) {
      if(params[i].set) {
        fprintf(logfile,"      %5d %15.5f %15.5f %15.5f %15.5f\n", i, params[i].X, params[i].U, params[i].V, params[i].p);
      }
    }
    fprintf(logfile, "\n");
  }

}


/* ----------------------------------------------------------------------
   Destructor
------------------------------------------------------------------------- */

FixBSCT::~FixBSCT()
{
  if(nparams > 0) {
    delete [] params;
  }
  nparams = 0;

  if(phi != NULL) memory->destroy(phi);

  atom->delete_callback(id, 0);
}


/* ----------------------------------------------------------------------
   Masks
------------------------------------------------------------------------- */

int FixBSCT::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}


/* ----------------------------------------------------------------------
   Init: Only initialize internal Coulomb classes
------------------------------------------------------------------------- */

void FixBSCT::init()
{
}


/* ----------------------------------------------------------------------
   Setup
------------------------------------------------------------------------- */

void FixBSCT::setup(int vflag)
{
  // setup for internal Coulomb classes

  // make sure arrays are resized if necessary
  grow_arrays(atom->nmax);

  // total charge
  adjust_total_charge();

  // mark first call and reset
  ncall = 0;
  lambda = 0.0;
}


/* ----------------------------------------------------------------------
   Pre-force
------------------------------------------------------------------------- */

void FixBSCT::pre_force(int vflag)
{
  // keep track of call
  ncall++;  // on first call, this becomes 1

  // nevery
  if((ncall-1)%nevery!=0) return;

  // update charges
  update_charges_anderson();

  // update kspace if needed
  update_kspace();
}


/* ----------------------------------------------------------------------
   Called to update k-space parameters after converging charges
------------------------------------------------------------------------- */

void FixBSCT::update_kspace() {

  if(pppm != NULL) {
    force->kspace->init();
    force->kspace->setup();

    force->pair->init_style();
  }

}


/* ----------------------------------------------------------------------
  Calculate energy from charge transfer (band structure, Hubbard-U and
  electronegativity/bias contributions.
------------------------------------------------------------------------- */

double FixBSCT::compute_scalar()
{

  int i;
  double esum = 0.0;
  double energy;

  if(nparams < 1) {
    error->all(FLERR, "fix bsct: Bug: Parameter data not allocated!");
  }

  for(i=0; i<atom->nlocal; i++) {
    if(atom->mask[i] & groupbit) {
      // params
      int itype = atom->type[i];
      if(itype >= nparams) error->all(FLERR, "fix bsct: Bug: Atom type number too big!");
      if(!params[itype].set) error->all(FLERR, "fix bsct: Parameters missing for atom type included in fix group!");
      double X = params[itype].X;
      double U = params[itype].U;
      double V = params[itype].V;
      double p = params[itype].p;

      esum += (V/2)*pow(fabs(atom->q[i]),p) + (U/2)*atom->q[i]*atom->q[i] + atom->q[i]*X;
    }
  }
  MPI_Allreduce(&esum, &energy, 1, MPI_DOUBLE, MPI_SUM, world);

  return energy;
}


/* ----------------------------------------------------------------------
   Calculate next charges from iterative formula: q_{i+1} = nextq(q_i)
------------------------------------------------------------------------- */

double FixBSCT::nextq(double q, double phi, double lambda, double X, double U, double V, double p) {

  double s;
  double temp;
  double newq;

  // --------------------------
  // --- Newton's iteration ---
  // --------------------------
  if(anderson_mode==1 || (anderson_mode==3 && ncall > 2)) {
    temp = phi + U*q + X + lambda;
    if(temp < 0.0) s = -1.0;
    else s = 1.0;

    temp = fabs(temp)*2.0/(p*V);

    newq = -s * pow(temp, 1.0/(p-1.0));
  }
  // ----------------
  // --- Original ---
  // ----------------
  else if(anderson_mode==2 || (anderson_mode==3 && ncall <= 2)) {
    if(q < 0.0) s = -1.0;
    else s = 1.0;

    newq = -(1.0/U)*(s*(p*V/2)*pow(fabs(q),p-1) + phi + X + lambda);
  }
  else {
    error->all(FLERR, "fix bsct: Unknown Anderson mixer mode.");
  }

  return newq;
}


/* ----------------------------------------------------------------------
   Calculate Lagrange multiplier, lambda
------------------------------------------------------------------------- */

void FixBSCT::get_lambda() {

  int i;                       // loops

  double s;                    // temp for sign
  double flam, fplam;          // f(lambda) and f'(lambda) for Newton's method and determining lambda
  double gflam, gfplam;        // global variable for MPI
  int newtonit;                // counter for Newton's method iterations
  double sum1, sum2;           // sums
  double gsum1, gsum2;         // globals for sums

  double qconv = 1e-6;         // convergence criterion for lambda
  int maxnewton = 1000;        // maximum iterations for Newton's method: hardcoded, just to prevent infinite loops

  // --------------------------
  // --- Newton's iteration ---
  // --------------------------
  if(anderson_mode==1 || (anderson_mode==3 && ncall > 2)) {
    //     f(lambda)  = Q + \sum 2/(pV) |phi + Uq + chi + lambda|^(1/(p-1)) sgn(phi + Uq + chi + lambda)
    //     f'(lambda) = \sum ( pV/2 (p-1) )^-1 |phi + Uq + chi + lambda|^((2-p)/(p-1))

    // Iteration: Newton's method for finding lambda
    if(log >= 2) fprintf(logfile,"                       Newton's method:    lambda       f(lambda)      f'(lambda)\n");
    newtonit = 0;
    while(newtonit <= maxnewton) {
      // compute sums
      flam = 0.0;
      fplam = 0.0;
      for(i=0; i<atom->nlocal; i++) {
        if(atom->mask[i] & groupbit) {
          double temp;
          // params
          int itype = atom->type[i];
          double X = params[itype].X;
          double U = params[itype].U;
          double V = params[itype].V;
          double p = params[itype].p;

          // f(lambda)
          flam += nextq(atom->q[i], this->phi[i], lambda, X, U, V, p);

          // f'(lambda)
          temp = this->phi[i] + U*atom->q[i] + X + lambda;
          if(temp > 0.0) s = 1.0;
          else s = -1.0;
          temp = fabs(temp);
          fplam += pow(2.0/(p*V), 1/(p-1.0)) * pow(temp, (2.0-p)/(p-1.0)) / (p-1.0);
        }
      }

      // compute f(lambda), f'(lambda)
      MPI_Allreduce(&flam, &gflam, 1, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(&fplam, &gfplam, 1, MPI_DOUBLE, MPI_SUM, world);
      gflam = qtot - gflam;

      // log
      if(log >= 2) fprintf(logfile,"                            %5d %15.5f %15.5f %15.5f\n", newtonit, lambda, gflam, gfplam);

      // check for convergence
      if(fabs(gflam) < qconv) break;

      // compute new lambda
      lambda = lambda - gflam/gfplam;
      newtonit++;

      // check that we're not in an infinite loop
      if(newtonit == maxnewton) {
        char buf[100];
        sprintf(buf, "fix bsct: Newton's iteration seems not to be converging after %d steps. f(lambda) = %e", maxnewton, gflam);
        error->all(FLERR, buf);
      }
    }  // end of Newton's method

    // log final result
    if(log >= 2) fprintf(logfile,"                            %5d %15.5f %15.5f %15.5f (converged)\n", newtonit, lambda, gflam, gfplam);

  }
  // ----------------
  // --- Original ---
  // ----------------
  else if(anderson_mode==2 || (anderson_mode==3 && ncall <= 2)) {

    sum1 = 0.0;
    sum2 = 0.0;
    for(i=0; i<atom->nlocal; i++) {
      if(atom->mask[i] & groupbit) {
        // params
        int itype = atom->type[i];
        double X = params[itype].X;
        double U = params[itype].U;
        double V = params[itype].V;
        double p = params[itype].p;

        if(atom->q[i] < 0.0) s = -1.0;
        else s = 1.0;

        // sum sum1, sum2 for lambda
        sum1 += 1/U;  // lambda
        sum2 += (s*(p*V/2.0)*pow(fabs(atom->q[i]),p-1.0) + this->phi[i] + X)/U;  // lambda
      }
    }

    // get global sums and max_dE
    MPI_Allreduce(&sum1, &gsum1, 1, MPI_DOUBLE, MPI_SUM, world);      // lambda
    MPI_Allreduce(&sum2, &gsum2, 1, MPI_DOUBLE, MPI_SUM, world);      // lambda

    // get the Lagrange multiplier, lambda
    lambda  = -(qtot + gsum2)/gsum1;

  }
  else {
    error->all(FLERR, "fix bsct: Unknown Anderson mixer mode.");
  }

}


/* ----------------------------------------------------------------------
   Calculate G, dG/dq_i and lambda
------------------------------------------------------------------------- */

void FixBSCT::FixBSCT_fdf(const gsl_vector *x, double *f, gsl_vector *g) {

  // --- variables
  int i;                       // loops

  double esum;                 // energy sum
  double qsum;                 // charge sum
  double qtsum;                // charge sum, also charges outside fix, for error checking
  size_t nand;                 // index for running over x and g
  double ecoul;                // Coulomb energy

  int initfrec = 10;           // frequency of calling PPPM::init()
  double qconv = 1e-6;         // warning threshold for non-zero charge (charge per atom)

  // ---
  // --- prepare
  // ---

  // unpack updated charges to q
  nand = 0;
  for(i=0; i<atom->nlocal; i++) {
    if(atom->mask[i] & groupbit) {
      atom->q[i] = gsl_vector_get(x, nand);
      nand++;
    }
  }

  // check
  if(nand != ncq) {
    error->all(FLERR, "fix bsct: Bug: Number of charges doesn't match!");
  }

  // MPI: Pass charges to all processors
  comm->forward_comm_fix(this);

  // ---
  // --- calculate electrostatic potential
  // ---

  // reset energy, phi
  ecoul = 0.0;
  for(i=0;i<atom->nlocal + atom->nghost;i++) {
    this->phi[i] = 0.0;
  }

  // calculate new potential
  if(pppm != NULL) {
    if(pcl==NULL) {
      error->all(FLERR, "fix bsct: Bug: Internal coul/long/bsct not setup correctly.");
    }

    // The frequency here is pretty arbitrary and the whole thing
    // could be cleanup up: Anyway, need to update g_ewald and PPPM
    // settings when charges change.
    if(iter == -1 || iter%initfrec==0) {
      pppm->init();
      pppm->setup();
    }
    else {
      pppm->qsum_update_flag = 1;
    }
    pcl->set_g_ewald(pppm->g_ewald); 
    pppm->compute(2,0);                    // energy to pppm->energy, including slab correction part
    pppm->phi(this->phi);                  // potential contribution to phi (including slab correction)
    pcl->phi(ecoul, this->phi);            // short range Coulombics
  }
  else {
    if(pcc==NULL) {
      error->all(FLERR, "fix bsct: Bug: Internal coul/cut/bsct not setup correctly.");
    }

    pcc->phi(ecoul, this->phi);            // short range Coulombics
  }

  // Get potential (phi) from all nodes, sum short range Coulomb part
  // and get total Coulomb energy.
  comm->reverse_comm_fix(this);
  {
    double temp;
    MPI_Allreduce(&ecoul, &temp, 1, MPI_DOUBLE, MPI_SUM, world);
    ecoul = temp;
  }
  if(pppm != NULL) {
    ecoul += pppm->energy;
  }

  // debug: Output Coulomb energy
  if(plog >= 3) {
    if(log >= 3) {
      if(pppm != NULL) {
        fprintf(logfile, "FixBSCT_fdf: short- and long-range Coulomb energies: %f %f\n", ecoul - pppm->energy, pppm->energy);
      }
      else {
        fprintf(logfile, "FixBSCT_fdf: Coulomb energy: %f\n", ecoul);
      }
    }

    // get energy from sum over q*phi
    double sum = 0.0;
    double gsum;
    for(i=0; i<atom->nlocal; i++) {
      sum += atom->q[i]*this->phi[i];
    }
    MPI_Allreduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, world);
    if(log >=3) fprintf(logfile, "FixBSCT_fdf: Coulomb energy from sum of q*phi: %f\n", 0.5*gsum);
  }

  // ---
  // --- calculate G and dG/dq_i without lambda
  // ---
  esum = 0.0;
  qsum = 0.0;
  qtsum = 0.0;
  nand = 0;
  for(i=0; i<atom->nlocal; i++) {
    // sum global q for checking that it's zero
    qtsum += atom->q[i];

    if(atom->mask[i] & groupbit) {
      // params
      int itype = atom->type[i];
      double X = params[itype].X;
      double U = params[itype].U;
      double V = params[itype].V;
      double p = params[itype].p;
      double s;  // temp for sign

      // sum q within fix
      qsum += atom->q[i];

      // gradient
      if(atom->q[i] < 0.0) s = -1.0;
      else s = 1.0;
      gsl_vector_set(g, nand, s*(p*V/2)*pow(fabs(atom->q[i]),p-1) + this->phi[i] + U*atom->q[i] + X);  // dG/dq_i without lambda

      // sum energy and sum1, sum2 for lambda
      esum += (V/2)*pow(fabs(atom->q[i]),p) + (U/2)*atom->q[i]*atom->q[i] + atom->q[i]*X;  // energy

      nand++;
    }
  }

  // check
  if(nand != ncq) {
    error->all(FLERR, "fix bsct: Bug: Number of charges doesn't match!");
  }

  // get global sums and max_dE
  {
    double temp;
    MPI_Allreduce(&esum, &temp, 1, MPI_DOUBLE, MPI_SUM, world);      // energy, band structure + Hubbard U + electrode chemical potential
    esum = temp;
    MPI_Allreduce(&qsum, &temp, 1, MPI_DOUBLE, MPI_SUM, world);      // charge within fix
    qsum = temp;
    MPI_Allreduce(&qtsum, &temp, 1, MPI_DOUBLE, MPI_SUM, world);    // total charge
    qtsum = temp;
  }

  // check
  if(ncq > 0) {
    if(qtsum/ncq > qconv) {
      char buf[100];
      sprintf(buf, "fix bsct: Total charge of system non-zero: %e", qtsum);
      error->warning(FLERR, buf);
    }
  }

  // get lambda
  get_lambda();

  // ---
  // --- finalize
  // ---

  // add lambda to the derivatives
  if(ncq > 0) {
    gsl_vector_add_constant(g, lambda);
  }

  // return total energy (BSCT + Coulomb)
  *f = esum + ecoul + lambda*(qsum - qtot);
}


/* ----------------------------------------------------------------------
   Update charges using Anderson mixing
   ---------------------------------------------------------------------- */

void FixBSCT::update_charges_anderson()
{

  // --- variables

  int i;                                    // loops
  int nand;                                 // number of local charges within fix
  double gmax;                              // maximum |dG/dq_i|
  double prevgmax;                          // previous value
  double gnorm;                             // |dG/dq|
  double temp;                              // temp
  double energyadd;                         // energy contribution to total energy (before updating charges, added to total energy of system)
  double energy;                            // energy contribution to total energy (after)
  int done;                                 // 1 - minimum found
  double qmax;                              // maximum |q_i|
  double phimax;                            // maximum |phi_i|
  double ebsct;                             // current BSCT energy

  gsl_vector *gsl_x;                        // charge vector for optimizers
  gsl_vector *gsl_g;                        // gradient vector

  // Anderson mixer
  AndersonMixer *amixer;                   // Anderson mixer
  double b = beta;                          // current mixing parameter

  // store forces and return after calculation
  //   (we don't want to zero forces since not all
  //    may be re-calculated by the calls below)
  int nbuf;
  double *fbuf;

  // --- store forces

  fbuf = new double[3*(atom->nlocal + atom->nghost) + 2];
  nbuf = 0;
  for(i=0;i<atom->nlocal + atom->nghost;i++) {
    fbuf[nbuf] = atom->f[i][0];
    nbuf++;
    fbuf[nbuf] = atom->f[i][1];
    nbuf++;
    fbuf[nbuf] = atom->f[i][2];
    nbuf++;
  }

  // --- allocate Anderson mixer

  amixer = new AndersonMixer(and_history);

  // --- setup stuff

  // Calculate contribution to total energy from band structure,
  // Hubbard U and electrode chemical potential terms, before updating
  // charges.
  //
  // also get number of local charges in fix
  {
    nand = 0;
    double esum = 0.0;
    for(i=0; i<atom->nlocal; i++) {
      if(atom->mask[i] & groupbit) {
        // params
        int itype = atom->type[i];
        if(itype >= nparams) error->all(FLERR, "fix bsct: Bug: Atom type number too big!");
        if(!params[itype].set) error->all(FLERR, "fix bsct: Parameters missing for atom type included in fix!");
        double X = params[itype].X;
        double U = params[itype].U;
        double V = params[itype].V;
        double p = params[itype].p;

        esum += (V/2)*pow(fabs(atom->q[i]),p) + (U/2)*atom->q[i]*atom->q[i] + atom->q[i]*X;
        nand++;
      }
    }
    MPI_Allreduce(&esum, &energyadd, 1, MPI_DOUBLE, MPI_SUM, world);
    ncq = nand;
  }

  // allocate arrays
  if(ncq > 0) {
    gsl_x = gsl_vector_alloc(ncq);  // charges, q_i
    gsl_g = gsl_vector_alloc(ncq);  // partial G / partial q_i
  }
  else {
    gsl_x = gsl_vector_alloc(1);
    gsl_g = gsl_vector_alloc(1);
  }

  // setup minimized charges
  nand = 0;
  for(i=0; i<atom->nlocal; i++) {
    if(atom->mask[i] & groupbit) {
      gsl_vector_set(gsl_x, nand, atom->q[i]);
      nand++;
    }
  }


  // --- main loop

  iter = 0;
  done = 0;

  if(log >= 1) fprintf(logfile,"                  CT   ------------------------------------\n");
  if(log >= 1) fprintf(logfile,"                  CT   iter               G    max|dG/dq_i|         |dG/dq|        max|q_i|      max|phi_i|            beta\n");

  // initial state for this step
  FixBSCT_fdf(gsl_x, &ebsct, gsl_g);  // also updates lambda

  // main loop
  gmax = -1.0;
  while(iter < max_iter && done==0) {
    // |dG/dq| and maximum |dG/dq_i|
    prevgmax = gmax;
    if(ncq > 0) {
      gsl_absmax_sqnorm(gsl_g, &gmax, &gnorm);
    }
    else {
      gmax = 0.0;
      gnorm = 0.0;
    }
    MPI_Allreduce(&gmax, &temp, 1, MPI_DOUBLE, MPI_MAX, world);
    gmax = temp;
    MPI_Allreduce(&gnorm, &temp, 1, MPI_DOUBLE, MPI_SUM, world);
    gnorm = sqrt(temp);

    // maximum |q_i|
    qmax = 0.0;
    phimax = 0.0;
    for(i=0; i<atom->nlocal; i++) {
      if(atom->mask[i] & groupbit) {
        if(fabs(atom->q[i]) > qmax) qmax = fabs(atom->q[i]);
        if(fabs(this->phi[i]) > phimax) phimax = fabs(this->phi[i]);
      }
    }
    MPI_Allreduce(&qmax, &temp, 1, MPI_DOUBLE, MPI_MAX, world);
    qmax = temp;
    MPI_Allreduce(&phimax, &temp, 1, MPI_DOUBLE, MPI_MAX, world);
    phimax = temp;

    // output
    if(log >= 1) fprintf(logfile,"                  CT  %5d %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f\n", iter, ebsct, gmax, gnorm, qmax, phimax, b);

    // test convergence
    if(fabs(gmax) < dGimax) done = 1;
    if(done==1) break;

    // new charges from iterative formula
    nand = 0;
    for(i=0; i<atom->nlocal; i++) {
      if(atom->mask[i] & groupbit) {
        // params
        double X = params[atom->type[i]].X;
        double U = params[atom->type[i]].U;
        double V = params[atom->type[i]].V;
        double p = params[atom->type[i]].p;

        gsl_vector_set(gsl_g, nand, nextq(atom->q[i], this->phi[i], lambda, X, U, V, p));
        nand++;
      }
    }
    if(nand != ncq) {
      error->all(FLERR, "fix bsct: Bug: Number of charges doesn't match!");
    }

    // mix
    amixer->mix(iter+1, nand, gsl_x, gsl_g, b, world, error);

    // Calculate new gradient (also extracts atom->q from gsl_x)
    temp = ebsct;  // store previous energy
    FixBSCT_fdf(gsl_x, &ebsct, gsl_g);

    // Update mixing
    //   normal mode: update if energy is going down
    if(beta_mul > 0.0) {
      if(ebsct < temp && fabs(gmax) < fabs(prevgmax)) b = b*beta_mul;
      if(prevgmax > 0 && fabs(gmax/prevgmax) > 2.0) b = beta;
      if(b > beta_max) b = beta_max;
    }
    //  for beta_mul < 0 => First two steps with beta, then beta_max
    else {
      if(ncall > 2) beta = beta_max;
    }

    iter++;
  }  // end of main loop

  // report if minimum found, stop with error if not
  if(done==1 && log >= 1) fprintf(logfile,"                  CT  %5d %15.5f %15.5f %15.5f  <- minimum found\n", iter, ebsct, gmax, gnorm);

  if(done==0) {
    error->warning(FLERR, "fix bsct: Charge transfer not converged!");
  }

  if(log >= 1) fprintf(logfile,"                  CT   ------------------------------------\n");


  // --- unallocate GSL

  gsl_vector_free(gsl_x);
  gsl_vector_free(gsl_g);

  // --- return forces

  nbuf = 0;
  for(i=0;i<atom->nlocal + atom->nghost;i++) {
    atom->f[i][0] = fbuf[nbuf];
    nbuf++;
    atom->f[i][1] = fbuf[nbuf];
    nbuf++;
    atom->f[i][2] = fbuf[nbuf];
    nbuf++;
  }

  // --- free memory

  delete [] fbuf;
  delete amixer;
}


/* ----------------------------------------------------------------------
   Norm and maximum value of a GSL vector
   ---------------------------------------------------------------------- */

void FixBSCT::gsl_absmax_sqnorm(gsl_vector *x, double *absmax, double *sqnorm) {

  int i;       // loops
  double sum;  // sum for norm squared
  double max;  // max
  double v;    // current value

  sum = 0.0;
  max = 0.0;
  for(i=0;i<x->size;i++) {
    v = gsl_vector_get(x, i);
    if(fabs(v) > max) max = fabs(v);
    sum += v*v;
  }

  *absmax = max;
  *sqnorm = sum;
}


/* ----------------------------------------------------------------------
     OTHER ROUTINES
   ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Adjust total charge to whatever value specified on input
   ---------------------------------------------------------------------- */

void FixBSCT::adjust_total_charge() {

  // --- variables
  int i;  // loops
  double local_qsum;  // sum of charges (within fix) locally
  double global_qsum;  // sum of charges (within fix) globally
  int local_n;  // number of charges summed locally
  int global_n;  // number of charges summed globally

  // --- adjust

  // sum local charges
  local_qsum = 0.0;
  local_n = 0;
  for(i=0;i<atom->nlocal;i++) {
    if(atom->mask[i] & groupbit) {
      local_qsum += atom->q[i];
      local_n++;
    }
  }

  // get sum over charges globally
  MPI_Allreduce(&local_qsum, &global_qsum, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&local_n, &global_n, 1, MPI_INT, MPI_SUM, world);

  // adjust charges
  for(i=0;i<atom->nlocal;i++) {
    if(atom->mask[i] & groupbit) {
      atom->q[i] = atom->q[i] - (global_qsum - qtot)/global_n;
    }
  }
}


/* ----------------------------------------------------------------------
   PER-ATOM DATA STORED IN FIX
   ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   Memory usage of local atom-based arrays
   ---------------------------------------------------------------------- */

double FixBSCT::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;

  // phi
  bytes += nmax * sizeof(double);

  return bytes;
}


/* ----------------------------------------------------------------------
   Allocate memory for local per-atom arrays
   ---------------------------------------------------------------------- */

void FixBSCT::grow_arrays(int nmax)
{
  // Always phi
  memory->grow(phi,nmax,"bsct:phi");
  vector_atom = phi;
}


/* ----------------------------------------------------------------------
   Copy values within local atom-based arrays
   ---------------------------------------------------------------------- */

void FixBSCT::copy_arrays(int i, int j)
{
  phi[j] = phi[i];
}


/* ----------------------------------------------------------------------
   Initialize one atom's array values, called when atom is created
   ---------------------------------------------------------------------- */

void FixBSCT::set_arrays(int i)
{
  phi[i] = 0.0;
}


/* ----------------------------------------------------------------------
   Pack values in local atom-based arrays for exchange with another proc
   ---------------------------------------------------------------------- */

int FixBSCT::pack_exchange(int i, double *buf)
{
  buf[0] = phi[i];

  return 1;
}


/* ----------------------------------------------------------------------
   Unpack values in local atom-based arrays from exchange with another proc
   ---------------------------------------------------------------------- */

int FixBSCT::unpack_exchange(int nlocal, double *buf)
{
  phi[nlocal] = buf[0];

  return 1;
}


/* ----------------------------------------------------------------------
   Forward communication: Charges (from Atom)
   ---------------------------------------------------------------------- */

int FixBSCT::pack_comm(int n, int *list, double *buf,
                       int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = atom->q[j];
  }
  return 1;
}


/* ----------------------------------------------------------------------
   Forward communication: unpack
   ---------------------------------------------------------------------- */

void FixBSCT::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    atom->q[i] = buf[m++];
  }
}


/* ----------------------------------------------------------------------
   Reverse communication: Electrostatic potentials
   ---------------------------------------------------------------------- */

int FixBSCT::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = this->phi[i];
  }
  return 1;
}


/* ----------------------------------------------------------------------
   Reverse communication: unpack
   ---------------------------------------------------------------------- */

void FixBSCT::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    this->phi[j] += buf[m++];
  }
}
