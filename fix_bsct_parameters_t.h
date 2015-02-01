
#ifndef LMP_FIX_BSCT_PARAMETERS_H
#define LMP_FIX_BSCT_PARAMETERS_H

namespace LAMMPS_NS {

  /*
    Collection of charge transfer model parameters for one atom type.
  */
  struct FixBSCT_parameters_t {

    bool set;         // parameters set?

    double X;         // chemical potential of electrode
    double U;         // Hubbard U
    double V;         // Band structure energy prefactor
    double p;         // Band structure energy exponential
  };
}

#endif
