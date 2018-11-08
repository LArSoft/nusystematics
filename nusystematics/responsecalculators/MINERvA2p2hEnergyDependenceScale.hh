#ifndef nusystematics_RESPONSE_CALCULATORS_MINERvA2p2EnergyDependencyScaling_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_MINERvA2p2EnergyDependencyScaling_HH_SEEN

#include "nusystematics/utility/GENIEUtils.hh"

namespace nusyst {

static double const A_OneSig = 0.9; // GeV2
static double const B_OneSig = 0.3; // GeV

static double const A_CV = 0;
static double const B_CV = 0;

inline double Get_MINERvA2p2h2EnergyDependencyScaling(int mode_simb, int is_CC,
                                                      double nu_Energy_GeV,
                                                      double A_val = 0,
                                                      double B_val = 0) {

  if (!is_CC) {
    return 1;
  }

  if (mode_simb != e2i(simb_mode_copy::kMEC)) {
    return 1;
  }

  double A = A_CV + A_val * A_OneSig;
  double B = B_CV + B_val * B_OneSig;

  double ScaleFactor = ((A / pow(nu_Energy_GeV, 2)) + (B / nu_Energy_GeV));

  return 1.0 / ScaleFactor;
}

} // namespace nusyst

#endif
