#ifndef nusystematics_RESPONSE_CALCULATORS_SPPLOWQ2SUPPRESSION_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_SPPLOWQ2SUPPRESSION_HH_SEEN

#include "nusystematics/utility/enumclas2int.hxx"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

///
///
///\note From Phys. Rev. D 91, 012005 (2015) Fig. 9
inline double GetSPPLowQ2SuppressionWeight(int mode_simb, int is_CC,
                                           double Q2_GeV,
                                           double parameter_value = 1) {

  static double const A = 1.01;
  static double const Q2_Max = 0.7;
  static double const Q0 = 0.156;
  static double const central_value = 0;
  static double const uncertainty = 1;

  if ((mode_simb != e2i(simb_mode_copy::kRES)) || !is_CC) {
    return 1;
  }

  if ((Q2_GeV > Q2_Max) || (Q2_GeV < 0)) {
    return 1;
  }

  return (central_value + parameter_value * uncertainty) *
         (A / (1 + exp(1 - (sqrt(Q2_GeV) / Q0))));
}
} // namespace nusyst

#endif
