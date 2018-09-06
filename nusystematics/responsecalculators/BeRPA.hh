#ifndef nusystematics_RESPONSE_CALCULATORS_BERPA_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_BERPA_HH_SEEN

#include <cmath>

namespace nusyst {

/// C. Wilkinson's fitted nominals
/// A = 0.59 +/- 20%
/// B = 1.05 +/- 20%
/// D = 1.13 +/- 15%
/// E = 0.88 +/- 40%
/// U = 1.2

double BeRPA(double Q2_GeV2, double A = 0.59, double B = 1.05,
               double D = 1.13, double E = 0.88, double U = 1.20) {
  // Kept for convenience
  double eRPA = 1.;

  // Q2 transition; if Q2 less than U -> polynominal, otherwise exponential
  // decay
  if (Q2_GeV2 < U) {
    // xprime as prescribed by Callum
    const double xprime = Q2 / U;
    const double one_minus_xprime = 1. - xprime;
    const double C = D + U * E * (D - 1) / 3.;
    eRPA = A * one_minus_xprime * one_minus_xprime * one_minus_xprime +
           3 * B * one_minus_xprime * one_minus_xprime * xprime +
           3 * C * one_minus_xprime * xprime * xprime +
           D * xprime * xprime * xprime;
  } else {
    eRPA = 1 + (D - 1) * std::exp(-E * (Q2 - U));
  }

  return eRPA;
}
} // namespace nusyst

#endif
