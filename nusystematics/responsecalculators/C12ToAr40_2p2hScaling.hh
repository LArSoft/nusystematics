#ifndef nusystematics_RESPONSE_CALCULATORS_CAr2p2hWeight_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CAr2p2hWeight_HH_SEEN

#include "nusystematics/utility/GENIEUtils.hh"

namespace nusyst {

inline double GetC_Ar2p2hScalingWeight(QELikeTarget_t mec_topology,
                                       double parameter_value = 0) {

  switch (mec_topology) {
  case kNN: {
    static double const central_value = 1.33;
    static double const uncertainty = 0.13;
    return (central_value + parameter_value * uncertainty);
  }
  case knp: {
    static double const central_value = 0.9;
    static double const uncertainty = 0.4;
    return (central_value + parameter_value * uncertainty);
  }
  default: { return 1; }
  }
}

} // namespace nusyst

#endif
