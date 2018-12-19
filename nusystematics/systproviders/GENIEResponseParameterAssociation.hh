#ifndef nusystematics_SYSTPROVIDERS_GENIERESPONSEPARAMETERASSOCIATION_SEEN
#define nusystematics_SYSTPROVIDERS_GENIERESPONSEPARAMETERASSOCIATION_SEEN

#include "systematicstools/utility/exceptions.hh"

// GENIE includes
#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "ReWeight/GReWeight.h"
#else
  // Use these for GENIE v3
  #include "RwFramework/GReWeight.h"
#endif

#include <memory>
#include <vector>

namespace nusyst {

typedef size_t parameter_idx_t;

struct GENIEResponseParameter {
  struct DependentParameter {
    genie::rew::GSyst_t gdial;
    parameter_idx_t pidx;
  };
  parameter_idx_t pidx;
  std::vector<DependentParameter> dependents;
  std::vector<std::unique_ptr<genie::rew::GReWeight>> Herg;
};

NEW_SYSTTOOLS_EXCEPT(invalid_GENIE_parameter_index);

template <typename T>
GENIEResponseParameter &GetResponseParameter(T &RespContainer, parameter_idx_t pidx) {
  for (GENIEResponseParameter &pr : RespContainer) {
    if (pr.pidx == pidx) {
      return pr;
    }
  }
  throw invalid_GENIE_parameter_index();
}

} // namespace nusyst
#endif
