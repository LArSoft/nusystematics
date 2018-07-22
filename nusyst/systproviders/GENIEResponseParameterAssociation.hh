#ifndef NUSYST_SYSTPROVIDERS_GENIERESPONSEPARAMETERASSOCIATION_SEEN
#define NUSYST_SYSTPROVIDERS_GENIERESPONSEPARAMETERASSOCIATION_SEEN

#include "larsyst/utility/exceptions.hh"

// GENIE
#include "ReWeight/GReWeight.h"

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

NEW_LARSYST_EXCEPT(invalid_GENIE_parameter_index);

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
