#ifndef nusystematics_MAKE_INSTANCE_SEEN
#define nusystematics_MAKE_INSTANCE_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/GENIEReWeight_tool.hh"
#include "nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh"
#include "nusystematics/systproviders/MKSinglePiTemplate_tool.hh"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(unknown_nusyst_systprovider);

inline std::unique_ptr<IGENIESystProvider_tool>
make_instance(fhicl::ParameterSet const &paramset) {
  std::string tool_type = paramset.get<std::string>("tool_type");

  if (tool_type == "GENIEReWeight") {
    return std::make_unique<GENIEReWeight>(paramset);
  } else if (tool_type == "MKSinglePiTemplate") {
    return std::make_unique<MKSinglePiTemplate>(paramset);
  } else if (tool_type == "MINERvAq0q3Weighting") {
    return std::make_unique<MINERvAq0q3Weighting>(paramset);
  } else {
    throw unknown_nusyst_systprovider()
        << "[ERROR]: Unknown tool type: " << std::quoted(tool_type);
  }
}

} // namespace nusyst

#endif
