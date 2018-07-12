#ifndef NUSYST_CONFIGURESYSTPROVIDERS_SEEN
#define NUSYST_CONFIGURESYSTPROVIDERS_SEEN

#include "nusyst/interface/IGENIESystProvider_tool.hh"
#include "nusyst/interface/types.hh"
#include "nusyst/systproviders/GENIEReWeight_tool.hh"
#include "nusyst/systproviders/MINERvAq0q3Weighting_tool.hh"
#include "nusyst/systproviders/MKSinglePiEnuq0q3_tool.hh"

#include "larsyst/interface/types.hh"

#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <string>

namespace nusyst {

inline std::unique_ptr<IGENIESystProvider_tool>
make_instance(fhicl::ParameterSet const &paramset) {
  std::string tool_type = paramset.get<std::string>("tool_type");

  if (tool_type == "GENIEReWeight") {
    return std::make_unique<GENIEReWeight>(paramset);
  } else if (tool_type == "MKSinglePiEnuq0q3") {
    return std::make_unique<MKSinglePiEnuq0q3>(paramset);
  } else if (tool_type == "MINERvAq0q3Weighting") {
    return std::make_unique<MINERvAq0q3Weighting>(paramset);
  } else {
    std::cout << "[ERROR]: Unknown tool type: " << std::quoted(tool_type)
              << std::endl;
    throw;
  }
}

inline genie_provider_map_t
load_syst_provider_configuration(fhicl::ParameterSet const &paramset,
                                 std::string const &key = "syst_providers") {

  genie_provider_map_t loaded_providers;

  size_t nparams = 0;
  std::cout << "[INFO]: Loading configured syst providers:" << std::endl;
  for (auto const &provider_name :
       paramset.get<std::vector<std::string>>(key)) {
    // Get fhicl config for provider
    std::cout << "[INFO]:\t Retrieving meta data for: \"" << provider_name
              << "\"..." << std::endl;
    auto const &provider_cfg = paramset.get<fhicl::ParameterSet>(provider_name);
    std::cout << " found!" << std::endl;

    // make an instance of the plugin
    std::cout << "[INFO]:\t Requesting provider instance..." << std::endl;
    std::unique_ptr<IGENIESystProvider_tool> is = make_instance(provider_cfg);
    std::cout << " success!" << std::endl;

    // configure the plugin
    std::cout << "[INFO]:\t Loading provider configuration..." << std::endl;

    if (is->ReadParameterHeaders(provider_cfg)) {
      std::cout << "[INFO]\t Success!" << std::endl;
    } else {
      std::cout << "[ERROR]:\t Failure." << std ::endl;
      throw;
    }
    nparams += is->GetSystSetConfiguration().headers.size();

    // build unique name
    std::string FQName = is->GetFullyQualifiedName();
    std::cout
        << "[INFO]:\t Attempting to register provider with fully qualifed "
           "name: \""
        << FQName << "\"..." << std::endl;

    // check that this unique name hasn't been used before.
    if (loaded_providers.find(FQName) != loaded_providers.end()) {
      std::cout << "failed." << std::endl
                << "[ERROR]:\t Provider with that name already exists, please "
                   "correct provider set (Hint: Use the 'unique_name' property "
                   "of the tool configuration table to dismabiguate multiple "
                   "uses of the same tool)."
                << std::endl;
      throw;
    }
    std::cout << "Success!" << std::endl;
    // std::cout << "[INFO]: Configured " << is->AsString() << std::endl;
    loaded_providers.emplace(std::move(FQName), std::move(is));
  }
  std::cout << "[INFO]: Loaded " << loaded_providers.size()
            << " systematic providers with " << nparams << " parameters."
            << std::endl;
  return loaded_providers;
}

} // namespace nusyst

#endif
