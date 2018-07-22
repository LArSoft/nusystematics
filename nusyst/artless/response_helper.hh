#ifndef NUSYST_RESPONSE_HELPER_SEEN
#define NUSYST_RESPONSE_HELPER_SEEN

#include "nusyst/interface/IGENIESystProvider_tool.hh"

#include "nusyst/systproviders/GENIEReWeight_tool.hh"

#include "larsyst/interface/SystParamHeader.hh"

#include "larsyst/interpreters/ParamHeaderHelper.hh"

#include "larsyst/utility/ParameterAndProviderConfigurationUtility.hh"

#include "fhiclcpp/make_ParameterSet.h"

#include "EVGCore/EventRecord.h"

#include <string>

namespace nusyst {

inline std::unique_ptr<IGENIESystProvider_tool>
make_instance(fhicl::ParameterSet const &paramset) {
  std::string tool_type = paramset.get<std::string>("tool_type");

  if (tool_type == "GENIEReWeight") {
    return std::make_unique<GENIEReWeight>(paramset);
    // } else if (tool_type == "MKSinglePiEnuq0q3") {
    //   return std::make_unique<MKSinglePiEnuq0q3>(paramset);
    // } else if (tool_type == "MINERvAq0q3Weighting") {
    //   return std::make_unique<MINERvAq0q3Weighting>(paramset);
  } else {
  std:
    std::cout << "[ERROR]: Unknown tool type: " << std::quoted(tool_type)
              << std::endl;
    throw;
  }
}

NEW_LARSYST_EXCEPT(response_helper_found_no_parameters);

class response_helper : public larsyst::ParamHeaderHelper {

private:
  std::string config_file;
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> syst_providers;
  larsyst::param_header_map_t configuredParameterHeaders;

  void LoadProvidersAndHeaders(fhicl::ParameterSet const &ps) {
    syst_providers = larsyst::ConfigureISystProvidersFromParameterHeaders<
        IGENIESystProvider_tool>(ps, make_instance);

    if (!syst_providers.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected to load some systematic providers from input: "
          << std::quoted(config_file);
    }

    configuredParameterHeaders = larsyst::BuildParameterHeaders(syst_providers);
    if (!configuredParameterHeaders.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected systematric providers loaded from input: "
          << std::quoted(config_file) << " to provide some parameter headers.";
    }

    SetHeaders(configuredParameterHeaders);
  }

public:
  response_helper() {}
  response_helper(std::string const &fhicl_config_filename) {
    LoadConfiguration(fhicl_config_filename);
  }

  void LoadConfiguration(std::string const &fhicl_config_filename) {
    config_file = fhicl_config_filename;

    fhicl::ParameterSet ps =
        fhicl::make_ParameterSet(config_file)
            .get<fhicl::ParameterSet>(
                "generated_systematic_provider_configuration");

    LoadProvidersAndHeaders(ps);
  }

  larsyst::event_unit_response_t
  GetEventResponses(genie::EventRecord &GenieGHep) {
    larsyst::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      larsyst::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }
};
} // namespace nusyst

#endif
