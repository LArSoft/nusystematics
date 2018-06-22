#ifndef NUSYST_RESPONSE_HELPER_SEEN
#define NUSYST_RESPONSE_HELPER_SEEN

#include "nusyst/artless/configure_nusyst_providers.hh"
#include "nusyst/interface/IGENIESystProvider_tool.hh"
#include "nusyst/interface/types.hh"

#include "larsyst/interface/types.hh"
#include "larsyst/interpreters/ParamHeaderHelper.hh"
#include "larsyst/interpreters/load_parameter_headers.hh"

#include "fhiclcpp/make_ParameterSet.h"

#include "EVGCore/EventRecord.h"

#include <string>

namespace nusyst {
class response_helper : public larsyst::ParamHeaderHelper {

private:
  std::string config_file;
  genie_provider_map_t syst_providers;
  larsyst::param_header_map_t configuredParameterHeaders;

  void LoadProviders(fhicl::ParameterSet const &ps) {
    syst_providers = nusyst::load_syst_provider_configuration(ps);

    if (!syst_providers.size()) {
      std::cout
          << "[ERROR]: Expected to load some systematic providers from input: "
          << std::quoted(config_file) << std::endl;
      throw;
    }
  }

  void LoadHeaders(fhicl::ParameterSet const &ps) {
    configuredParameterHeaders = larsyst::load_syst_provider_headers(ps);
    if (!configuredParameterHeaders.size()) {
      std::cout << "[ERROR]: Expected to find some parameter headers in "
                << std::quoted(config_file) << std::endl;
      throw;
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

    LoadProviders(ps);
    LoadHeaders(ps);
  }

  larsyst::event_unit_response_t
  GetEventResponses(genie::EventRecord &GenieGHep) {
    larsyst::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      larsyst::event_unit_response_t prov_response =
          sp.second->GetEventResponse(GenieGHep);
      for (auto &&er : prov_response) {
        response.insert(std::move(er));
      }
    }
    return response;
  }
};
} // namespace nusyst

#endif
