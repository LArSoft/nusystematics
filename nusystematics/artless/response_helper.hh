#ifndef nusystematics_RESPONSE_HELPER_SEEN
#define nusystematics_RESPONSE_HELPER_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/GENIEReWeight_tool.hh"
#include "nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh"
#include "nusystematics/systproviders/MKSinglePiEnuq0q3_tool.hh"

#include "systematicstools/interface/SystParamHeader.hh"

#include "systematicstools/interpreters/ParamHeaderHelper.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"

#include "fhiclcpp/make_ParameterSet.h"

#include "EVGCore/EventRecord.h"

#include "TFile.h"
#include "TTree.h"

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

NEW_SYSTTOOLS_EXCEPT(response_helper_found_no_parameters);

class response_helper : public systtools::ParamHeaderHelper {

private:
  constexpr static size_t Order = 5;
  constexpr static size_t NCoeffs = Order + 1;

  std::string config_file;
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> syst_providers;

  void LoadProvidersAndHeaders(fhicl::ParameterSet const &ps) {
    syst_providers = systtools::ConfigureISystProvidersFromParameterHeaders<
        IGENIESystProvider_tool>(ps, make_instance);

    if (!syst_providers.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected to load some systematic providers from input: "
          << std::quoted(config_file);
    }

    systtools::param_header_map_t configuredParameterHeaders =
        systtools::BuildParameterHeaders(syst_providers);
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

    fhicl::ParameterSet ps = fhicl::make_ParameterSet(config_file);
    LoadProvidersAndHeaders(ps.get<fhicl::ParameterSet>(
        "generated_systematic_provider_configuration"));
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord &GenieGHep) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord &GenieGHep, systtools::paramId_t i) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep, i);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  double GetEventWeightResponse(genie::EventRecord &GenieGHep,
                                systtools::param_value_list_t const &vals) {
    double weight = 1;
    for (auto &sp : syst_providers) {
      weight *= sp->GetEventWeightResponse(GenieGHep, vals);
    }
    return weight;
  }
};
} // namespace nusyst

#endif
