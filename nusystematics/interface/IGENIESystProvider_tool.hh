#ifndef nusystematics_INTERFACE_IGENIESYSTPROVIDER_TOOL_SEEN
#define nusystematics_INTERFACE_IGENIESYSTPROVIDER_TOOL_SEEN

#include "systematicstools/interface/ISystProviderTool.hh"

#ifndef NO_ART
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#endif

#include "fhiclcpp/ParameterSet.h"

// GENIE includes
#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "EVGCore/EventRecord.h"
#else
  // Use these for GENIE v3
  #include "Framework/EventGen/EventRecord.h"
#endif

namespace nusyst {

class IGENIESystProvider_tool : public systtools::ISystProviderTool {
public:
  IGENIESystProvider_tool(fhicl::ParameterSet const &ps)
      : ISystProviderTool(ps), fGENIEModuleLabel(ps.get<std::string>(
                                   "genie_module_label", "generator")) {}

  NEW_SYSTTOOLS_EXCEPT(invalid_response);

#ifndef NO_ART
  std::unique_ptr<systtools::EventResponse>
  GetEventResponse(art::Event const &ev) {
    std::unique_ptr<systtools::EventResponse> er =
        std::make_unique<systtools::EventResponse>();

    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    art::Handle<std::vector<simb::GTruth>> gTruthHandle;
    ev.getByLabel(fGENIEModuleLabel, mcTruthHandle);
    ev.getByLabel(fGENIEModuleLabel, gTruthHandle);

    size_t NEventUnits = mcTruthHandle->size();
    if (mcTruthHandle->size() != gTruthHandle->size()) {
      NEventUnits = std::min(mcTruthHandle->size(), gTruthHandle->size());
    }

    std::vector<std::unique_ptr<genie::EventRecord>> gheps;
    for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
      gheps.emplace_back(evgb::RetrieveGHEP(mcTruthHandle->at(eu_it),
                                            gTruthHandle->at(eu_it)));
    }

    er->resize(NEventUnits);
    for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
      er->push_back(GetEventResponse(*gheps[eu_it]));
    }
    return er;
  }
#endif

  /// Calculates configured response for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &) = 0;

  systtools::event_unit_response_w_cv_t
  GetEventVariationAndCVResponse(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_w_cv_t responseandCV;

    systtools::event_unit_response_t prov_response =
        GetEventResponse(GenieGHep);

    // Foreach param
    for (systtools::ParamResponses &pr : prov_response) {
      // Get CV resp
      systtools::SystParamHeader const &hdr =
          GetParam(GetSystMetaData(), pr.pid);

      if (pr.responses.size() != hdr.paramVariations.size()) {
        throw invalid_response()
            << "[ERROR]: Parameter: " << hdr.prettyName << ", with "
            << hdr.paramVariations.size() << " parameter variations, returned "
            << pr.responses.size() << " responses.";
      }

      double CVResp = hdr.isWeightSystematicVariation ? 1 : 0;
      size_t NVars = hdr.paramVariations.size();

      double cv_param_val = 0;
      if (hdr.centralParamValue != systtools::kDefaultDouble) {
        cv_param_val = hdr.centralParamValue;
      }
      for (size_t idx = 0; idx < NVars; ++idx) {
        if (fabs(cv_param_val - hdr.paramVariations[idx]) <=
            std::numeric_limits<float>::epsilon()) {
          CVResp = pr.responses[idx];
          break;
        }
      }
      // if we didn't find it, the CVResp stays as 1/0 depending on whether it
      // is a weight or not.
      for (size_t idx = 0; idx < NVars; ++idx) {
        if (hdr.isWeightSystematicVariation) {
          pr.responses[idx] /= CVResp;
        } else {
          pr.responses[idx] -= CVResp;
        }
      }

      responseandCV.push_back({pr.pid, CVResp, pr.responses});
    } // end for parameter response

    return responseandCV;
  }

  /// Calculates the response to a single parameter for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &, systtools::paramId_t) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement systtools::event_unit_response_t "
           "GetEventResponse(genie::EventRecord &, systtools::paramId_t).";
  }

  /// Calculates the multiplicatively combined responses for a given set of
  /// parameter--value pairs.
  ///
  /// \note This convenience method should only be used for weight responses.
  virtual double GetEventWeightResponse(genie::EventRecord const &,
                                        systtools::param_value_list_t const &) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement double "
           "GetEventWeightResponse(genie::EventRecord "
           "&,systtools::param_value_list_t const &).";
  }

  std::string fGENIEModuleLabel;
};
} // namespace nusyst

#endif
