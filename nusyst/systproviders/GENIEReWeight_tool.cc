#include "nusyst/systproviders/GENIEReWeight_tool.hh"

#include "nusyst/systproviders/GENIEReWeightEngineConfig.hh"
#include "nusyst/systproviders/GENIEReWeightParamConfig.hh"

#ifndef NO_ART
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#endif

#include "larsyst/utility/printers.hh"
#include "larsyst/utility/string_parsers.hh"

#ifndef NO_ART
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "art/Utilities/ToolMacros.h"
#endif

// GENIE
#include "Messenger/Messenger.h"

#include <chrono>
#include <sstream>

using namespace fhicl;
using namespace larsyst;
using namespace nusyst;

GENIEReWeight::GENIEReWeight(ParameterSet const &params)
    : IGENIESystProvider_tool(params) {}

std::string GENIEReWeight::AsString() {
  CheckHaveMetaData();
  return "";
}

SystMetaData GENIEReWeight::BuildSystMetaData(ParameterSet const &params,
                                              paramId_t firstParamId) {

  tool_options = fhicl::ParameterSet();

  bool ignore_parameter_dependence =
      params.get<bool>("ignore_parameter_dependence", false);

  tool_options.put("ignore_parameter_dependence", ignore_parameter_dependence);

  SystMetaData QEmd =
      ConfigureQEParameterHeaders(params, firstParamId, tool_options);
  firstParamId += QEmd.size();

#ifndef GRWTEST
  SystMetaData NCELmd =
      ConfigureNCELParameterHeaders(params, firstParamId, tool_options);
  firstParamId += NCELmd.size();
  SystMetaData RESmd =
      ConfigureRESParameterHeaders(params, firstParamId, tool_options);
  firstParamId += RESmd.size();
  SystMetaData COHmd =
      ConfigureCOHParameterHeaders(params, firstParamId, tool_options);
  firstParamId += COHmd.size();
  SystMetaData DISmd =
      ConfigureDISParameterHeaders(params, firstParamId, tool_options);
  firstParamId += DISmd.size();
  SystMetaData FSImd =
      ConfigureFSIParameterHeaders(params, firstParamId, tool_options);
  firstParamId += FSImd.size();
  SystMetaData Othermd =
      ConfigureOtherParameterHeaders(params, firstParamId, tool_options);
  firstParamId += Othermd.size();

  // Don't extend inline to make firstParamId incrementing more clear.
  extend_SystMetaData(QEmd, NCELmd);
  extend_SystMetaData(QEmd, RESmd);
  extend_SystMetaData(QEmd, COHmd);
  extend_SystMetaData(QEmd, DISmd);
  extend_SystMetaData(QEmd, FSImd);
  extend_SystMetaData(QEmd, Othermd);
#endif

  return QEmd;
}

void GENIEReWeight::extend_ResponseToGENIEParameters(
    std::vector<GENIEResponseParameter> &&other) {
  for (auto &&o : other) {
    for (auto const &configured : ResponseToGENIEParameters) {
      if (configured.pidx == o.pidx) {
        std::cout << "[ERROR]: Attempted to merge GENIE GSyst response map, "
                     "but found duplicate response parameter index: "
                  << o.pidx << ", which corresponds to parameter: "
                  << std::quoted(GetSystMetaData()[o.pidx].prettyName)
                  << std::endl;
        throw;
      }
    }
    ResponseToGENIEParameters.push_back(std::move(o));
  }
}

bool GENIEReWeight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  extend_ResponseToGENIEParameters(
      ConfigureQEWeightEngine(GetSystMetaData(), tool_options));

#ifndef GRWTEST
  extend_ResponseToGENIEParameters(
      ConfigureNCELWeightEngine(GetSystMetaData(), tool_options));
  extend_ResponseToGENIEParameters(
      ConfigureRESWeightEngine(GetSystMetaData(), tool_options));
  extend_ResponseToGENIEParameters(
      ConfigureCOHWeightEngine(GetSystMetaData(), tool_options));
  extend_ResponseToGENIEParameters(
      ConfigureDISWeightEngine(GetSystMetaData(), tool_options));
  extend_ResponseToGENIEParameters(
      ConfigureFSIWeightEngine(GetSystMetaData(), tool_options));
  extend_ResponseToGENIEParameters(
      ConfigureOtherWeightEngine(GetSystMetaData(), tool_options));
#endif

  return true;
}

#ifndef NO_ART
std::unique_ptr<EventResponse> GENIEReWeight::GetEventResponse(art::Event &e) {
  std::unique_ptr<EventResponse> er = std::make_unique<EventResponse>();

  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
  art::Handle<std::vector<simb::GTruth>> gTruthHandle;
  e.getByLabel(fGENIEModuleLabel, mcTruthHandle);
  e.getByLabel(fGENIEModuleLabel, gTruthHandle);

  size_t NEventUnits = mcTruthHandle->size();
  if (mcTruthHandle->size() != gTruthHandle->size()) {
    NEventUnits = std::min(mcTruthHandle->size(), gTruthHandle->size());
  }

  std::vector<std::unique_ptr<genie::EventRecord>> gheps;
  for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
    gheps.emplace_back(
        evgb::RetrieveGHEP(mcTruthHandle->at(eu_it), gTruthHandle->at(eu_it)));
  }

  er->responses.resize(NEventUnits);
  for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
    er->responses.push_back(GetEventResponse(*gheps[eu_it]));
  }
  return er;
}

DEFINE_ART_CLASS_TOOL(GENIEReWeight)

#endif

larsyst::event_unit_response_t
GENIEReWeight::GetEventResponse(genie::EventRecord &gev) {

  larsyst::event_unit_response_t event_responses;

  for (auto &GENIEResponse : ResponseToGENIEParameters) {
    larsyst::SystParamHeader const &hdr = GetSystMetaData()[GENIEResponse.pidx];
    size_t NVars = hdr.isCorrection ? 1 : hdr.paramVariations.size();
#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
    std::cout << "[INFO]: Resp dial: " << hdr.prettyName << " with " << NVars
              << " variations of " << GENIEResponse.dependents.size()
              << std::endl;
#endif

    // Have one GENIEReWeight per response rather than per variation, must
    // reconfigure.
    bool IsReducedHERG = (NVars > GENIEResponse.Herg.size());
    event_responses.push_back({hdr.systParamId, {}});
    for (size_t var_it = 0; var_it < NVars; ++var_it) {

      if (IsReducedHERG) { // Need a reconfigure for each variation
        for (auto const &dep : GENIEResponse.dependents) {
          GENIEResponse.Herg.front()->Systematics().Set(
              dep.gdial, GetSystMetaData()[dep.pidx].paramVariations[var_it]);
        }
        GENIEResponse.Herg.front()->Reconfigure();
        event_responses.back().responses.push_back(
            GENIEResponse.Herg.front()->CalcWeight(gev));
      } else {
        event_responses.back().responses.push_back(
            GENIEResponse.Herg[var_it]->CalcWeight(gev));
      }

#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
      std::cout << "\t Var = " << hdr.paramVariations[var_it] << " GDial: "
                << genie::rew::GSyst::AsString(
                       GENIEResponse.dependents.front().gdial)
                << " at "
                << GENIEResponse.Herg[var_it]
                       ->Systematics()
                       .Info(GENIEResponse.dependents.front().gdial)
                       ->CurValue
                << " -> " << event_responses.back().responses.back()
                << std::endl;
#endif
    }
  }
  return event_responses;
}
