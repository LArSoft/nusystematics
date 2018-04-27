#include "GENIEReWeight_config.hh"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larsyst/interface/ISystProvider_tool.hh"
#include "larsyst/utility/printers.hh"
#include "larsyst/utility/string_parsers.hh"

#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"

#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <chrono>
#include <sstream>

using namespace larsyst;
using namespace fhicl;

class GENIEReWeight : public larsyst::ISystProvider_tool {
public:
  explicit GENIEReWeight(ParameterSet const &);

  SystMetaData ConfigureFromFHICL(ParameterSet const &, paramId_t);

  bool Configure();
  std::unique_ptr<EventResponse> GetEventResponse(art::Event &);
  std::string AsString();

private:
  std::unique_ptr<genie::rew::GReWeight> GReWeightEngine;

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
      ResponseToGENIEParameters;
  std::set<genie::rew::GSyst_t> GENIEEngineDials;
  std::string fGENIEModuleLabel = "generator";
};

GENIEReWeight::GENIEReWeight(ParameterSet const &params)
    : ISystProvider_tool(params), GReWeightEngine{nullptr} {}

std::string GENIEReWeight::AsString() {
  CheckHaveMetaData();
  return "";
}

SystMetaData GENIEReWeight::ConfigureFromFHICL(ParameterSet const &params,
                                               paramId_t firstParamId) {

  std::cout << "[INFO]: Configuring GENIEReWeight" << std::endl;

  Table<GENIEReWeightConfig> cfg{
      params, std::set<std::string>{"tool_type", "uniqueName"}};

  SystMetaData QEmd = ConfigureQEParameterHeaders(cfg, firstParamId);
  firstParamId += QEmd.headers.size();
  return QEmd;
}

bool GENIEReWeight::Configure() {
  GReWeightEngine = std::make_unique<genie::rew::GReWeight>();
  genie::rew::GSystSet &gdials = GReWeightEngine->Systematics();
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  ResponseToGENIEParameters =
      ConfigureQEWeightEngine(fMetaData, GReWeightEngine);

  std::cout << "[INFO]: NuSyst provider "
            << std::quoted(GetFullyQualifiedName()) << " configured with "
            << ResponseToGENIEParameters.size() << " response parameters. "
            << std::endl;

  for (auto const &resp : ResponseToGENIEParameters) {
    std::cout << "\tParameter #" << resp.first << ":"
              << std::quoted(fMetaData.headers[resp.first].prettyName)
              << " uses " << resp.second.size()
              << " GENIE parameters: " << std::endl;
    for (auto const &gparpair : resp.second) {
      std::cout << "\t\t" << genie::rew::GSyst::AsString(gparpair.first) << ", "
                << fMetaData.headers[gparpair.second].systParamId << ":"
                << std::quoted(fMetaData.headers[gparpair.second].prettyName)
                << std::endl;
      gdials.Init(gparpair.first);
      GENIEEngineDials.insert(gparpair.first);
    }
  }

  // larsyst::extend(ResponseToGENIEParameters,
  // ConfigureRESWeightEngine(fMetaData, GReWeightEngine);

  return true;
}

std::unique_ptr<EventResponse> GENIEReWeight::GetEventResponse(art::Event &e) {
  std::unique_ptr<EventResponse> er = std::make_unique<EventResponse>();

  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
  art::Handle<std::vector<simb::GTruth>> gTruthHandle;
  e.getByLabel(fGENIEModuleLabel, mcTruthHandle);
  e.getByLabel(fGENIEModuleLabel, gTruthHandle);

  size_t NEventUnits = mcTruthHandle->size();
  if (mcTruthHandle->size() != gTruthHandle->size()) {
    std::cout << "[WARN]: Found " << mcTruthHandle->size()
              << " MC truth instances, and " << gTruthHandle->size()
              << " GENIE truth instances in event " << e.event() << std::endl;
    NEventUnits = std::min(mcTruthHandle->size(), gTruthHandle->size());
  }

  std::vector<std::unique_ptr<genie::EventRecord>> gheps;
  for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
    gheps.emplace_back(
        evgb::RetrieveGHEP(mcTruthHandle->at(eu_it), gTruthHandle->at(eu_it)));
    std::cout << "[INFO]: GENIE Interaction: "
              << gheps.back()->Summary()->AsString() << std::endl;
  }

  er->responses.resize(NEventUnits);

  // Calculate response for each dial-set
  for (auto const &resp : ResponseToGENIEParameters) {

    paramId_t respParamId = fMetaData.headers[resp.first].systParamId;
    // Reset all dials
    for (auto const &gpar : GENIEEngineDials) {
      GReWeightEngine->Systematics().Set(gpar, 0);
    }

    if (fMetaData.headers[resp.first].isCorrection) {

      for (auto const &gparpair : resp.second) {
        GReWeightEngine->Systematics().Set(
            gparpair.first,
            fMetaData.headers[gparpair.second].centralParamValue);
      }
      // std::cout << "[INFO]: GDials: " << std::endl;
      // for (auto const &gpar : GENIEEngineDials) {
      //   std::cout << "\t" << genie::rew::GSyst::AsString(gpar) << " = "
      //             << GReWeightEngine->Systematics().Info(gpar)->CurValue
      //             << std::endl;
      // }

      // auto preconfigure = std::chrono::high_resolution_clock::now();
      GReWeightEngine->Reconfigure();
      // auto postconfigure = std::chrono::high_resolution_clock::now();
      // auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
      //     postconfigure - preconfigure);
      // std::cout << "Reconfigure took " << milliseconds.count() << " ms"
      //           << std::endl;
      for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
        er->responses[eu_it][respParamId].push_back(
            GReWeightEngine->CalcWeight(*gheps[eu_it]));
      }

    } else {
      size_t NVariations = fMetaData.headers[resp.first].paramVariations.size();
      for (size_t var_it = 0; var_it < NVariations; ++var_it) {
        // Calculate weights for each.
        for (auto const &gparpair : resp.second) {
          GReWeightEngine->Systematics().Set(
              gparpair.first,
              fMetaData.headers[gparpair.second].paramVariations[var_it]);
        }
        // std::cout << "[INFO]: GDials: " << std::endl;
        // for (auto const &gpar : GENIEEngineDials) {
        //   std::cout << "\t" << genie::rew::GSyst::AsString(gpar) << " = "
        //             << GReWeightEngine->Systematics().Info(gpar)->CurValue
        //             << std::endl;
        // }
        // auto preconfigure = std::chrono::high_resolution_clock::now();
        GReWeightEngine->Reconfigure();
        // auto postconfigure = std::chrono::high_resolution_clock::now();
        // auto milliseconds =
        //     std::chrono::duration_cast<std::chrono::milliseconds>(
        //         postconfigure - preconfigure);
        // std::cout << "Reconfigure took " << milliseconds.count() << " ms"
        //           << std::endl;

        for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
          // auto preweight = std::chrono::high_resolution_clock::now();
          er->responses[eu_it][respParamId].push_back(
              GReWeightEngine->CalcWeight(*gheps[eu_it]));
          // auto postweight = std::chrono::high_resolution_clock::now();
          // auto milliseconds =
          //     std::chrono::duration_cast<std::chrono::milliseconds>(postweight -
          //                                                           preweight);
          // std::cout << "Reweight[" << eu_it << "](" << e.event() << ") took "
          //           << milliseconds.count() << " ms, got "
          //           << er->responses[eu_it][respParamId].back() << std::endl;
        }
      }
    }
  }
  return er;
}

DEFINE_ART_CLASS_TOOL(GENIEReWeight)
