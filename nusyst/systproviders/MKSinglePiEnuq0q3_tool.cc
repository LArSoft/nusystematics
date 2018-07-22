#include "MKSinglePiEnuq0q3_tool.hh"

#ifndef NO_ART
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#endif

#include "nusyst/utility/GENIEUtils.hh"
#include "nusyst/utility/exceptions.hh"

using namespace larsyst;
using namespace nusyst;
using namespace fhicl;

// #define DEBUG_MKSINGLEPI

MKSinglePiEnuq0q3::MKSinglePiEnuq0q3(ParameterSet const &params)
    : IGENIESystProvider_tool(params), templateReweighter(nullptr) {}

namespace {
struct channel_id {
  std::string name;
  genie::SppChannel_t channel;
};
}

#ifndef NO_ART
std::unique_ptr<EventResponse>
MKSinglePiEnuq0q3::GetEventResponse(art::Event &e) {
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
  for (size_t eu_it = 0; eu_it < NEventUnits; ++eu_it) {
    er->responses.push_back(GetEventResponse(*gheps[eu_it]));
  }
  return er;
}
SystMetaData MKSinglePiEnuq0q3::ConfigureFromFHICL(ParameterSet const &ps,
                                                   paramId_t firstId) {

  SystMetaData smd;

  SystParamHeader resp;
  resp.prettyName = "MKSPP_Enuq0q3_response";
  resp.systParamId = firstId++;

  resp.centralParamValue = 1;
  resp.isCorrection = true;

  resp.opts.push_back(std::string("InputManifest=") +
                      ps.get<std::string>("InputManifest"));

  smd.push_back(resp);

  for (std::string const &chName :
       ps.get<std::vector<std::string>>("channels")) {
    SystParamHeader channel;
    channel.prettyName = chName;
    channel.systParamId = firstId++;

    channel.centralParamValue = 1;
    channel.isCorrection = true;

    channel.isResponselessParam = true;
    channel.responseParamId = resp.systParamId;
    smd.push_back(channel);
  }

  return smd;
}
#endif

bool MKSinglePiEnuq0q3::Configure() {

  ResponseParameterId = GetParamId(fMetaData, "MKSPP_Enuq0q3_response");

  if (ResponseParameterId == kParamUnhandled<paramId_t>) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("MKSPP_Enuq0q3_response");
  }

  SystParamHeader hdr = GetParam(fMetaData, ResponseParameterId);

  if (!SystHasOptKV(fMetaData, ResponseParameterId, "InputManifest")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find option "
        << std::quoted("InputManifest=/path/to/input/manifest.fcl");
  }

  std::string inputs_manifest =
      SystGetOptKV(fMetaData, ResponseParameterId, "InputManifest");

  std::map<std::string, larsyst::paramId_t> ParamNames;
  for (channel_id const &ch :
       std::vector<channel_id>{{{"NumuPPiPlus", genie::kSpp_vp_cc_10100},
                                {"NumuPPi0", genie::kSpp_vn_cc_10010},
                                {"NumuNPiPlus", genie::kSpp_vn_cc_01100},

                                {"NumuBNPiMinus", genie::kSpp_vbn_cc_01001},
                                {"NumuBNPi0", genie::kSpp_vbp_cc_01010},
                                {"NumuBPPiMinus", genie::kSpp_vbp_cc_10001}}}) {

    if (HasParam(fMetaData, ch.name)) {
      larsyst::paramId_t pid = GetParamId(fMetaData, ch.name);
      ParamNames[ch.name] = pid;
      ChannelParameterMapping[ch.channel] = pid;
    }
  }

  templateReweighter =
      std::make_unique<MKSinglePiEnuq0q3_ReWeight>(ParamNames, inputs_manifest);
}

event_unit_response_t
MKSinglePiEnuq0q3::GetEventResponse(genie::EventRecord &ev) {

  event_unit_response_t resp{{ResponseParameterId, std::vector<double>{{1}}}};

  if (!ev.Summary()->ProcInfo().IsResonant() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  genie::SppChannel_t chan = SPPChannelFromGHep(ev);

  if ((chan == genie::kSppNull) ||
      (ChannelParameterMapping.find(chan) == ChannelParameterMapping.end())) {
    return resp;
  }

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  if (!tgt.HitNucIsSet()) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to get hit nucleon kinematics as it was not "
           "included in this GHep event. This is a fatal error.";
  }

  TLorentzVector NucP4 = tgt.HitNucP4();
  TLorentzVector FSLepP4 = ev.Summary()->Kine().FSLeptonP4();
  TLorentzVector ISLepP4 = *ev.Summary()->InitState().GetProbeP4(genie::kRfLab);

#ifdef DEBUG_MKSINGLEPI
  std::cout << "[INFO]: Lab frame ENu: " << ISLepP4.E()
            << ", q0: " << (ISLepP4 - FSLepP4).E()
            << ", q3: " << (ISLepP4 - FSLepP4).Vect().Mag()
            << ", Target Nucleon: [ " << NucP4[0] << ", " << NucP4[1] << ", "
            << NucP4[2] << ", M: " << NucP4.M() << "]" << std::endl;
#endif

  FSLepP4.Boost(-NucP4.BoostVector());
  ISLepP4.Boost(-NucP4.BoostVector());
  NucP4.Boost(-NucP4.BoostVector());

#ifdef DEBUG_MKSINGLEPI
  std::cout << "[INFO]: Post-boost: Target Nucleon: [ " << NucP4[0] << ", "
            << NucP4[1] << ", " << NucP4[2] << ", M: " << NucP4.M() << "]"
            << std::endl;
#endif

  double Enu = ISLepP4.E();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
  double q0 = emTransfer.E();
  double q3 = emTransfer.Vect().Mag();

  double weight = 1;

  resp[ResponseParameterId].responses[0] = templateReweighter->GetVariation(
      ChannelParameterMapping[chan], 1, std::array<double, 3>{{Enu, q0, q3}});

  return resp;
}

std::string MKSinglePiEnuq0q3::AsString() { return ""; }
