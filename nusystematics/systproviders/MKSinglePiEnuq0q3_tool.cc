#include "MKSinglePiEnuq0q3_tool.hh"

#ifndef NO_ART
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#endif

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/exceptions.hh"

using namespace systtools;
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
} // namespace

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
#endif

SystMetaData MKSinglePiEnuq0q3::BuildSystMetaData(ParameterSet const &cfg,
                                                  paramId_t firstId) {

  SystMetaData smd;

  SystParamHeader resp;
  resp.prettyName = "MKSPP_Enuq0q3_response";
  resp.systParamId = firstId++;

  resp.centralParamValue = 1;
  resp.isCorrection = true;
  smd.push_back(resp);

  if (!cfg.has_key("MKSPP_Enuq0q3_input_manifest") ||
      !cfg.is_key_to_table("MKSPP_Enuq0q3_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for "
           "MKSPP_Enuq0q3, expected to find a FHiCL table keyed by "
           "MKSPP_Enuq0q3_input_manifest describing the location of the "
           "histogram inputs. See "
           "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh "
           "for the layout.";
  }

  fhicl::ParameterSet ps =
      cfg.get<fhicl::ParameterSet>("MKSPP_Enuq0q3_input_manifest");
  tool_options.put("MKSPP_Enuq0q3_input_manifest", ps);
  size_t NChannels = 0;
  for (fhicl::ParameterSet const &ch_ps :
       ps.get<std::vector<fhicl::ParameterSet>>("InputTemplates")) {

    SystParamHeader channel;
    channel.prettyName = std::string("MKSPP_Enuq0q3_") +
                         ch_ps.get<std::string>("parameter_name");
    channel.systParamId = firstId++;

    channel.centralParamValue = 1;
    channel.isCorrection = true;

    channel.isResponselessParam = true;
    channel.responseParamId = resp.systParamId;
    smd.push_back(channel);
    NChannels++;
  }

  if (!NChannels) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring a MKSPP_Enuq0q3 reweighting instance, "
           "failed to find any configured channels. Input templates must be "
           "described by in a table keyed MKSPP_Enuq0q3_input_manifest with "
           "the layout follows that consumed by "
           "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh";
  }

  return smd;
}

bool MKSinglePiEnuq0q3::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  if (!HasParam(GetSystMetaData(), "MKSPP_Enuq0q3_response")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("MKSPP_Enuq0q3_response");
  }

  ResponseParameterId = GetParamId(GetSystMetaData(), "MKSPP_Enuq0q3_response");

  SystParamHeader hdr = GetParam(GetSystMetaData(), ResponseParameterId);

  std::map<std::string, systtools::paramId_t> ParamNames;
  for (channel_id const &ch :
       std::vector<channel_id>{{{"NumuPPiPlus", genie::kSpp_vp_cc_10100},
                                {"NumuPPi0", genie::kSpp_vn_cc_10010},
                                {"NumuNPiPlus", genie::kSpp_vn_cc_01100},

                                {"NumuBNPiMinus", genie::kSpp_vbn_cc_01001},
                                {"NumuBNPi0", genie::kSpp_vbp_cc_01010},
                                {"NumuBPPiMinus", genie::kSpp_vbp_cc_10001}}}) {

    if (HasParam(GetSystMetaData(), std::string("MKSPP_Enuq0q3_") + ch.name)) {
      systtools::paramId_t pid = GetParamId(
          GetSystMetaData(), std::string("MKSPP_Enuq0q3_") + ch.name);
      ParamNames[ch.name] = pid;
      ChannelParameterMapping[ch.channel] = pid;
    }
  }

  if (!tool_options.has_key("MKSPP_Enuq0q3_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: MKSPP_Enuq0q3_response parameter exists in the "
           "SystMetaData, "
           "but no MKSPP_Enuq0q3_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "MKSinglePiEnuq0q3::BuildSystMetaData, but wasn't, this is a "
           "bug, please report to the maintiner.";
  }

  templateReweighter = std::make_unique<MKSinglePiEnuq0q3_ReWeight>(
      ParamNames,
      tool_options.get<fhicl::ParameterSet>("MKSPP_Enuq0q3_input_manifest"));
}

event_unit_response_t
MKSinglePiEnuq0q3::GetEventResponse(genie::EventRecord &ev) {

  event_unit_response_t resp{{{ResponseParameterId, std::vector<double>{{1}}}}};

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

  resp.push_back({ResponseParameterId,
                  std::vector<double>{{templateReweighter->GetVariation(
                      ChannelParameterMapping[chan], 1,
                      std::array<double, 3>{{Enu, q0, q3}})}}});

  return resp;
}

std::string MKSinglePiEnuq0q3::AsString() { return ""; }
