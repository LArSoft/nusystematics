#include "nusystematics/systproviders/FSILikeEAvailSmearing_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/exceptions.hh"

#ifndef NO_ART
#include "art/Utilities/ToolMacros.h"
#endif

#include "GHEP/GHepParticle.h"

#include "Messenger/Messenger.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

// #define DEBUG_MKSINGLEPI

FSILikeEAvailSmearing::FSILikeEAvailSmearing(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>) {}

#ifndef NO_ART
DEFINE_ART_CLASS_TOOL(FSILikeEAvailSmearing)
#endif

SystMetaData FSILikeEAvailSmearing::BuildSystMetaData(ParameterSet const &cfg,
                                                      paramId_t firstId) {

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseFHiCLSimpleToolConfigurationParameter(cfg, "FSILikeEAvailSmearing",
                                                 phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet templateManifest =
      cfg.get<fhicl::ParameterSet>("FSILikeEAvailSmearing_input_manifest");

  if (!cfg.has_key("FSILikeEAvailSmearing_input_manifest") ||
      !cfg.is_key_to_table("FSILikeEAvailSmearing_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for "
           "FSILikeEAvailSmearing, expected to find a FHiCL table keyed by "
           "FSILikeEAvailSmearing_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }

  size_t NChannels = 0;
  for (std::string const &ch : {"CC", "NC"}) {

    if (!templateManifest.has_key(ch)) {
      continue;
    }

    NChannels++;
  }

  if (!NChannels) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring a FSILikeEAvailSmearing reweighting "
           "instance, failed to find any configured channels. Input templates "
           "must be described by in a table keyed "
           "FSILikeEAvailSmearing_input_manifest with the layout follows that "
           "consumed by "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh";
  }

  tool_options.put("FSILikeEAvailSmearing_input_manifest", templateManifest);

  return smd;
}

namespace {
struct channel_id {
  std::string name;
  FSILikeEAvailSmearing::chan channel;
};
} // namespace

bool FSILikeEAvailSmearing::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  if (!HasParam(GetSystMetaData(), "FSILikeEAvailSmearing")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("FSILikeEAvailSmearing");
  }

  if (!tool_options.has_key("FSILikeEAvailSmearing_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: FSILikeEAvailSmearing parameter exists in the "
           "SystMetaData, but "
           "no FSILikeEAvailSmearing_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "FSILikeEAvailSmearing::BuildSystMetaData, but wasn't, this is a "
           "bug, "
           "please report to the maintiner.";
  }

  fhicl::ParameterSet const &templateManifest =
      tool_options.get<fhicl::ParameterSet>(
          "FSILikeEAvailSmearing_input_manifest");

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "FSILikeEAvailSmearing");

  for (channel_id const &ch :
       std::vector<channel_id>{{"CC", chan::kCC}, {"NC", chan::kNC}}) {

    if (!templateManifest.has_key(ch.name)) {
      continue;
    }

    TemplateHelper th;
    th.Template = std::make_unique<FSILikeEAvailSmearing_ReWeight>();
    th.Template->LoadInputHistograms(
        templateManifest.get<fhicl::ParameterSet>(ch.name));
    th.ZeroIsValid = th.Template->IsValidVariation(0);

    ChannelParameterMapping.emplace(ch.channel, std::move(th));
  }

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");
  genie::Messenger::Instance()->SetPriorityLevel("GHepUtils",
                                                 log4cpp::Priority::FATAL);

  return true;
}

event_unit_response_t
FSILikeEAvailSmearing::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;

  chan evch = ev.Summary()->ProcInfo().IsWeakCC() ? chan::kCC : chan::kNC;

  if (ChannelParameterMapping.find(evch) == ChannelParameterMapping.end()) {
    return resp;
  }

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  std::array<double, 3> kinematics;
  kinematics[0] = emTransfer.Vect().Mag();
  kinematics[1] = emTransfer[0];
  kinematics[2] = GetErecoil_MINERvA_LowRecoil(ev) / emTransfer[0];

  resp.push_back({hdr.systParamId, {}});
  for (double val : hdr.paramVariations) {

    if ((val == 0) && !ChannelParameterMapping[evch].ZeroIsValid) {
      resp.back().responses.push_back(1);
    } else {
      resp.back().responses.push_back(
          ChannelParameterMapping[evch].Template->GetVariation(val,
                                                               kinematics));
    }
  }
  return resp;
}

std::string FSILikeEAvailSmearing::AsString() { return ""; }

FSILikeEAvailSmearing::~FSILikeEAvailSmearing() {}
