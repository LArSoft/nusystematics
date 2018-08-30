#include "MKSinglePiEnuq0q3_tool.hh"

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

MKSinglePiEnuq0q3::MKSinglePiEnuq0q3(ParameterSet const &params)
    : IGENIESystProvider_tool(params), templateReweighter(nullptr),
      valid_file(nullptr), valid_tree(nullptr) {}

namespace {
struct channel_id {
  std::string name;
  genie::SppChannel_t channel;
};
} // namespace

#ifndef NO_ART
DEFINE_ART_CLASS_TOOL(MKSinglePiEnuq0q3)
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
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh "
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
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh";
  }

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool MKSinglePiEnuq0q3::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

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

  fill_valid_tree = tool_options.get("fill_valid_tree", false);
  std::cout << tool_options.to_indented_string() << std::endl;
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
MKSinglePiEnuq0q3::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp{{{ResponseParameterId, std::vector<double>{{1}}}}};

  if (!(ev.Summary()->ProcInfo().IsResonant() ||
        ev.Summary()->ProcInfo().IsDeepInelastic()) ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  if (ev.Summary()->Kine().W(true) > 1.7) {
    return resp;
  }

  genie::SppChannel_t chan = genie::kSppNull;

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  if (!tgt.HitNucIsSet()) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to get hit nucleon kinematics as it was not "
           "included in this GHep event. This is a fatal error.";
  }

  TLorentzVector NucP4 = tgt.HitNucP4();

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();

  FSLepP4.Boost(-NucP4.BoostVector());
  ISLepP4.Boost(-NucP4.BoostVector());
  NucP4.Boost(-NucP4.BoostVector());

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double wght = 1;

  if (ev.Summary()->ProcInfo().IsResonant()) {
    chan = SPPChannelFromGHep(ev);

    if ((chan == genie::kSppNull) ||
        (ChannelParameterMapping.find(chan) == ChannelParameterMapping.end())) {
      return resp;
    }

#ifdef DEBUG_MKSINGLEPI
    int NEUTCh = genie::utils::ghep::NeutReactionCode(&ev);

    int NPi0 = 0, NPip = 0, NP = 0, NN = 0;

    bool nuclear_target = ev.Summary()->InitState().Tgt().IsNucleus();

    TIter event_iter(&ev);
    genie::GHepParticle *p = 0;

    while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
      genie::GHepStatus_t ghep_ist = (genie::GHepStatus_t)p->Status();
      int ghep_pdgc = p->Pdg();
      int ghep_fm = p->FirstMother();
      int ghep_fmpdgc = (ghep_fm == -1) ? 0 : ev.Particle(ghep_fm)->Pdg();

      // For nuclear targets use hadrons marked as 'hadron in the nucleus'
      // which are the ones passed in the intranuclear rescattering
      // For free nucleon targets use particles marked as 'final state'
      // but make an exception for decayed pi0's,eta's (count them and not their
      // daughters)

      bool decayed =
          (ghep_ist == genie::kIStDecayedState &&
           (ghep_pdgc == genie::kPdgPi0 || ghep_pdgc == genie::kPdgEta));
      bool parent_included =
          (ghep_fmpdgc == genie::kPdgPi0 || ghep_fmpdgc == genie::kPdgEta);

      bool count_it =
          (nuclear_target && ghep_ist == genie::kIStHadronInTheNucleus) ||
          (!nuclear_target && decayed) ||
          (!nuclear_target && ghep_ist == genie::kIStStableFinalState &&
           !parent_included);

      if (!count_it)
        continue;
      if (ghep_pdgc == genie::kPdgPiP) {
        NPip++; // pi+
      } else if (ghep_pdgc == genie::kPdgPi0) {
        NPi0++; // pi0NPim = 0
      } else if (ghep_pdgc == genie::kPdgProton) {
        NP++; // pi0NPim = 0
      } else if (ghep_pdgc == genie::kPdgNeutron) {
        NN++; // pi0NPim = 0
      }
    }

    std::cout << "[INFO]: Event, NEUTChannel: " << NEUTCh
              << ", GENIE SPPChannel: " << chan
              << ", W = " << ev.Summary()->Kine().W(true)
              << ", PID: " << ChannelParameterMapping[chan]
              << ", NPi0 = " << NPi0 << ", NPip = " << NPip << ", NP = " << NP
              << ", NN = " << NN << std::endl;
#endif

#ifdef DEBUG_MKSINGLEPI
    std::cout << "[INFO]: Lab frame ENu: " << ISLepP4.E()
              << ", q0: " << (ISLepP4 - FSLepP4).E()
              << ", q3: " << (ISLepP4 - FSLepP4).Vect().Mag()
              << ", Target Nucleon: [ " << NucP4[0] << ", " << NucP4[1] << ", "
              << NucP4[2] << ", M: " << NucP4.M() << "]" << std::endl;
#endif

#ifdef DEBUG_MKSINGLEPI
    std::cout << "[INFO]: Post-boost: Target Nucleon: [ " << NucP4[0] << ", "
              << NucP4[1] << ", " << NucP4[2] << ", M: " << NucP4.M() << "]"
              << std::endl;
#endif

    wght = templateReweighter->GetVariation(
        ChannelParameterMapping[chan], 1,
        std::array<double, 3>{
            {ISLepP4.E(), emTransfer.E(), emTransfer.Vect().Mag()}});

  } else { // GENIE non-resonant gets weighted to 0.
    wght = 0;
  }

  resp.push_back({ResponseParameterId, std::vector<double>{{wght}}});

  if (fill_valid_tree) {

    q0_nuc_rest_frame = emTransfer.E();
    q3_nuc_rest_frame = emTransfer.Vect().Mag();

    ISLepP4 = *ISLep->P4();
    FSLepP4 = *FSLep->P4();
    emTransfer = (ISLepP4 - FSLepP4);

    pdgfslep = ev.FinalStatePrimaryLepton()->Pdg();
    momfslep = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();

    Pdgnu = ISLep->Pdg();
    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

    SppChannel = chan;

    momhmfspi = 0;
    size_t NParts = ev.GetEntries();
    for (size_t it = 0; it < NParts; ++it) {
      genie::GHepParticle *part = ev.Particle(it);
      if (part->Status() != genie::kIStStableFinalState) {
        continue;
      }
      if ((part->Pdg() != 211) && (part->Pdg() != 111)) {
        continue;
      }
      if (momhmfspi < part->P4()->Vect().Mag()) {
        pdghmfspi = part->Pdg();
        momhmfspi = part->P4()->Vect().Mag();
        cthetahmfspi = part->P4()->Vect().CosTheta();
      }
    }

    weight = wght;

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    OutOfReWeightPS = (templateReweighter->GetBin(
                           ChannelParameterMapping[chan],
                           std::array<double, 3>{{ISLepP4.E(), emTransfer.E(),
                                                  emTransfer.Vect().Mag()}}) ==
                       templateReweighter->kBinOutsideRange);

    valid_tree->Fill();
  }

  return resp;
}

std::string MKSinglePiEnuq0q3::AsString() { return ""; }

void MKSinglePiEnuq0q3::InitValidTree() {
  valid_file = new TFile("MKSPPEnuq0q3WeightValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("SppChannel", &SppChannel);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Pdg_HMFSPi", &pdghmfspi);
  valid_tree->Branch("P_HMFSPi", &momhmfspi);
  valid_tree->Branch("CosTheta_HMFSPi", &cthetahmfspi);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("q0_nuc_rest_frame", &q0_nuc_rest_frame);
  valid_tree->Branch("q3_nuc_rest_frame", &q3_nuc_rest_frame);
  valid_tree->Branch("weight", &weight);
  valid_tree->Branch("OutOfReWeightPS", &OutOfReWeightPS);
}

MKSinglePiEnuq0q3::~MKSinglePiEnuq0q3() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
