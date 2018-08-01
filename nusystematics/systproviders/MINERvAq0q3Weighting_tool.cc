#include "MINERvAq0q3Weighting_tool.hh"

#ifndef NO_ART
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#endif

#include "nusystematics/responsecalculators/MINERvA2p2hq0q3.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

MINERvAq0q3Weighting::MINERvAq0q3Weighting(ParameterSet const &params)
    : IGENIESystProvider_tool(params), RPATemplateReweighter(nullptr),
      valid_file(nullptr), valid_tree(nullptr), ApplyRPAToSPP(false),
      ApplyRPAToRES(false) {}

#ifndef NO_ART
std::unique_ptr<EventResponse>
MINERvAq0q3Weighting::GetEventResponse(art::Event &e) {
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

SystMetaData MINERvAq0q3Weighting::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  SystMetaData smd;
  { // RPA
    systtools::SystParamHeader param;
    if (ParseFHiCLSimpleToolConfigurationParameter(cfg, "MINERvATune_RPA",
                                                   param)) {
      param.systParamId = firstId++;

      if (param.isRandomlyThrown || param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: When configuring calculated variations for "
               "MINERvATune_RPA, found MINERvATune_RPA_variation_descriptor "
               "that described a spline or a set of random throws. This "
               "parameter can only be set to: -1, 0, 1 for three different "
               "tunes of the RPA response.";
      }

      if (!cfg.has_key("MINERvATune_RPA_input_manifest") ||
          !cfg.is_key_to_table("MINERvATune_RPA_input_manifest")) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: When configuring calculated variations for "
               "MINERvATune_RPA, expected to find a FHiCL table keyed by "
               "MINERvATune_RPA_input_manifest describing the location of the "
               "histogram inputs. See "
               "nusystematics/responsecalculators/"
               "TemplateResponseCalculatorBase.hh "
               "for the layout.";
      }

      fhicl::ParameterSet ps =
          cfg.get<fhicl::ParameterSet>("MINERvATune_RPA_input_manifest");
      tool_options.put("MINERvATune_RPA_input_manifest", ps);

      ApplyRPAToSPP = cfg.get("apply_RPA_to_SPP", false);
      ApplyRPAToRES = cfg.get("apply_RPA_to_RES", false);
      if (ApplyRPAToSPP) {
        tool_options.put("apply_RPA_to_SPP", true);
      }
      if (ApplyRPAToRES) {
        tool_options.put("apply_RPA_to_RES", true);
      }

      smd.push_back(param);
    }
  }
  { // 2p2h
    systtools::SystParamHeader param;
    if (ParseFHiCLSimpleToolConfigurationParameter(
            cfg, "MINERvATune_2p2hGaussEnhancement", param)) {
      param.systParamId = firstId++;

      if (param.isRandomlyThrown || param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: When configuring calculated variations for "
               "MINERvATune_2p2hGaussEnhancement, found "
               "MINERvATune_2p2hGaussEnhancement_variation_descriptor that "
               "described a spline or a set of random throws. This parameter "
               "can only be set to: 0, 1, 2, 3 for four different tunes of the "
               "2p2h enhancement.";
      }

      smd.push_back(param);
    }
  }

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool MINERvAq0q3Weighting::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  if (HasParam(GetSystMetaData(), "MINERvATune_RPA")) {
    ConfiguredParameters[param_t::kMINERvARPA] =
        GetParamId(GetSystMetaData(), "MINERvATune_RPA");

    if (!tool_options.has_key("MINERvATune_RPA_input_manifest")) {
      throw systtools::invalid_ToolOptions()
          << "[ERROR]: MINERvATune_RPA parameter exists in the SystMetaData, "
             "but no MINERvATune_RPA_input_manifest key can be found on the "
             "tool_options table. This reweighting requires input histograms "
             "that must be specified. This should have been caught by  "
             "MINERvAq0q3Weighting::BuildSystMetaData, but wasn't, this is a "
             "bug, please report to the maintiner.";
    }

    RPATemplateReweighter = std::make_unique<MINERvARPAq0q3_ReWeight>(
        std::map<std::string, systtools::paramId_t>{
            {{"MINERvATune_RPA", ConfiguredParameters[param_t::kMINERvARPA]}}},
        tool_options.get<fhicl::ParameterSet>(
            "MINERvATune_RPA_input_manifest"));

    ApplyRPAToSPP = tool_options.get("apply_RPA_to_SPP", false);
    ApplyRPAToRES = tool_options.get("apply_RPA_to_RES", false);
  }

  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement")) {
    ConfiguredParameters[param_t::kMINERvA2p2h] =
        GetParamId(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement");
  }

  fill_valid_tree = tool_options.get("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }
  return true;
}

double MINERvAq0q3Weighting::GetMINERvARPATuneWeight(double val, double q0,
                                                     double q3) {
  MINERvARPAq0q3_ReWeight::RPATweak_t tval;
  if (val == 0) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kCV;
  } else if (val == 1) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kPlus1;
  } else if (val == -1) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kMinus1;
  } else {
    throw invalid_parameter_value()
        << "[ERROR]: When applying MINERvA RPA tune expected to find "
           "parameter values of [ 0 == CV, 1 == Plus1, -1 == Minus1 ], "
           "but found "
        << val;
  }

  return RPATemplateReweighter->GetWeight(q0, q3, tval);
}

double
MINERvAq0q3Weighting::GetMINERvA2p2hTuneWeight(double val, double q0, double q3,
                                               QELikeTarget_t QELTarget) {
  std::array<double, 6> const *GaussParams;
  if (val == 0) {
    if (!((QELTarget == QELikeTarget_t::kNN) ||
          (QELTarget == QELikeTarget_t::knp))) {
      return 1;
    }

    GaussParams = &Gauss2DParams_CV;
  } else if (val == 1) {

    if (QELTarget != QELikeTarget_t::kNN) {
      return 1;
    }

    GaussParams = &Gauss2DParams_NNOnly;
  } else if (val == 2) {

    if (QELTarget != QELikeTarget_t::knp) {
      return 1;
    }

    GaussParams = &Gauss2DParams_npOnly;
  } else if (val == 3) {

    if (QELTarget != QELikeTarget_t::kQE) {
      return 1;
    }

    GaussParams = &Gauss2DParams_1p1hOnly;
  } else {
    throw invalid_parameter_value()
        << "[ERROR]: When applying MINERvA 2p2h tune expected to find "
           "parameter values of [ 0 == CV, 1 == NN, 2 == np, 3 == 1p1h ], "
           "but found "
        << val;
  }
  return Gaussian2D(q0, q3, *GaussParams) + 1;
}

event_unit_response_t
MINERvAq0q3Weighting::GetEventResponse(genie::EventRecord &ev) {

  // make default response for configured parameter
  event_unit_response_t resp;

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  if (!(ev.Summary()->ProcInfo().IsQuasiElastic() ||
        ev.Summary()->ProcInfo().IsMEC()) &&
      !(ApplyRPAToSPP || ApplyRPAToRES)) {
    return resp;
  } else if (ApplyRPAToSPP || ApplyRPAToRES) {
    bool GoodSPP = false;
    if (ApplyRPAToSPP) {
      genie::SppChannel_t sppch = nusyst::SPPChannelFromGHep(ev);
      GoodSPP = (sppch != genie::kSppNull);
    }
    bool GoodRES = ApplyRPAToRES && !ev.Summary()->ProcInfo().IsResonant();
    if (!(GoodRES || GoodSPP)) {
      return resp;
    }
  }

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
  std::array<double, 2> q0q3{{emTransfer.E(), emTransfer.Vect().Mag()}};

  // Can apply to RES/SPP events
  if (ConfiguredParameters.find(param_t::kMINERvARPA) !=
      ConfiguredParameters.end()) {

    SystParamHeader hdr =
        GetParam(GetSystMetaData(), ConfiguredParameters[param_t::kMINERvARPA]);

    resp.push_back({hdr.systParamId, {}});
    if (hdr.isCorrection) {
      resp.back().responses.push_back(
          GetMINERvARPATuneWeight(hdr.centralParamValue, q0q3[0], q0q3[1]));
    } else {
      for (double var : hdr.paramVariations) {
        resp.back().responses.push_back(
            GetMINERvARPATuneWeight(var, q0q3[0], q0q3[1]));
      }
    }
  }

  // Only ever applies to 2p2h/qe events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h) !=
       ConfiguredParameters.end()) &&
      (ev.Summary()->ProcInfo().IsQuasiElastic() ||
       ev.Summary()->ProcInfo().IsMEC())) {

    SystParamHeader hdr = GetParam(GetSystMetaData(),
                                   ConfiguredParameters[param_t::kMINERvA2p2h]);

    resp.push_back({hdr.systParamId, {}});
    if (hdr.isCorrection) {
      resp.back().responses.push_back(GetMINERvA2p2hTuneWeight(
          hdr.centralParamValue, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
    } else {
      for (double var : hdr.paramVariations) {
        resp.back().responses.push_back(GetMINERvA2p2hTuneWeight(
            var, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
      }
    }
  }

  if (fill_valid_tree) {

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

    QELTarget = e2i(GetQELikeTarget(ev));

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    RPA_weights.clear();
    paramId_t RPA_param = GetParamId(GetSystMetaData(), "MINERvATune_RPA");
    if (RPA_param != kParamUnhandled<paramId_t>) {
      RPA_weights = GetParamElementFromContainer(resp, RPA_param).responses;
    }

    MEC_weights.clear();
    paramId_t MEC_param =
        GetParamId(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement");
    if (MEC_param != kParamUnhandled<paramId_t>) {
      MEC_weights = GetParamElementFromContainer(resp, MEC_param).responses;
    }
    valid_tree->Fill();
  }

  return resp;
}

std::string MINERvAq0q3Weighting::AsString() { return ""; }

void MINERvAq0q3Weighting::InitValidTree() {

  valid_file = new TFile("MINERvAq0q3WeightValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("QELTarget", &QELTarget);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("RPA_weights", &RPA_weights);
  valid_tree->Branch("MEC_weights", &MEC_weights);
}

MINERvAq0q3Weighting::~MINERvAq0q3Weighting() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
