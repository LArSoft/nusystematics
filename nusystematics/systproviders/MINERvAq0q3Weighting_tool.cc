#include "nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh"

#include "nusystematics/responsecalculators/MINERvA2p2hq0q3.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#ifndef NO_ART
#include "art/Utilities/ToolMacros.h"
#endif

#include "GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

MINERvAq0q3Weighting::MINERvAq0q3Weighting(ParameterSet const &params)
    : IGENIESystProvider_tool(params), RPATemplateReweighter(nullptr),
      valid_file(nullptr), valid_tree(nullptr) {}

#ifndef NO_ART
DEFINE_ART_CLASS_TOOL(MINERvAq0q3Weighting)
#endif

SystMetaData MINERvAq0q3Weighting::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  SystMetaData smd;
  if (cfg.get<bool>("use_MINERvA_RPA_tunes", false)) { // RPA
    systtools::SystParamHeader param;
    param.systParamId = firstId++;

    param.prettyName = "MINERvATune_RPA";
    param.centralParamValue = 0;
    param.paramVariations = std::vector<double>{-1, 0, 1};

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

    smd.push_back(param);
  }

  if (cfg.get<bool>("use_MINERvA_2p2h_tunes", false)) { // 2p2h

    parameter_per_2p2h_universe =
        cfg.get<bool>("parameter_per_2p2h_universe", false);

    tool_options.put("parameter_per_2p2h_universe",
                     parameter_per_2p2h_universe);

    if (parameter_per_2p2h_universe) {
      systtools::SystParamHeader param_CV;
      param_CV.systParamId = firstId++;
      param_CV.isCorrection = true;
      param_CV.centralParamValue = 1;
      param_CV.prettyName = "MINERvATune_2p2hGaussEnhancement_CV";
      smd.push_back(param_CV);
      systtools::SystParamHeader param_NN;
      param_NN.systParamId = firstId++;
      param_NN.isCorrection = true;
      param_NN.centralParamValue = 1;
      param_NN.prettyName = "MINERvATune_2p2hGaussEnhancement_NN";
      smd.push_back(param_NN);
      systtools::SystParamHeader param_np;
      param_np.systParamId = firstId++;
      param_np.isCorrection = true;
      param_np.centralParamValue = 1;
      param_np.prettyName = "MINERvATune_2p2hGaussEnhancement_np";
      smd.push_back(param_np);
      systtools::SystParamHeader param_QE;
      param_QE.systParamId = firstId++;
      param_QE.isCorrection = true;
      param_QE.centralParamValue = 1;
      param_QE.prettyName = "MINERvATune_2p2hGaussEnhancement_QE";
      smd.push_back(param_QE);
    } else {
      systtools::SystParamHeader param;
      param.systParamId = firstId++;
      param.centralParamValue = 1;
      param.paramVariations = {1, 2, 3, 4};
      param.prettyName = "MINERvATune_2p2hGaussEnhancement";
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
        GetParamIndex(GetSystMetaData(), "MINERvATune_RPA");

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
        tool_options.get<fhicl::ParameterSet>(
            "MINERvATune_RPA_input_manifest"));
  }

  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement")) {
    ConfiguredParameters[param_t::kMINERvA2p2h] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement");
  }
  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_CV")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_CV] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_CV");
  }
  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_NN")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_NN] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_NN");
  }
  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_np")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_np] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_np");
  }
  if (HasParam(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_QE")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_QE] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_2p2hGaussEnhancement_QE");
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
  if (val == 1) {
    if (!((QELTarget == QELikeTarget_t::kNN) ||
          (QELTarget == QELikeTarget_t::knp))) {
      return 1;
    }

    GaussParams = &Gauss2DParams_CV;
  } else if (val == 2) {

    if (QELTarget != QELikeTarget_t::kNN) {
      return 1;
    }

    GaussParams = &Gauss2DParams_NNOnly;
  } else if (val == 3) {

    if (QELTarget != QELikeTarget_t::knp) {
      return 1;
    }

    GaussParams = &Gauss2DParams_npOnly;
  } else if (val == 4) {

    if (QELTarget != QELikeTarget_t::kQE) {
      return 1;
    }

    GaussParams = &Gauss2DParams_1p1hOnly;
  } else {
    throw invalid_parameter_value()
        << "[ERROR]: When applying MINERvA 2p2h tune expected to find "
           "parameter values of [ 1 == CV, 2 == NN, 3 == np, 4 == 1p1h ], "
           "but found "
        << val;
  }
  return 1 + Gaussian2D(q0, q3, *GaussParams);
}

event_unit_response_t
MINERvAq0q3Weighting::GetEventResponse(genie::EventRecord const &ev) {

  // make default response for configured parameter
  event_unit_response_t resp;

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  if (!(ev.Summary()->ProcInfo().IsQuasiElastic() ||
        ev.Summary()->ProcInfo().IsMEC())) {
    return resp;
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

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvARPA]];

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

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h]];

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

  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_CV) !=
       ConfiguredParameters.end()) &&
      ev.Summary()->ProcInfo().IsMEC()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_CV]];

    resp.push_back({hdr.systParamId, {}});
    resp.back().responses.push_back(
        GetMINERvA2p2hTuneWeight(1, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
  }
  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_NN) !=
       ConfiguredParameters.end()) &&
      ev.Summary()->ProcInfo().IsMEC()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_NN]];

    resp.push_back({hdr.systParamId, {}});
    resp.back().responses.push_back(
        GetMINERvA2p2hTuneWeight(2, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
  }
  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_np) !=
       ConfiguredParameters.end()) &&
      ev.Summary()->ProcInfo().IsMEC()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_np]];

    resp.push_back({hdr.systParamId, {}});
    resp.back().responses.push_back(
        GetMINERvA2p2hTuneWeight(3, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
  }
  // Only ever applies to qe events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_QE) !=
       ConfiguredParameters.end()) &&
      ev.Summary()->ProcInfo().IsQuasiElastic()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_QE]];

    resp.push_back({hdr.systParamId, {}});
    resp.back().responses.push_back(
        GetMINERvA2p2hTuneWeight(4, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
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
    MEC_weights.clear();

    if (ConfiguredParameters.find(param_t::kMINERvARPA) !=
        ConfiguredParameters.end()) {
      paramId_t RPA_param =
          GetSystMetaData()[ConfiguredParameters[param_t::kMINERvARPA]]
              .systParamId;
      RPA_weights = GetParamElementFromContainer(resp, RPA_param).responses;
    }
    if ((ConfiguredParameters.find(param_t::kMINERvA2p2h) !=
         ConfiguredParameters.end()) &&
        (ev.Summary()->ProcInfo().IsQuasiElastic() ||
         ev.Summary()->ProcInfo().IsMEC())) {
      paramId_t MEC_param =
          GetSystMetaData()[ConfiguredParameters[param_t::kMINERvARPA]]
              .systParamId;
      MEC_weights = GetParamElementFromContainer(resp, MEC_param).responses;
    }
    for (param_t tune_2p2h_universe :
         {param_t::kMINERvA2p2h_CV, param_t::kMINERvA2p2h_NN,
          param_t::kMINERvA2p2h_np, param_t::kMINERvA2p2h_QE}) {
      if ((ConfiguredParameters.find(tune_2p2h_universe) !=
           ConfiguredParameters.end()) &&
          (ev.Summary()->ProcInfo().IsQuasiElastic() ||
           ev.Summary()->ProcInfo().IsMEC())) {
        paramId_t MEC_param =
            GetSystMetaData()[ConfiguredParameters[tune_2p2h_universe]]
                .systParamId;
        size_t idx = GetParamContainerIndex(resp, MEC_param);
        if (idx != kParamUnhandled<size_t>) {
          MEC_weights.push_back(resp[idx].responses.front());
        }
      }
    }

    nRPA_weights = RPA_weights.size();
    nMEC_weights = MEC_weights.size();
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
  valid_tree->Branch("nRPA_weights", &nRPA_weights);
  valid_tree->Branch("nMEC_weights", &nMEC_weights);
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
