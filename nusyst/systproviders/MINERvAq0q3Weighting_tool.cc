#include "MINERvAq0q3Weighting_tool.hh"

#ifndef NO_ART
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#endif

#include "nusyst/responsecalculators/MINERvA2p2hq0q3.hh"

#include "nusyst/utility/exceptions.hh"

#include "TLorentzVector.h"

using namespace larsyst;
using namespace nusyst;
using namespace fhicl;

MINERvAq0q3Weighting::MINERvAq0q3Weighting(ParameterSet const &params)
    : IGENIESystProvider_tool(params), RPATemplateReweighter(nullptr) {}

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

SystMetaData MINERvAq0q3Weighting::ConfigureFromFHICL(ParameterSet const &ps,
                                                      paramId_t firstId) {

  SystMetaData smd;
  // { // RPA
  //   SystParamHeader resp_mnvRPA;
  //   resp_mnvRPA.prettyName = "MINERvATune_RPA";
  //   resp_mnvRPA.systParamId = firstId++;
  //
  //   resp_mnvRPA.opts.push_back(std::string("InputFileRPA=") +
  //                          ps.get<std::string>("InputFileRPA"));
  //
  //   smd.push_back(resp_mnvRPA);
  // }
  { // 2p2h
    SystParamHeader resp_mnv2p2h;
    resp_mnv2p2h.prettyName = "MINERvATune_2p2hGaussEnhancement";
    resp_mnv2p2h.systParamId = firstId++;

    smd.push_back(resp_mnv2p2h);
  }
  { // Gaussian2D
    SystParamHeader resp_Gaus;
    resp_2p2h.prettyName = "2p2hGaussEnhancement_response";
    resp_2p2h.systParamId = firstId++;
    smd.push_back(resp_2p2h);

    SystParamHeader Gaus_Norm;
    Gaus_Norm.prettyName = "2p2hGaussEnhancement_Gaus_Norm";
    Gaus_Norm.isResponselessParam = true;
    Gaus_Norm.systParamId = firstId++;
    Gaus_Norm.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_Norm);

    SystParamHeader Gaus_MeanQ0;
    Gaus_MeanQ0.prettyName = "2p2hGaussEnhancement_Gaus_MeanQ0";
    Gaus_MeanQ0.isResponselessParam = true;
    Gaus_MeanQ0.systParamId = firstId++;
    Gaus_MeanQ0.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_MeanQ0);

    SystParamHeader Gaus_MeanQ3;
    Gaus_MeanQ3.prettyName = "2p2hGaussEnhancement_Gaus_MeanQ3";
    Gaus_MeanQ3.isResponselessParam = true;
    Gaus_MeanQ3.systParamId = firstId++;
    Gaus_MeanQ3.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_MeanQ3);

    SystParamHeader Gaus_SigmaQ0;
    Gaus_SigmaQ0.prettyName = "2p2hGaussEnhancement_Gaus_SigmaQ0";
    Gaus_SigmaQ0.isResponselessParam = true;
    Gaus_SigmaQ0.systParamId = firstId++;
    Gaus_SigmaQ0.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_SigmaQ0);

    SystParamHeader Gaus_SigmaQ3;
    Gaus_SigmaQ3.prettyName = "2p2hGaussEnhancement_Gaus_SigmaQ3";
    Gaus_SigmaQ3.isResponselessParam = true;
    Gaus_SigmaQ3.systParamId = firstId++;
    Gaus_SigmaQ3.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_SigmaQ3);

    SystParamHeader Gaus_Correlation;
    Gaus_Correlation.prettyName = "2p2hGaussEnhancement_Gaus_Correlation";
    Gaus_Correlation.isResponselessParam = true;
    Gaus_Correlation.systParamId = firstId++;
    Gaus_Correlation.responseParamId = resp_2p2h.systParamId;
    smd.push_back(Gaus_Correlation);
  }

  return smd;
}
#endif

bool MINERvAq0q3Weighting::Configure() {

  if (HasParam(fMetaData, "MINERvATune_RPA")) {
    ConfiguredParameters[param_t::kMINERvARPA] =
        GetParamId(fMetaData, "MINERvATune_RPA");

    if (!SystHasOptKV(fMetaData, ConfiguredParameters[param_t::kMINERvARPA],
                      "InputFileRPA")) {
      throw incorrectly_configured()
          << "[ERROR]: Expected to find option "
          << std::quoted("InputFileRPA=/path/to/input/manifest.fcl");
    }

    std::string inputs_manifest = SystGetOptKV(
        fMetaData, ConfiguredParameters[param_t::kMINERvARPA], "InputFileRPA");

    RPATemplateReweighter = std::make_unique<MINERvARPAq0q3_ReWeight>(
        std::map<std::string, larsyst::paramId_t>{
            {{"MINERvATune_RPA", ConfiguredParameters[param_t::kMINERvARPA]}}},
        inputs_manifest);
  }

  if (HasParam(fMetaData, "MINERvATune_2p2hGaussEnhancement")) {
    ConfiguredParameters[param_t::kMINERvA2p2h] =
        GetParamId(fMetaData, "MINERvATune_2p2hGaussEnhancement");
  }

  // if (HasParam(fMetaData, "2p2hGaussEnhancement_response")) {
  //   ConfiguredParameters[param_t::kGaussResponse] =
  //       GetParamId(fMetaData, "2p2hGaussEnhancement_response");
  //
  //   for (param_t_name const &param : std::vector<param_t_name>{
  //            {{"2p2hGaussEnhancement_Gaus_Norm", param_t::kGaussNorm},
  //             {"2p2hGaussEnhancement_Gaus_MeanQ0", param_t::kGaussMeanQ0},
  //             {"2p2hGaussEnhancement_Gaus_MeanQ3", param_t::kGaussMeanQ3},
  //
  //             {"2p2hGaussEnhancement_Gaus_SigmaQ0", param_t::kGaussSigmaQ0},
  //             {"2p2hGaussEnhancement_Gaus_SigmaQ3", param_t::kGaussSigmaQ3},
  //             {"2p2hGaussEnhancement_Gaus_Correlation",
  //              param_t::kGaussCorrelation}}}) {
  //
  //     if (HasParam(fMetaData, param.name)) {
  //       ConfiguredParameters[param.lid] = GetParamId(fMetaData, param.name);
  //     }
  //   }
  // }
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
  return Gaussian2D(q0, q3, *GaussParams);
}

event_unit_response_t
MINERvAq0q3Weighting::GetEventResponse(genie::EventRecord &ev) {

  // make default response for configured parameter
  event_unit_response_t resp;
  for (auto &param : ConfiguredParameters) {
    if (GetParam(fMetaData, param.second).isResponselessParam) {
      continue;
    }
    resp[param.second].responses = std::vector<double>{1};
  }

  if (!ev.Summary()->ProcInfo().IsWeakCC() ||
      (!ev.Summary()->ProcInfo().IsQuasiElastic() &&
       !ev.Summary()->ProcInfo().IsMEC())) {
    return resp;
  }

  TLorentzVector FSLepP4 = ev.Summary()->Kine().FSLeptonP4();
  TLorentzVector ISLepP4 = *ev.Summary()->InitState().GetProbeP4(genie::kRfLab);
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  std::array<double, 2> q0q3{{emTransfer.E(), emTransfer.Vect().Mag()}};

  if (ConfiguredParameters.find(param_t::kMINERvARPA) !=
      ConfiguredParameters.end()) {

    SystParamHeader hdr =
        GetParam(fMetaData, ConfiguredParameters[param_t::kMINERvARPA]);

    resp[hdr.systParamId].responses.clear();
    if (hdr.isCorrection) {
      resp[hdr.systParamId].responses.push_back(
          GetMINERvARPATuneWeight(hdr.centralParamValue, q0q3[0], q0q3[1]));
    } else {
      for (double var : hdr.paramVariations) {
        resp[hdr.systParamId].responses.push_back(
            GetMINERvARPATuneWeight(var, q0q3[0], q0q3[1]));
      }
    }
  }

  if (ConfiguredParameters.find(param_t::kMINERvA2p2h) !=
      ConfiguredParameters.end()) {

    SystParamHeader hdr =
        GetParam(fMetaData, ConfiguredParameters[param_t::kMINERvA2p2h]);

    resp[hdr.systParamId].clear();
    if (hdr.isCorrection) {
      resp[hdr.systParamId].push_back(GetMINERvA2p2hTuneWeight(
          hdr.centralParamValue, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
    } else {
      for (double var : hdr.paramVariations) {
        resp[hdr.systParamId].push_back(GetMINERvA2p2hTuneWeight(
            var, q0q3[0], q0q3[1], GetQELikeTarget(ev)));
      }
    }
  }

  // if (ConfiguredParameters.find(param_t::kGaussResponse) !=
  //     ConfiguredParameters.end()) {
  //
  //   SystParamHeader hdr =
  //       GetParam(fMetaData, ConfiguredParameters[param_t::kGaussResponse]);
  //
  //   resp[hdr.systParamId].clear();
  //   for (double var : hdr.paramVariations) {
  //     resp[hdr.systParamId].push_back(GetMINERvATune_RPAWeight(ev, var));
  //   }
  //
  //     genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  //
  //     // get the targetNucleon variable (dunno the default GENIE variable
  //     right now)
  //     //
  //     // subtract 2000000200
  //     // if ==0 or ==2 it is nn or pp
  //     // if ==1 it is np
  // }

  return resp;
}

std::string MINERvAq0q3Weighting::AsString() { return ""; }
