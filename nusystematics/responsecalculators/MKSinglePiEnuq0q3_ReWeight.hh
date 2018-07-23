#ifndef nusystematics_RESPONSE_CALCULATORS_MKSINGLEPIENUQ0Q3_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_MKSINGLEPIENUQ0Q3_REWEIGHT_HH_SEEN

#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

#include "systematicstools/utility/ROOTUtility.hh"

#include "fhiclcpp/make_ParameterSet.h"

NEW_SYSTTOOLS_EXCEPT(is_flow_bin);

// #define MKSinglePiEnuq0q3_ReWeight_DEBUG

namespace nusyst {
class MKSinglePiEnuq0q3_ReWeight
    : public nusyst::TemplateResponseCalculatorBase<3, false> {

  enum bin_indices { kIndex_Enu = 0, kIndex_q0 = 1, kIndex_q3 = 2 };

  std::map<std::string, systtools::paramId_t> ParamNames;

public:
  MKSinglePiEnuq0q3_ReWeight(std::map<std::string, systtools::paramId_t> params,
                             fhicl::ParameterSet const &InputManifest)
      : ParamNames(params) {

    LoadInputHistograms(InputManifest, ParamNames);
    ValidateInputHistograms();
  }

  virtual bin_it_t GetBin(systtools::paramId_t pId,
                          std::array<double, 3> const &kinematics) {
    if (BinnedResponses.find(pId) == BinnedResponses.end()) {
      throw systtools::invalid_parameter_Id()
          << "[ERROR]: Template reweight bin requested for parameter " << pId
          << ", but " << GetCalculatorName()
          << " does not handle this parameter.";
    }

    TH3D *firstHist = BinnedResponses.at(pId).begin()->second.get();

    Int_t XBin = firstHist->GetXaxis()->FindFixBin(kinematics[kIndex_q0]);
    if (IsFlowBin(firstHist->GetXaxis(), XBin)) {
      throw is_flow_bin() << "[WARN]: q0: " << kinematics[kIndex_q0]
                          << " outside normal bin range: "
                          << firstHist->GetXaxis()->GetBinLowEdge(1) << " -- "
                          << firstHist->GetXaxis()->GetBinUpEdge(
                                 firstHist->GetXaxis()->GetNbins());
    }
    Int_t YBin = firstHist->GetYaxis()->FindFixBin(kinematics[kIndex_q3]);
    if (IsFlowBin(firstHist->GetYaxis(), YBin)) {
      throw is_flow_bin() << "[WARN]: q3: " << kinematics[kIndex_q3]
                          << " outside normal bin range: "
                          << firstHist->GetYaxis()->GetBinLowEdge(1) << " -- "
                          << firstHist->GetYaxis()->GetBinUpEdge(
                                 firstHist->GetYaxis()->GetNbins());
    }
    Int_t ZBin = firstHist->GetZaxis()->FindFixBin(kinematics[kIndex_Enu]);
    if (IsFlowBin(firstHist->GetZaxis(), ZBin)) {
      throw is_flow_bin() << "[WARN]: Enu: " << kinematics[kIndex_Enu]
                          << " outside normal bin range: "
                          << firstHist->GetZaxis()->GetBinLowEdge(1) << " -- "
                          << firstHist->GetZaxis()->GetBinUpEdge(
                                 firstHist->GetZaxis()->GetNbins());
    }

#ifdef MKSinglePiEnuq0q3_ReWeight_DEBUG
    std::cout << "[INFO]: Getting bin Enu: " << kinematics[kIndex_Enu]
              << ", q0: " << kinematics[kIndex_q0]
              << ", q3: " << kinematics[kIndex_q3] << " = { " << XBin << ", "
              << YBin << ", " << ZBin
              << " } = " << firstHist->GetBin(XBin, YBin, ZBin) << std::endl;
#endif
    return firstHist->GetBin(XBin, YBin, ZBin);
  }

  std::string GetCalculatorName() { return "MKSinglePiEnuq0q3_ReWeight"; }
};
} // namespace nusyst

#endif
