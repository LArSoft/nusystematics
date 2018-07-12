#ifndef NUSYST_RESPONSE_CALCULATORS_MINERVARPAQ0Q3_REWEIGHT_HH_SEEN
#define NUSYST_RESPONSE_CALCULATORS_MINERVARPAQ0Q3_REWEIGHT_HH_SEEN

#include "nusyst/responsecalculators/TemplateResponseCalculatorBase.hh"

#include "nusyst/responsecalculators/MINERvAq0q3Weighting_data.hh"
#include "nusyst/utility/enumclass2int.hh"
#include "nusyst/utility/exceptions.hh"

#include "fhiclcpp/make_ParameterSet.h"

#include <ostream>

#define MINERvARPAq0q3_ReWeight_DEBUG

NEW_LARSYST_EXCEPT(invalid_MINERvA_RPA_tweak);

namespace nusyst {
class MINERvARPAq0q3_ReWeight
    : private nusyst::TemplateResponseCalculatorBase<2, false> {

  enum bin_indices { kIndex_q0 = 0, kIndex_q3 = 1 };

  larsyst::paramId_t ResponseParamId;

  constexpr static double const q0_offsetValenciaGENIE_GeV = 1E-2;

  std::array<double, 2> const Q2Lims{{0, 9}};
  std::array<double, 2> const WeightLims{{1E-3, 2}};

public:
  enum class RPATweak_t { kCV = 0, kPlus1 = 1, kMinus1 = -1 };

  MINERvARPAq0q3_ReWeight(std::map<std::string, larsyst::paramId_t> params,
                          std::string const &InputManifest_fhicl) {

    fhicl::ParameterSet ps;
#ifndef NO_ART
    std::unique_ptr<cet::filepath_maker> fm =
        std::make_unique<cet::filepath_lookup>(ev);
    fhicl::make_ParameterSet(InputManifest_fhicl, *fm, ps);
#else
    ps = fhicl::make_ParameterSet(InputManifest_fhicl);
#endif

    if (!params.size()) {
      throw incorrectly_generated()
          << "[ERROR]: MINERvARPAq0q3_ReWeight expected to be passed the "
             "paramId for its response.";
    }

    ResponseParamId = params.begin()->second;

    LoadInputHistograms(ps, params);
    ValidateInputHistograms();
  }

  virtual bin_it_t GetBin(larsyst::paramId_t pId,
                          std::array<double, 2> const &kinematics) {
    if (pId != ResponseParamId) {
      throw invalid_parameter_Id()
          << "[ERROR]: Template reweight bin requested for parameter " << pId
          << ", but " << GetCalculatorName()
          << " does not handle this parameter, it handles: " << ResponseParamId
          << ".";
    }

    std::array<double, 2> kinematics_var = kinematics;

    // Hold events outside of the Valencia calculation phase space at the
    // closest valid bin.
    if (kinematics_var[kIndex_q0] < 0.018) {
      kinematics_var[kIndex_q0] = 0.018 + q0_offsetValenciaGENIE_GeV;
    }

    TH2D *firstHist = BinnedResponses.at(ResponseParamId).begin()->second.get();

    Int_t XBin =
        firstHist->GetXaxis()->FindFixBin(kinematics[kIndex_q3]);
    // Hold events outside of the Valencia calculation phase space at the
    // closest valid bin.
    if (IsFlowBin(firstHist->GetXaxis(), XBin)) {
      XBin = (XBin == 0) ? XBin + 1 : XBin - 1;
    }
#ifdef MINERvARPAq0q3_ReWeight_DEBUG
    std::cout << "\t\tXBin: " << XBin << " from "
              << kinematics[kIndex_q3] << std::endl;
#endif

    Int_t YBin = firstHist->GetYaxis()->FindFixBin(
        kinematics_var[kIndex_q0] - q0_offsetValenciaGENIE_GeV);
    // Hold events outside of the Valencia calculation phase space at the
    // closest valid bin.
    if (IsFlowBin(firstHist->GetYaxis(), YBin)) {
      YBin = (YBin == 0) ? YBin + 1 : YBin - 1;
    }

#ifdef MINERvARPAq0q3_ReWeight_DEBUG
    std::cout << "\t\tYBin: " << YBin << " from "
              << (kinematics_var[kIndex_q0] - q0_offsetValenciaGENIE_GeV)
              << std::endl;
#endif

    return firstHist->GetBin(XBin, YBin);
  }

  double GetWeightQ2(const double Q2_GeV2, RPATweak_t tweak = RPATweak_t::kCV) {

    if (Q2Lims[0] < 0.0) {
      return 1.0;
    }
    if (Q2Lims[1] > 9.0) {
      return 1.0;
    }

    double powerQ2 = 1.0;
    double weight = 0.0;
    for (int ii = 0; ii < 10; ii++) {

      switch (tweak) {
      case RPATweak_t::kCV: {
        weight += nusyst::RPAPolyQ2_CV[ii] * powerQ2;
      }
      case RPATweak_t::kPlus1: {
        weight += nusyst::RPAPolyQ2_Plus1[ii] * powerQ2;
      }
      case RPATweak_t::kMinus1: {
        weight += nusyst::RPAPolyQ2_Minus1[ii] * powerQ2;
      }
      default: { throw invalid_MINERvA_RPA_tweak(); }
      }

      powerQ2 *= Q2_GeV2;
    }
    return weight;
  }

  double GetWeight(double q0_GeV, double q3_GeV,
                   RPATweak_t tweak = RPATweak_t::kCV) {

    double weight = 1;
    double Q2_GeV2 = (q3_GeV * q3_GeV) - (q0_GeV * q0_GeV);

#ifdef MINERvARPAq0q3_ReWeight_DEBUG
    std::cout << "[MINERvARPAq0q3_ReWeight]: Get weight for q0: " << q0_GeV
              << ", "
              << "q3: " << q3_GeV << ", Q2: " << Q2_GeV2 << std::endl;
#endif
    if (Q2_GeV2 < Q2Lims[1]) {
      if (Q2_GeV2 > 3.0) {
        weight = GetWeightQ2(Q2_GeV2, tweak);
      } else {
        int bin2d =
            GetBin(ResponseParamId, std::array<double, 2>{{q0_GeV, q3_GeV}});
#ifdef MINERvARPAq0q3_ReWeight_DEBUG
        std::cout << "\t\tGot bin: " << bin2d << std::endl;
#endif
        weight = GetVariation(ResponseParamId, e2i(tweak), bin2d);

        // now trap bogus entries.  Not sure why they happen, but set to 1.0 not
        // 0.0
        if (weight <= 0.001) {
          weight = 1.0;
        }

        // events in genie but not in valencia should get a weight
        // related to a similar q0 from the bulk distribution.
        if ((q0_GeV < 0.15) && (weight > 0.9)) {
          bin2d = GetBin(ResponseParamId,
                         std::array<double, 2>{{q0_GeV, q3_GeV + 0.15}});
#ifdef MINERvARPAq0q3_ReWeight_DEBUG
          std::cout << "\t\t[INFO]: Moved to bulk bin: " << bin2d << std::endl;
#endif
          weight = GetVariation(ResponseParamId, e2i(tweak), bin2d);
        }
      }
    }

#ifdef MINERvARPAq0q3_ReWeight_DEBUG
    std::cout << "\t\t[INFO]: Final weight: " << weight << std::endl;
#endif

    if ((weight < WeightLims[0]) || (weight > WeightLims[1])) {
      weight = 1.0;
    }

    return weight;
  }

  std::string GetCalculatorName() { return "MINERvARPAq0q3_ReWeight"; }
};
} // namespace nusyst

#endif
