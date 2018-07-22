#include "nusyst/systproviders/GENIEReWeightEngineConfig.hh"

// GENIE
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEaxial.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecNCEL.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightResonanceDecay.h"

using namespace larsyst;

namespace nusyst {

std::vector<GENIEResponseParameter>
ConfigureQEWeightEngine(SystMetaData const &QEmd,
                        fhicl::ParameterSet const &tool_options) {

  std::vector<GENIEResponseParameter> param_map;

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  bool UseFullHERG = tool_options.get<bool>("UseFullHERG", false);

  if (HasAnyParams(QEmd, {"NormCCQE", "MaCCQE"})) {

    if (HasParam(QEmd, "NormCCQE")) {
      size_t nqeidx = GetParamIndex(QEmd, "NormCCQE");
      GENIEResponseParameter nqer;
      nqer.pidx = nqeidx;
      nqer.dependents.push_back({genie::rew::kXSecTwkDial_NormCCQE, nqeidx});

      for (double var : QEmd[nqeidx].paramVariations) {
        std::unique_ptr<genie::rew::GReWeight> grw =
            std::make_unique<genie::rew::GReWeight>();

        grw->AdoptWghtCalc("xsec_ccqe_axFF",
                           new genie::rew::GReWeightNuXSecCCQE);
        grw->Systematics().Init(genie::rew::kXSecTwkDial_NormCCQE, var);

        grw->Reconfigure();
        nqer.Herg.push_back(std::move(grw));
        if (!UseFullHERG) {
          break;
        }
      }

      param_map.push_back(std::move(nqer));
    }

    if (HasParam(QEmd, "MaCCQE")) {
      size_t maqeidx = GetParamIndex(QEmd, "MaCCQE");
      bool IsShape = tool_options.get<bool>("MAQEIsShapeOnly", false);
      GENIEResponseParameter maqer;
      maqer.pidx = maqeidx;
      maqer.dependents.push_back({IsShape ? genie::rew::kXSecTwkDial_MaCCQEshape
                                          : genie::rew::kXSecTwkDial_MaCCQE,
                                  maqeidx});

      for (double var : QEmd[maqeidx].paramVariations) {
        std::unique_ptr<genie::rew::GReWeight> grw =
            std::make_unique<genie::rew::GReWeight>();

        genie::rew::GReWeightNuXSecCCQE *rwccqe =
            new genie::rew::GReWeightNuXSecCCQE();

        rwccqe->SetMode(
            IsShape ? genie::rew::GReWeightNuXSecCCQE::kModeNormAndMaShape
                    : genie::rew::GReWeightNuXSecCCQE::kModeMa);

        grw->AdoptWghtCalc("xsec_ccqe_axFF", rwccqe);

        grw->Systematics().Init(IsShape ? genie::rew::kXSecTwkDial_MaCCQEshape
                                        : genie::rew::kXSecTwkDial_MaCCQE,
                                var);

        grw->Reconfigure();
        maqer.Herg.push_back(std::move(grw));
        if (!UseFullHERG) {
          break;
        }
      }

      param_map.push_back(std::move(maqer));
    }
  }

  if (HasParam(QEmd, "AxFFCCQEshape")) {
    size_t zaxidx = GetParamIndex(QEmd, "AxFFCCQEshape");

    GENIEResponseParameter axffqer;
    axffqer.pidx = zaxidx;
    axffqer.dependents.push_back(
        {genie::rew::kXSecTwkDial_AxFFCCQEshape, zaxidx});

    std::unique_ptr<genie::rew::GReWeight> grw =
        std::make_unique<genie::rew::GReWeight>();

    grw->AdoptWghtCalc("xsec_ccqe_axFF_form",
                       new genie::rew::GReWeightNuXSecCCQEaxial());

    grw->Systematics().Init(genie::rew::kXSecTwkDial_AxFFCCQEshape, 1);

    grw->Reconfigure();

    axffqer.Herg.push_back(std::move(grw));

    param_map.push_back(std::move(axffqer));
  }

  if (HasAnyParams(QEmd, {"ZNormCCQE", "ZExpAVariationResponse", "ZExpA1CCQE",
                          "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE"})) {

    if (HasParam(QEmd, "ZNormCCQE")) {
      size_t znormidx = GetParamIndex(QEmd, "ZNormCCQE");
      GENIEResponseParameter znr;
      znr.pidx = znormidx;
      znr.dependents.push_back({genie::rew::kXSecTwkDial_ZNormCCQE, znormidx});
      for (double var : QEmd[znormidx].paramVariations) {
        std::unique_ptr<genie::rew::GReWeight> grw =
            std::make_unique<genie::rew::GReWeight>();

        grw->AdoptWghtCalc("xsec_ccqe_axFF",
                           new genie::rew::GReWeightNuXSecCCQE);
        grw->Systematics().Init(genie::rew::kXSecTwkDial_ZNormCCQE, var);

        grw->Reconfigure();
        znr.Herg.push_back(std::move(grw));
        if (!UseFullHERG) {
          break;
        }
      }
    }

    // If not ignoring dependence then have to hook up multiple parameter
    // variations to the response parameter.
    if (HasParam(QEmd, "ZExpAVariationResponse") &&
        HasAnyParams(
            QEmd, {"ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE"})) {

      size_t zexpridx = GetParamIndex(QEmd, "ZExpAVariationResponse");
      GENIEResponseParameter Zexpr;
      Zexpr.pidx = zexpridx;

      for (auto const &ZExpANameDial :
           std::vector<std::pair<std::string, genie::rew::GSyst_t>>{
               {"ZExpA1CCQE", genie::rew::kXSecTwkDial_ZExpA1CCQE},
               {"ZExpA2CCQE", genie::rew::kXSecTwkDial_ZExpA2CCQE},
               {"ZExpA3CCQE", genie::rew::kXSecTwkDial_ZExpA3CCQE},
               {"ZExpA4CCQE", genie::rew::kXSecTwkDial_ZExpA4CCQE}}) {
        if (HasParam(QEmd, ZExpANameDial.first)) {
          Zexpr.dependents.push_back(
              {ZExpANameDial.second, GetParamIndex(QEmd, ZExpANameDial.first)});
        }
      }

      for (size_t i = 0; i < QEmd[zexpridx].paramVariations.size(); ++i) {
        std::unique_ptr<genie::rew::GReWeight> grw =
            std::make_unique<genie::rew::GReWeight>();

        grw->AdoptWghtCalc("xsec_ccqe_axFF",
                           new genie::rew::GReWeightNuXSecCCQE);

        for (auto const &dep : Zexpr.dependents) {
          grw->Systematics().Init(dep.gdial, QEmd[dep.pidx].paramVariations[i]);
        }

        grw->Reconfigure();
        Zexpr.Herg.push_back(std::move(grw));
        if (!UseFullHERG) {
          break;
        }
      }
      // We are ignoring the dependence of the zexpansion parameters.
    } else if (HasAnyParams(QEmd, {"ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE",
                                   "ZExpA4CCQE"})) {

      for (auto const &ZExpANameDial :
           std::vector<std::pair<std::string, genie::rew::GSyst_t>>{
               {"ZExpA1CCQE", genie::rew::kXSecTwkDial_ZExpA1CCQE},
               {"ZExpA2CCQE", genie::rew::kXSecTwkDial_ZExpA2CCQE},
               {"ZExpA3CCQE", genie::rew::kXSecTwkDial_ZExpA3CCQE},
               {"ZExpA4CCQE", genie::rew::kXSecTwkDial_ZExpA4CCQE}}) {

        if (HasParam(QEmd, ZExpANameDial.first)) {
          size_t zdialidx = GetParamIndex(QEmd, ZExpANameDial.first);
          GENIEResponseParameter Zexpar;
          Zexpar.pidx = zdialidx;
          Zexpar.dependents.push_back({ZExpANameDial.second, zdialidx});
          for (double var : QEmd[zdialidx].paramVariations) {
            std::unique_ptr<genie::rew::GReWeight> grw =
                std::make_unique<genie::rew::GReWeight>();

            grw->AdoptWghtCalc("xsec_ccqe_axFF",
                               new genie::rew::GReWeightNuXSecCCQE);
            grw->Systematics().Init(ZExpANameDial.second, var);

            grw->Reconfigure();
            Zexpar.Herg.push_back(std::move(grw));
            if (!UseFullHERG) {
              break;
            }
          }
        }
      }
    }
  }

  if (HasParam(QEmd, "VecFFCCQEshape")) {
    size_t vecffidx = GetParamIndex(QEmd, "VecFFCCQEshape");

    GENIEResponseParameter vecffqer;
    vecffqer.pidx = vecffidx;
    vecffqer.dependents.push_back(
        {genie::rew::kXSecTwkDial_VecFFCCQEshape, vecffidx});

    std::unique_ptr<genie::rew::GReWeight> grw =
        std::make_unique<genie::rew::GReWeight>();

    grw->AdoptWghtCalc("xsec_ccqe_vecFF",
                       new genie::rew::GReWeightNuXSecCCQEvec);

    grw->Systematics().Init(genie::rew::kXSecTwkDial_VecFFCCQEshape,
                            QEmd[vecffidx].centralParamValue);
    grw->Reconfigure();

    vecffqer.Herg.push_back(std::move(grw));

    param_map.push_back(std::move(vecffqer));
  }
  return param_map;
} // namespace nusyst

#ifndef GRWTEST

std::vector<GENIEResponseParameter> ConfigureNCELWeightEngine(
    SystMetaData const &NCELmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::vector<GENIEResponseParameter> param_map;

  bool UseNCEL = HasParam(NCELmd, "MaNCEL") || HasParam(NCELmd, "EtaNCEL");

  if (UseNCEL) {
    GReWeightEngine->AdoptWghtCalc("xsec_ncel",
                                   new genie::rew::GReWeightNuXSecNCEL);

    size_t NCELrespidx = GetParamIndex(NCELmd, "NCELVariationResponse");
    size_t MaNCELidx = GetParamIndex(NCELmd, "MaNCEL");
    size_t EtaNCELidx = GetParamIndex(NCELmd, "EtaNCEL");

    if (IndexIsHandled(NCELmd, MaNCELidx)) {
      param_map[NCELrespidx].insert(
          {genie::rew::kXSecTwkDial_MaNCEL, MaNCELidx});
    }

    if (IndexIsHandled(NCELmd, EtaNCELidx)) {
      param_map[NCELrespidx].insert(
          {genie::rew::kXSecTwkDial_EtaNCEL, EtaNCELidx});
    }
  }
  return param_map;
}

std::vector<GENIEResponseParameter> ConfigureRESWeightEngine(
    SystMetaData const &RESmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::vector<GENIEResponseParameter> param_map;

  bool UseCCRES = HasParam(RESmd, "NormCCRES") || HasParam(RESmd, "MaCCRES") ||
                  HasParam(RESmd, "MvCCRES");

  if (UseCCRES) {
    GReWeightEngine->AdoptWghtCalc("xsec_CCRES",
                                   new genie::rew::GReWeightNuXSecCCRES);

    size_t CCRESrespidx = GetParamIndex(RESmd, "CCRESVariationResponse");
    size_t NormCCRESidx = GetParamIndex(RESmd, "NormCCRES");
    size_t MaCCRESidx = GetParamIndex(RESmd, "MaCCRES");
    size_t MvCCRESidx = GetParamIndex(RESmd, "MvCCRES");

    genie::rew::GReWeightNuXSecCCRES *rwccres =
        dynamic_cast<genie::rew::GReWeightNuXSecCCRES *>(
            GReWeightEngine->WghtCalc("xsec_CCRES"));

    bool MaCCRESHasShapeOpt = SystHasOpt(RESmd, "MaCCRES", "shape");
    bool MvCCRESHasShapeOpt = SystHasOpt(RESmd, "MvCCRES", "shape");
    bool IsShape = MaCCRESHasShapeOpt || MvCCRESHasShapeOpt ||
                   IndexIsHandled(RESmd, NormCCRESidx);

    rwccres->SetMode(
        IsShape ? genie::rew::GReWeightNuXSecCCRES::kModeNormAndMaMvShape
                : genie::rew::GReWeightNuXSecCCRES::kModeMaMv);

    if (IndexIsHandled(RESmd, NormCCRESidx)) {
      param_map[NormCCRESidx].insert(
          {genie::rew::kXSecTwkDial_NormCCRES, NormCCRESidx});
    }

    if (IndexIsHandled(RESmd, MaCCRESidx)) {
      param_map[CCRESrespidx].insert(
          {IsShape ? genie::rew::kXSecTwkDial_MaCCRESshape
                   : genie::rew::kXSecTwkDial_MaCCRES,
           MaCCRESidx});
    }

    if (IndexIsHandled(RESmd, MvCCRESidx)) {
      param_map[CCRESrespidx].insert(
          {IsShape ? genie::rew::kXSecTwkDial_MvCCRESshape
                   : genie::rew::kXSecTwkDial_MvCCRES,
           MvCCRESidx});
    }
  }

  bool UseNCRES = HasParam(RESmd, "NormNCRES") || HasParam(RESmd, "MaNCRES") ||
                  HasParam(RESmd, "MvNCRES");

  if (UseNCRES) {
    GReWeightEngine->AdoptWghtCalc("xsec_NCRES",
                                   new genie::rew::GReWeightNuXSecNCRES);

    size_t NCRESrespidx = GetParamIndex(RESmd, "NCRESVariationResponse");
    size_t MaNCRESidx = GetParamIndex(RESmd, "MaNCRES");
    size_t MvNCRESidx = GetParamIndex(RESmd, "MvNCRES");
    size_t NormNCRESidx = GetParamIndex(RESmd, "NormNCRES");

    genie::rew::GReWeightNuXSecNCRES *rwncres =
        dynamic_cast<genie::rew::GReWeightNuXSecNCRES *>(
            GReWeightEngine->WghtCalc("xsec_NCRES"));

    bool MaNCRESHasShapeOpt = SystHasOpt(RESmd, "MaNCRES", "shape");
    bool MvNCRESHasShapeOpt = SystHasOpt(RESmd, "MvNCRES", "shape");
    bool IsShape = MaNCRESHasShapeOpt || MvNCRESHasShapeOpt ||
                   IndexIsHandled(RESmd, NormNCRESidx);

    rwncres->SetMode(
        IsShape ? genie::rew::GReWeightNuXSecNCRES::kModeNormAndMaMvShape
                : genie::rew::GReWeightNuXSecNCRES::kModeMaMv);

    if (IndexIsHandled(RESmd, NormNCRESidx)) {
      param_map[NormNCRESidx].insert(
          {genie::rew::kXSecTwkDial_NormNCRES, NormNCRESidx});
    }

    if (IndexIsHandled(RESmd, MaNCRESidx)) {
      param_map[NCRESrespidx].insert(
          {IsShape ? genie::rew::kXSecTwkDial_MaNCRESshape
                   : genie::rew::kXSecTwkDial_MaNCRES,
           MaNCRESidx});
    }

    if (IndexIsHandled(RESmd, MvNCRESidx)) {
      param_map[NCRESrespidx].insert(
          {IsShape ? genie::rew::kXSecTwkDial_MvNCRESshape
                   : genie::rew::kXSecTwkDial_MvNCRES,
           MvNCRESidx});
    }
  }

  bool AdoptedBkgCalc = false;
  for (std::string const &pname :
       {"RvpCC1pi", "RvpCC2pi", "RvpNC1pi", "RvpNC2pi", "RvnCC1pi", "RvnCC2pi",
        "RvnNC1pi", "RvnNC2pi", "RvbarpCC1pi", "RvbarpCC2pi", "RvbarpNC1pi",
        "RvbarpNC2pi", "RvbarnCC1pi", "RvbarnCC2pi", "RvbarnNC1pi",
        "RvbarnNC2pi"}) {
    if (!HasParam(RESmd, pname)) {
      continue;
    }
    if (!AdoptedBkgCalc) {
      GReWeightEngine->AdoptWghtCalc("xsec_NonResBkg",
                                     new genie::rew::GReWeightNonResonanceBkg);
      AdoptedBkgCalc = true;
    }

    size_t NRRESidx = GetParamIndex(RESmd, pname);
    std::string GENIEFromStringName = "NonRESBG" + pname.substr(1);
    genie::rew::GSyst_t syst =
        genie::rew::GSyst::FromString(GENIEFromStringName);
    if (syst == genie::rew::kNullSystematic) {
      std::cout << "[ERROR]: Expected to be able to parse "
                << std::quoted(pname)
                << " as a GENIE Systematic dial, but failed." << std::endl;
      throw;
    }
    param_map[NRRESidx].insert({syst, NRRESidx});
  }
  bool AdoptedResDecayCalc = false;
  for (std::string const &pname : {"BR1gamma", "BR1eta", "Theta_Delta2Npi"}) {
    if (!HasParam(RESmd, pname)) {
      continue;
    }
    if (!AdoptedResDecayCalc) {
      GReWeightEngine->AdoptWghtCalc("xsec_ResDecay",
                                     new genie::rew::GReWeightResonanceDecay);
      AdoptedResDecayCalc = true;
    }

    size_t ResDecayidx = GetParamIndex(RESmd, pname);
    std::string GENIEFromStringName =
        (pname[0] == 'B') ? std::string("RDec") + pname : pname;
    genie::rew::GSyst_t syst =
        genie::rew::GSyst::FromString(GENIEFromStringName);
    if (syst == genie::rew::kNullSystematic) {
      std::cout << "[ERROR]: Expected to be able to parse "
                << std::quoted(pname)
                << " as a GENIE Systematic dial, but failed." << std::endl;
      throw;
    }
    param_map[ResDecayidx].insert({syst, ResDecayidx});
  }

  return param_map;
}

std::vector<GENIEResponseParameter> ConfigureCOHWeightEngine(
    SystMetaData const &COHmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::vector<GENIEResponseParameter> param_map;
  bool UseCOH = HasParam(COHmd, "MaCOHpi") || HasParam(COHmd, "R0COHpi");

  if (UseCOH) {
    GReWeightEngine->AdoptWghtCalc("xsec_COH",
                                   new genie::rew::GReWeightNuXSecCOH);

    size_t COHrespidx = GetParamIndex(COHmd, "COHVariationResponse");
    size_t MaCOHpiidx = GetParamIndex(COHmd, "MaCOHpi");
    size_t R0COHpiidx = GetParamIndex(COHmd, "R0COHpi");

    if (IndexIsHandled(COHmd, MaCOHpiidx)) {
      param_map[COHrespidx].insert(
          {genie::rew::kXSecTwkDial_MaCOHpi, MaCOHpiidx});
    }

    if (IndexIsHandled(COHmd, R0COHpiidx)) {
      param_map[COHrespidx].insert(
          {genie::rew::kXSecTwkDial_R0COHpi, R0COHpiidx});
    }
  }
  return param_map;
}

std::vector<GENIEResponseParameter> ConfigureDISWeightEngine(
    SystMetaData const &DISmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::vector<GENIEResponseParameter> param_map;

  bool UseDISXSec = HasParam(DISmd, "AhtBY") || HasParam(DISmd, "BhtBY") ||
                    HasParam(DISmd, "CV1uBY") || HasParam(DISmd, "CV2uBY");

  if (UseDISXSec) {
    GReWeightEngine->AdoptWghtCalc("xsec_DIS",
                                   new genie::rew::GReWeightNuXSecDIS);

    size_t DISrespidx = GetParamIndex(DISmd, "DISVariationResponse");
    size_t AhtBYidx = GetParamIndex(DISmd, "AhtBY");
    size_t BhtBYidx = GetParamIndex(DISmd, "BhtBY");
    size_t CV1uBYidx = GetParamIndex(DISmd, "CV1uBY");
    size_t CV2uBYidx = GetParamIndex(DISmd, "CV2uBY");

    bool AhtBYHasShapeOpt = SystHasOpt(DISmd, "AhtBY", "shape");
    bool BhtBYHasShapeOpt = SystHasOpt(DISmd, "BhtBY", "shape");
    bool CV1uBYHasShapeOpt = SystHasOpt(DISmd, "CV1uBY", "shape");
    bool CV2uBYHasShapeOpt = SystHasOpt(DISmd, "CV2uBY", "shape");
    bool IsShape = AhtBYHasShapeOpt || BhtBYHasShapeOpt || CV1uBYHasShapeOpt ||
                   CV2uBYHasShapeOpt;

    genie::rew::GReWeightNuXSecDIS *rwdis =
        dynamic_cast<genie::rew::GReWeightNuXSecDIS *>(
            GReWeightEngine->WghtCalc("xsec_DIS"));
    rwdis->SetMode(IsShape ? genie::rew::GReWeightNuXSecDIS::kModeABCV12uShape
                           : genie::rew::GReWeightNuXSecDIS::kModeABCV12u);

    if (IndexIsHandled(DISmd, AhtBYidx)) {
      param_map[DISrespidx].insert(
          {{IsShape ? genie::rew::kXSecTwkDial_AhtBYshape
                    : genie::rew::kXSecTwkDial_AhtBY},
           AhtBYidx});
    }

    if (IndexIsHandled(DISmd, BhtBYidx)) {
      param_map[DISrespidx].insert(
          {{IsShape ? genie::rew::kXSecTwkDial_BhtBYshape
                    : genie::rew::kXSecTwkDial_BhtBY},
           BhtBYidx});
    }

    if (IndexIsHandled(DISmd, CV1uBYidx)) {
      param_map[DISrespidx].insert(
          {{IsShape ? genie::rew::kXSecTwkDial_CV1uBYshape
                    : genie::rew::kXSecTwkDial_CV1uBY},
           CV1uBYidx});
    }

    if (IndexIsHandled(DISmd, CV2uBYidx)) {
      param_map[DISrespidx].insert(
          {{IsShape ? genie::rew::kXSecTwkDial_CV2uBYshape
                    : genie::rew::kXSecTwkDial_CV2uBY},
           CV2uBYidx});
    }
  }

  bool UseAGKY = HasParam(DISmd, "AGKY_xF1pi") || HasParam(DISmd, "AGKY_pT1pi");

  if (UseAGKY) {
    GReWeightEngine->AdoptWghtCalc("xsec_AGKY", new genie::rew::GReWeightAGKY);

    size_t AGKYrespidx = GetParamIndex(DISmd, "AGKYVariationResponse");
    size_t AGKY_xF1piidx = GetParamIndex(DISmd, "AGKY_xF1pi");
    size_t AGKY_pT1piidx = GetParamIndex(DISmd, "AGKY_pT1pi");

    if (IndexIsHandled(DISmd, AGKY_xF1piidx)) {
      param_map[AGKYrespidx].insert(
          {genie::rew::kHadrAGKYTwkDial_xF1pi, AGKY_xF1piidx});
    }

    if (IndexIsHandled(DISmd, AGKY_pT1piidx)) {
      param_map[AGKYrespidx].insert(
          {genie::rew::kHadrAGKYTwkDial_pT1pi, AGKY_pT1piidx});
    }
  }

  size_t FormZoneidx = GetParamIndex(DISmd, "FormZone");
  if (IndexIsHandled(DISmd, FormZoneidx)) {
    GReWeightEngine->AdoptWghtCalc("xsec_FormZone",
                                   new genie::rew::GReWeightFZone);
    param_map[FormZoneidx].insert(
        {genie::rew::kHadrNuclTwkDial_FormZone, FormZoneidx});
  }

  return param_map;
}

std::vector<GENIEResponseParameter> ConfigureFSIWeightEngine(
    SystMetaData const &FSImd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::vector<GENIEResponseParameter> ConfigureOtherWeightEngine(
    SystMetaData const &Othermd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

#endif

} // namespace nusyst
