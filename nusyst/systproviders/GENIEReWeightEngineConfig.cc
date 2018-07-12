#include "GENIEReWeightEngineConfig.hh"

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

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> ConfigureQEWeightEngine(
    SystMetaData const &QEmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> param_map;

  bool UseDipoleAxial = HasParam(QEmd, "NormCCQE") || HasParam(QEmd, "MaCCQE");
  bool UseZExp =
      HasParam(QEmd, "ZNormCCQE") || HasParam(QEmd, "ZExpAVariationResponse");

  if (UseDipoleAxial || UseZExp) {
    GReWeightEngine->AdoptWghtCalc("xsec_ccqe_axFF",
                                   new genie::rew::GReWeightNuXSecCCQE);
    genie::rew::GReWeightNuXSecCCQE *rwccqe =
        dynamic_cast<genie::rew::GReWeightNuXSecCCQE *>(
            GReWeightEngine->WghtCalc("xsec_ccqe_axFF"));

    if (UseDipoleAxial) {
      size_t nqeidx = GetParamIndex(QEmd, "NormCCQE");
      size_t maqeidx = GetParamIndex(QEmd, "MaCCQE");
      bool maqeHasShapeOpt = SystHasOpt(QEmd, "MaCCQE", "shape");

      bool IsShape = HasParam(QEmd, "NormCCQE") || maqeHasShapeOpt;

      rwccqe->SetMode(IsShape
                          ? genie::rew::GReWeightNuXSecCCQE::kModeNormAndMaShape
                          : genie::rew::GReWeightNuXSecCCQE::kModeMa);

      if (HasParam(QEmd, "NormCCQE")) {
        param_map[nqeidx].insert({genie::rew::kXSecTwkDial_NormCCQE, nqeidx});
      }

      if (HasParam(QEmd, "MaCCQE")) {
        param_map[maqeidx].insert({IsShape
                                       ? genie::rew::kXSecTwkDial_MaCCQEshape
                                       : genie::rew::kXSecTwkDial_MaCCQE,
                                   maqeidx});
      }

    } else {
      rwccqe->SetMode(genie::rew::GReWeightNuXSecCCQE::kModeZExp);

      if (HasParam(QEmd, "AxFFCCQEshape")) {
        size_t zaxidx = GetParamIndex(QEmd, "AxFFCCQEshape");
        GReWeightEngine->AdoptWghtCalc(
            "xsec_ccqe_ZExp", new genie::rew::GReWeightNuXSecCCQEaxial);
        param_map[zaxidx].insert(
            {genie::rew::kXSecTwkDial_AxFFCCQEshape, zaxidx});
      }

      if (HasParam(QEmd, "ZNormCCQE")) {
        size_t znormidx = GetParamIndex(QEmd, "ZNormCCQE");
        param_map[znormidx].insert(
            {genie::rew::kXSecTwkDial_ZNormCCQE, znormidx});
      }

      // Turns of each zexpansion parameter effect a response via a single
      // parameter.
      size_t zrespidx = GetParamIndex(QEmd, "ZExpAVariationResponse");
      if (HasParam(QEmd, "ZExpA1CCQE")) {
        size_t za1idx = GetParamIndex(QEmd, "ZExpA1CCQE");
        param_map[zrespidx].insert(
            {genie::rew::kXSecTwkDial_ZExpA1CCQE, za1idx});
      }
      if (HasParam(QEmd, "ZExpA2CCQE")) {
        size_t za2idx = GetParamIndex(QEmd, "ZExpA2CCQE");
        param_map[zrespidx].insert(
            {genie::rew::kXSecTwkDial_ZExpA2CCQE, za2idx});
      }
      if (HasParam(QEmd, "ZExpA3CCQE")) {
        size_t za3idx = GetParamIndex(QEmd, "ZExpA3CCQE");
        param_map[zrespidx].insert(
            {genie::rew::kXSecTwkDial_ZExpA3CCQE, za3idx});
      }
      if (HasParam(QEmd, "ZExpA4CCQE")) {
        size_t za4idx = GetParamIndex(QEmd, "ZExpA4CCQE");
        param_map[zrespidx].insert(
            {genie::rew::kXSecTwkDial_ZExpA4CCQE, za4idx});
      }
    }
  }

  if (HasParam(QEmd, "VecFFCCQEshape")) {
    GReWeightEngine->AdoptWghtCalc("xsec_ccqe_vecFF",
                                   new genie::rew::GReWeightNuXSecCCQEvec);
    size_t vecffidx = GetParamIndex(QEmd, "VecFFCCQEshape");
    param_map[vecffidx].insert(
        {genie::rew::kXSecTwkDial_VecFFCCQEshape, vecffidx});
  }
  return param_map;
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureNCELWeightEngine(
    SystMetaData const &NCELmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> param_map;

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

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureRESWeightEngine(
    SystMetaData const &RESmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> param_map;

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

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureCOHWeightEngine(
    SystMetaData const &COHmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> param_map;
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

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureDISWeightEngine(
    SystMetaData const &DISmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> param_map;

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

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureFSIWeightEngine(
    SystMetaData const &FSImd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureOtherWeightEngine(
    SystMetaData const &Othermd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

} // namespace nusyst
