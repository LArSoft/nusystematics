#include "GENIEReWeightEngineConfig.hh"

#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEaxial.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"

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

    genie::rew::GReWeightNuXSecCCRes *rwccres =
        dynamic_cast<genie::rew::GReWeightNuXSecCCRes *>(
            GReWeightEngine->WghtCalc("xsec_CCRES"));

    bool MaCCRESHasShapeOpt = SystHasOpt(RESmd, "MaCCRES", "shape");
    bool MvCCRESHasShapeOpt = SystHasOpt(RESmd, "MvCCRES", "shape");
    bool IsShape = MaCCRESHasShapeOpt || MvCCRESHasShapeOpt ||
                   IndexIsHandled(RESmd, NormCCRESidx);

    rwccres->SetMode(
        IsShape ? genie::rew::GReWeightNuXSecCCRes::kModeNormAndMaMvShape
                : genie::rew::GReWeightNuXSecCCRes::kModeMaMv);

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
                                   new genie::rew::GReWeightNuXSecNCRRES);

    size_t NCRESrespidx = GetParamIndex(RESmd, "NCRESVariationResponse");
    size_t MaNCRESidx = GetParamIndex(RESmd, "MaNCRES");
    size_t MvNCRESidx = GetParamIndex(RESmd, "MvNCRES");
    size_t NormNCRESidx = GetParamIndex(RESmd, "NormNCRES");

    genie::rew::GReWeightNuXSecNCRes *rwccres =
        dynamic_cast<genie::rew::GReWeightNuXSecNCRes *>(
            GReWeightEngine->WghtCalc("xsec_NCRES"));

    bool MaNCRESHasShapeOpt = SystHasOpt(RESmd, "MaNCRES", "shape");
    bool MvNCRESHasShapeOpt = SystHasOpt(RESmd, "MvNCRES", "shape");
    bool IsShape = MaNCRESHasShapeOpt || MvNCRESHasShapeOpt ||
                   IndexIsHandled(RESmd, NormNCRESidx);

    rwccres->SetMode(
        IsShape ? genie::rew::GReWeightNuXSecNCRes::kModeNormAndMaMvShape
                : genie::rew::GReWeightNuXSecNCRes::kModeMaMv);

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
    genie::rew::GSyst_t syst = genie::rew::GSyst::FromString(pname);
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
    genie::rew::GSyst_t syst = genie::rew::GSyst::FromString(pname);
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
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureDISWeightEngine(
    SystMetaData const &DISmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
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
