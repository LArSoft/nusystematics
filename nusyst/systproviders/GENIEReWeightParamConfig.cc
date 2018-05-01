#include "GENIEReWeightParamConfig.hh"

#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "larsyst/utility/string_parsers.hh"

#include <iomanip>
#include <iostream>
#include <memory>

// Ugly macros because I'm lazy
#define PARAM_NAME_HELPER(PARAM_NAME, b) PARAM_NAME##b

#define PARAM_IS_USED_BY_CFG(PARAM_NAME)                                       \
  (cfg().PARAM_NAME_HELPER(PARAM_NAME, NominalValue)() != 0xdeadb33f) ||       \
      cfg().PARAM_NAME_HELPER(PARAM_NAME, TweakDefinition)().size();

#define ADD_PARAM_TO_SYST(PARAM_NAME, SYSTDATA)                                \
  larsyst::SystParamHeader param = BuildHeaderFromNomAndTweakDefintion(        \
      cfg().PARAM_NAME_HELPER(PARAM_NAME, NominalValue)(),                     \
      cfg().PARAM_NAME_HELPER(PARAM_NAME, TweakDefinition)()());               \
  param.prettyName = #PARAM_NAME;                                              \
  param.unitsAreNatural = false;                                               \
  param.systParamId = firstParamId++;                                          \
  SYSTDATA.headers.push_back(std::move(param))

#define CHECK_USED_ADD_PARAM_TO_SYST(PARAM_NAME, SYSTDATA)                     \
  bool PARAM_NAME_HELPER(PARAM_NAME, IsUsed) =                                 \
      PARAM_IS_USED_BY_CFG(PARAM_NAME);                                        \
  if (PARAM_NAME_HELPER(PARAM_NAME, IsUsed)) {                                 \
    larsyst::SystParamHeader param = BuildHeaderFromNomAndTweakDefintion(      \
        cfg().PARAM_NAME_HELPER(PARAM_NAME, NominalValue)(),                   \
        cfg().PARAM_NAME_HELPER(PARAM_NAME, TweakDefinition)()());             \
    param.prettyName = #PARAM_NAME;                                            \
    param.unitsAreNatural = false;                                             \
    param.systParamId = firstParamId++;                                        \
    SYSTDATA.headers.push_back(std::move(param))                               \
  }

#define ADD_PARAM_TO_SYST_RESPONSELESS(PARAM_NAME, SYSTDATA,                   \
                                       RESPONSE_PARAM_ID)                      \
  larsyst::SystParamHeader param = BuildHeaderFromNomAndTweakDefintion(        \
      cfg().PARAM_NAME_HELPER(PARAM_NAME, NominalValue)(),                     \
      cfg().PARAM_NAME_HELPER(PARAM_NAME, TweakDefinition)()());               \
  param.prettyName = #PARAM_NAME;                                              \
  param.unitsAreNatural = false;                                               \
  param.systParamId = firstParamId++;                                          \
  param.isResponselessParam = true;                                            \
  param.responseParamId = RESPONSE_PARAM_ID;                                   \
  if (param.isSplineable) {                                                    \
    std::cout                                                                  \
        << "[ERROR]: Attempted to build spline from parameter " #PARAM_NAME    \
           ", which enters into an intrinsically multi-parameter response "    \
           "calculation. As multi-dimensional splines are not yet supported, " \
           "this parameter must be used as a multi-sim or with a set of "      \
           "hand-picked offsets that are the same for each parameter that "    \
           "goes into the response calculation."                               \
        << std::endl;                                                          \
    throw;                                                                     \
  }                                                                            \
  SYSTDATA.headers.push_back(std::move(param))

#define GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(PARAMETER_NAME_LIST,        \
                                                   SYSTDATA, NVARPARAM)        \
  NVARPARAM = 0;                                                               \
  for (auto &name : PARAMETER_NAME_LIST) {                                     \
    if (!HasParam(SYSTDATA, name)) {                                           \
      continue;                                                                \
    }                                                                          \
    if (!NVARPARAM) {                                                          \
      NVARPARAM = GetParam(SYSTDATA, name).paramVariations.size();             \
      continue;                                                                \
    }                                                                          \
    if (NVARPARAM != GetParam(SYSTDATA, name).paramVariations.size()) {        \
      std::cout                                                                \
          << "[ERROR]: NCEL configuration error. Each form factor "            \
             "parameter specified must have the same number of variations, "   \
          << name << " specifies "                                             \
          << GetParam(SYSTDATA, name).paramVariations.size()                   \
          << ", but a previous parameter specified " << NVARPARAM << "."       \
          << std::endl;                                                        \
      throw;                                                                   \
    }                                                                          \
  }

using namespace larsyst;

namespace nusyst {

void MakeThrowsIfNeeded(SystParamHeader &sph,
                        std::unique_ptr<CLHEP::RandGaussQ> &RNJesus,
                        uint64_t NThrows) {
  if (sph.isRandomlyThrown) {
    double cv =
        (sph.centralParamValue == 0xdeadb33f) ? 0 : sph.centralParamValue;
    for (uint64_t t = 0; t < NThrows; ++t) {
      double thr = RNJesus->fire(0, 1);
      double shift = fabs(thr) * ((thr < 0) ? sph.oneSigmaShifts[0]
                                            : sph.oneSigmaShifts[1]);
      sph.paramVariations.push_back(cv + shift);
    }
  }
}

SystParamHeader
BuildHeaderFromNomAndTweakDefintion(double nominal = 0xdeadb33f,
                                    std::string tweakDefinition = "") {
  SystParamHeader sph;
  sph.centralParamValue = (nominal == 0xdeadb33f) ? 0 : nominal;
  trim(tweakDefinition);

  if (tweakDefinition.size()) {
    char fchar = tweakDefinition.front();
    tweakDefinition = tweakDefinition.substr(1, tweakDefinition.length() - 2);
    trim(tweakDefinition);
    if (fchar == '(') { // Spline knots
      sph.paramVariations = BuildDoubleList(tweakDefinition);
      sph.isSplineable = true;
    } else if (fchar == '[') { // Discrete tweaks
      sph.paramVariations = ParseToVect<double>(tweakDefinition, ",");
    } else if (fchar == '{') { // OneSigmaShifts
      std::vector<double> sigShifts = ParseToVect<double>(tweakDefinition, ",");
      if (sigShifts.size() == 1) {
        sph.oneSigmaShifts[0] = -sigShifts.front();
        sph.oneSigmaShifts[1] = sigShifts.front();
      } else if (sigShifts.size() == 2) {
        sph.oneSigmaShifts[0] = sigShifts.front();
        sph.oneSigmaShifts[1] = sigShifts.back();
      } else {
        std::cout << "[ERROR]: When parsing sigma shifts found "
                  << std::quoted(tweakDefinition)
                  << ", but expected {sigma_both_natural_units}, or "
                     "{sigma_low_natural_units, sigma_up_natural_units}."
                  << std::endl;
      }
      sph.isRandomlyThrown = true;
    } else {
      std::cout << "[ERROR]: Found tweak definition "
                << std::quoted(tweakDefinition)
                << ", but expected to find either, \"{sigma_low_natural_units, "
                   "sigma_up_natural_units}\" or \"[spline knot 1, spline knot "
                   "2, spline knot 3,...]\""
                << std::endl;
      throw;
    }
  } else { // Just use the central value every time
    sph.isCorrection = true;
  }
  return sph;
}

SystMetaData
ConfigureQEParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &cfg,
                            paramId_t firstParamId) {
  SystMetaData QEmd;

  // Axial FFs
  bool DipoleNormCCQEIsUsed = PARAM_IS_USED_BY_CFG(NormCCQE);
  bool DipoleIsShapeOnly = cfg().MAQEIsShapeOnly() || DipoleNormCCQEIsUsed;
  bool DipoleMaCCQEIsUsed = PARAM_IS_USED_BY_CFG(MaCCQE);

  bool IsDipoleReWeight =
      DipoleIsShapeOnly || DipoleNormCCQEIsUsed || DipoleMaIsUsed;

  bool ZNormIsUsed = PARAM_IS_USED_BY_CFG(ZNormCCQE);
  bool ZExpA1IsUsed = PARAM_IS_USED_BY_CFG(ZExpA1CCQE);
  bool ZExpA2IsUsed = PARAM_IS_USED_BY_CFG(ZExpA2CCQE);
  bool ZExpA3IsUsed = PARAM_IS_USED_BY_CFG(ZExpA3CCQE);
  bool ZExpA4IsUsed = PARAM_IS_USED_BY_CFG(ZExpA4CCQE);

  bool IsZExpReWeight = ZNormIsUsed || ZExpA1IsUsed || ZExpA2IsUsed ||
                        ZExpA3IsUsed || ZExpA4IsUsed;

  if (IsDipoleReWeight && IsZExpReWeight) {
    std::cout << "[ERROR]: Both dipole and Z-expansion axial form factor dials "
                 "are specified. This is an incompatible configuration."
              << std::endl;
    throw;
  }

  if (IsDipoleReWeight) {
    if (DipoleNormCCQEIsUsed) {
      ADD_PARAM_TO_SYST(NormCCQE, QEmd);
    }
    if (DipoleMaCCQEIsUsed) {
      ADD_PARAM_TO_SYST(MaCCQE, QEmd);
      if (DipoleIsShapeOnly) {
        GetParam(QEmd, "MaCCQE").opts.push_back("shape");
      }
    }
  } else if (IsZExpReWeight) {
    // ZNorm enters in linearly to the weight calculation so doesn't need to
    // output its response via the a meta-parameter.
    if (ZNormIsUsed) {
      ADD_PARAM_TO_SYST(ZNormCCQE, QEmd);
    }
    if (ZExpA1IsUsed || ZExpA2IsUsed || ZExpA3IsUsed || ZExpA4IsUsed) {
      SystParamHeader ZExp;
      ZExp.prettyName = "ZExpAVariationResponse";
      ZExp.systParamId = firstParamId++;

      if (ZExpA1IsUsed) {
        ADD_PARAM_TO_SYST_RESPONSELESS(ZExpA1CCQE, QEmd, ZExp.systParamId);
      }
      if (ZExpA2IsUsed) {
        ADD_PARAM_TO_SYST_RESPONSELESS(ZExpA2CCQE, QEmd, ZExp.systParamId);
      }
      if (ZExpA3IsUsed) {
        ADD_PARAM_TO_SYST_RESPONSELESS(ZExpA3CCQE, QEmd, ZExp.systParamId);
      }
      if (ZExpA4IsUsed) {
        ADD_PARAM_TO_SYST_RESPONSELESS(ZExpA4CCQE, QEmd, ZExp.systParamId);
      }
      QEmd.headers.push_back(std::move(ZExp));
    }
  }

  bool DipoleVecFFShapeIsUsed = !cfg().VecFFCCQEIsBBA();
  if (DipoleVecFFShapeIsUsed) {
    SystParamHeader vecFFQE;
    vecFFQE.prettyName = "VecFFCCQEshape";
    vecFFQE.systParamId = firstParamId++;
    vecFFQE.isCorrection = true;
    vecFFQE.centralParamValue = 1;
    QEmd.headers.push_back(std::move(vecFFQE));
  }

  bool AxFFCCQEDipoleToZExpIsUsed =
      (cfg().AxFFCCQEDipoleToZExp() == 1) || // If explcitly enabled
      (IsZExpReWeight && (cfg().AxFFCCQEDipoleToZExp() ==
                          -1)); // or if defaulted and zexpansion is enabled.
  if (AxFFCCQEDipoleToZExpIsUsed) {
    SystParamHeader axFFQE;
    axFFQE.prettyName = "AxFFCCQEshape";
    axFFQE.systParamId = firstParamId++;
    axFFQE.isCorrection = true;
    axFFQE.centralParamValue = 1;
    QEmd.headers.push_back(std::move(axFFQE));
  }

  std::unique_ptr<CLHEP::MTwistEngine> RNgine =
      std::make_unique<CLHEP::MTwistEngine>(0);
  std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
      std::make_unique<CLHEP::RandGaussQ>(*RNgine);

  for (auto &hdr : QEmd.headers) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }
  if (IsZExpReWeight) {
    uint64_t NVariations = 0;

    auto const &zexpParamNames = {"ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE",
                                  "ZExpA4CCQE"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(zexpParamNames, QEmd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(QEmd, "ZExpAVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return QEmd;
}

SystMetaData
ConfigureNCELParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &cfg,
                              paramId_t firstParamId) {
  SystMetaData NCELmd;

  bool MaNCELIsUsed = PARAM_IS_USED_BY_CFG(MaNCEL);
  bool EtaNCELIsUsed = PARAM_IS_USED_BY_CFG(EtaNCEL);

  if (MaNCELIsUsed || EtaNCELIsUsed) {
    SystParamHeader NCELResp;
    NCELResp.prettyName = "NCELVariationResponse";
    NCELResp.systParamId = firstParamId++;

    if (MaNCELIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(MaNCEL, NCELmd, NCELResp.systParamId);
    }
    if (MaNCELEtaNCELIsUsedIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(EtaNCEL, NCELmd, NCELResp.systParamId);
    }
    NCELmd.headers.push_back(std::move(NCELResp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : NCELmd.headers) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &NCELParamNames = {"MaNCEL", "EtaNCEL"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(NCELParamNames, NCELmd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(NCELmd, "NCELVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return NCELmd;
}

SystMetaData ConfigureRESParameterHeaders(
    fhicl::Table<GENIEReWeightParamConfig> const &cfg) {
  SystMetaData RESmd;

  //************* CCRES
  bool NormCCRESIsUsed = PARAM_IS_USED_BY_CFG(NormCCRES);
  bool CCRESIsShapeOnly = cfg().CCRESIsShapeOnly() || NormCCRESIsUsed;
  bool CCRESMaIsUsed = PARAM_IS_USED_BY_CFG(CCRESMa);
  bool CCRESMvIsUsed = PARAM_IS_USED_BY_CFG(CCRESMv);

  if (NormCCRESIsUsed) {
    ADD_PARAM_TO_SYST(NormCCRES, RESmd);
  }

  if (CCRESMaIsUsed || CCRESMvIsUsed) {
    SystParamHeader CCRESResp;
    CCRESResp.prettyName = "CCRESVariationResponse";
    CCRESResp.systParamId = firstParamId++;

    if (CCRESMaIsUsed) {
      ADD_PARAM_TO_SYST(MaCCRES, RESmd);
      if (CCRESIsShapeOnly) {
        GetParam(RESmd, "MaCCRES").opts.push_back("shape");
      }
    }
    if (CCRESMvIsUsed) {
      ADD_PARAM_TO_SYST(MvCCRES, RESmd);
      if (CCRESIsShapeOnly) {
        GetParam(RESmd, "MvCCRES").opts.push_back("shape");
      }
    }
    RESmd.headers.push_back(std::move(CCRESResp));
  }

  //************* NCRES
  bool NormNCRESIsUsed = PARAM_IS_USED_BY_CFG(NormNCRES);
  bool NCRESIsShapeOnly = cfg().NCRESIsShapeOnly() || NormNCRESIsUsed;
  bool NCRESMaIsUsed = PARAM_IS_USED_BY_CFG(NCRESMa);
  bool NCRESMvIsUsed = PARAM_IS_USED_BY_CFG(NCRESMv);

  if (NormNCRESIsUsed) {
    ADD_PARAM_TO_SYST(NormNCRES, RESmd);
  }

  if (NCRESMaIsUsed || NCRESMvIsUsed) {
    SystParamHeader NCRESResp;
    NCRESResp.prettyName = "NCRESVariationResponse";
    NCRESResp.systParamId = firstParamId++;

    if (NCRESMaIsUsed) {
      ADD_PARAM_TO_SYST(MaNCRES, RESmd);
      if (NCRESIsShapeOnly) {
        GetParam(RESmd, "MaNCRES").opts.push_back("shape");
      }
    }
    if (NCRESMvIsUsed) {
      ADD_PARAM_TO_SYST(MvNCRES, RESmd);
      if (NCRESIsShapeOnly) {
        GetParam(RESmd, "MvNCRES").opts.push_back("shape");
      }
    }
    RESmd.headers.push_back(std::move(NCRESResp));
  }

  // These are all independent and based upon the channel that was generated
  CHECK_USED_ADD_PARAM_TO_SYST(RvpCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvpCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvpNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvpNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvnCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvnCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvnNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvnNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarpCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarpCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarpNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarpNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarnCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarnCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarnNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RvbarnNC2pi, RESmd);

  CHECK_USED_ADD_PARAM_TO_SYST(BR1gamma, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(BR1eta, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(Theta_Delta2Npi, RESmd);

  std::unique_ptr<CLHEP::MTwistEngine> RNgine =
      std::make_unique<CLHEP::MTwistEngine>(0);
  std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
      std::make_unique<CLHEP::RandGaussQ>(*RNgine);

  for (auto &hdr : RESmd.headers) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }

  uint64_t NCCRESVariations = 0;

  auto const &CCRESParamNames = {"MaCCRES", "MvCCRES"};
  GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(CCRESParamNames, CCRESmd,
                                             NCCRESVariations);
  std::vector<double> dummyParamVars;
  for (size_t i = 0; i < NCCRESVariations; ++i) {
    dummyParamVars.push_back(i);
  }

  GetParam(RESmd, "CCRESVariationResponse").paramVariations =
      std::move(dummyParamVars);

  uint64_t NNCRESVariations = 0;

  auto const &NCRESParamNames = {"MaNCRES", "MvNCRES"};
  GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(NCRESParamNames, NCRESmd,
                                             NNCRESVariations);
  std::vector<double> dummyParamVars;
  for (size_t i = 0; i < NNCRESVariations; ++i) {
    dummyParamVars.push_back(i);
  }

  GetParam(RESmd, "NCRESVariationResponse").paramVariations =
      std::move(dummyParamVars);

  return RESmd;
}

SystMetaData ConfigureCOHParameterHeaders(
    fhicl::Table<GENIEReWeightParamConfig> const &cfg) {
  SystMetaData COHmd;

  bool MaCOHIsUsed = PARAM_IS_USED_BY_CFG(MaCOH);
  bool R0COHIsUsed = PARAM_IS_USED_BY_CFG(R0COH);

  if (MaCOHIsUsed || R0COHIsUsed) {
    SystParamHeader COHResp;
    COHResp.prettyName = "COHVariationResponse";
    COHResp.systParamId = firstParamId++;

    if (MaCOHIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(MaCOH, COHmd, COHResp.systParamId);
    }
    if (R0COHIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(R0COH, COHmd, COHResp.systParamId);
    }
    COHmd.headers.push_back(std::move(COHResp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : COHmd.headers) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &COHParamNames = {"MaCOH", "R0COH"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(COHParamNames, COHmd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(COHmd, "COHVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return COHmd;
}

SystMetaData ConfigureDISParameterHeaders(
    fhicl::Table<GENIEReWeightParamConfig> const &cfg) {
  SystMetaData DISmd;

  bool DISIsShapeOnly = cfg().DISIsShapeOnly();

  bool AhtBYIsUsed = PARAM_IS_USED_BY_CFG(AhtBY);
  bool BhtBYIsUsed = PARAM_IS_USED_BY_CFG(BhtBY);
  bool CV1uBYIsUsed = PARAM_IS_USED_BY_CFG(CV1uBY);
  bool CV2uBYIsUsed = PARAM_IS_USED_BY_CFG(CV2uBY);

  if (AhtBYIsUsed || BhtBYIsUsed || CV1uBYIsUsed || CV2uBYIsUsed) {
    SystParamHeader DISResp;
    DISResp.prettyName = "DISVariationResponse";
    DISResp.systParamId = firstParamId++;

    if (AhtBYIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(AhtBY, DISmd, DISResp.systParamId);
      if (DISIsShapeOnly) {
        GetParam(DISmd, "AhtBY").opts.push_back("shape");
      }
    }
    if (BhtBYIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(BhtBY, DISmd, DISResp.systParamId);
      if (DISIsShapeOnly) {
        GetParam(DISmd, "BhtBY").opts.push_back("shape");
      }
    }
    if (CV1uBYIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(CV1uBY, DISmd, DISResp.systParamId);
      if (DISIsShapeOnly) {
        GetParam(DISmd, "CV1uBY").opts.push_back("shape");
      }
    }
    if (CV2uBYIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(CV2uBY, DISmd, DISResp.systParamId);
      if (DISIsShapeOnly) {
        GetParam(DISmd, "CV2uBY").opts.push_back("shape");
      }
    }
    DISmd.headers.push_back(std::move(DISResp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : DISmd.headers) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &DISParamNames = {"AhtBY", "BhtBY", "CV1uBY", "CV2uBY"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(DISParamNames, DISmd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(DISmd, "DISVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return DISmd;
}

SystMetaData ConfigureFSIParameterHeaders(
    fhicl::Table<GENIEReWeightParamConfig> const &cfg) {
  SystMetaData FSImd;

  bool MFP_piIsUsed = PARAM_IS_USED_BY_CFG(MFP_pi);
  bool FrCEx_piIsUsed = PARAM_IS_USED_BY_CFG(FrCEx_pi);
  bool FrElas_piIsUsed = PARAM_IS_USED_BY_CFG(FrElas_pi);
  bool FrInel_piIsUsed = PARAM_IS_USED_BY_CFG(FrInel_pi);
  bool FrAbs_piIsUsed = PARAM_IS_USED_BY_CFG(FrAbs_pi);
  bool FrPiProd_piIsUsed = PARAM_IS_USED_BY_CFG(FrPiProd_pi);

  if (MFP_pi || FrCEx_pi || FrElas_pi || FrInel_pi || FrAbs_pi || FrPiProd_pi) {
    SystParamHeader FSI_pi_Resp;
    FSI_pi_Resp.prettyName = "FSI_pi_VariationResponse";
    FSI_pi_Resp.systParamId = firstParamId++;

    if (MFP_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(MFP_pi, FSImd, FSI_pi_Resp.systParamId);
    }
    if (FrCEx_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrCEx_pi, FSImd, FSI_pi_Resp.systParamId);
    }
    if (FrElas_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrElas_pi, FSImd, FSI_pi_Resp.systParamId);
    }
    if (FrInel_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrInel_pi, FSImd, FSI_pi_Resp.systParamId);
    }
    if (FrAbs_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrAbs_pi, FSImd, FSI_pi_Resp.systParamId);
    }
    if (FrPiProd_piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrPiProd_pi, FSImd,
                                     FSI_pi_Resp.systParamId);
    }
    FSImd.headers.push_back(std::move(FSI_pi_Resp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : FSImd.headers) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &FSIParamNames = {"MFP_pi",    "FrCEx_pi", "FrElas_pi",
                                 "FrInel_pi", "FrAbs_pi", "FrPiProd_pi"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(FSIParamNames, FSImd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(FSImd, "FSI_pi_VariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  bool MFP_NIsUsed = PARAM_IS_USED_BY_CFG(MFP_N);
  bool FrCEx_NIsUsed = PARAM_IS_USED_BY_CFG(FrCEx_N);
  bool FrElas_NIsUsed = PARAM_IS_USED_BY_CFG(FrElas_N);
  bool FrInel_NIsUsed = PARAM_IS_USED_BY_CFG(FrInel_N);
  bool FrAbs_NIsUsed = PARAM_IS_USED_BY_CFG(FrAbs_N);
  bool FrPiProd_NIsUsed = PARAM_IS_USED_BY_CFG(FrPiProd_N);

  if (MFP_N || FrCEx_N || FrElas_N || FrInel_N || FrAbs_N || FrPiProd_N) {
    SystParamHeader FSI_N_Resp;
    FSI_N_Resp.prettyName = "FSI_N_VariationResponse";
    FSI_N_Resp.systParamId = firstParamId++;

    if (MFP_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(MFP_N, FSImd, FSI_N_Resp.systParamId);
    }
    if (FrCEx_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrCEx_N, FSImd, FSI_N_Resp.systParamId);
    }
    if (FrElas_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrElas_N, FSImd, FSI_N_Resp.systParamId);
    }
    if (FrInel_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrInel_N, FSImd, FSI_N_Resp.systParamId);
    }
    if (FrAbs_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrAbs_N, FSImd, FSI_N_Resp.systParamId);
    }
    if (FrPiProd_NIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(FrPiProd_N, FSImd, FSI_N_Resp.systParamId);
    }
    FSImd.headers.push_back(std::move(FSI_N_Resp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : FSImd.headers) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &FSIParamNames = {"MFP_N",    "FrCEx_N", "FrElas_N",
                                 "FrInel_N", "FrAbs_N", "FrPiProd_N"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(FSIParamNames, FSImd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(FSImd, "FSI_N_VariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return FSImd;
}

SystMetaData ConfigureOtherParameterHeaders(
    fhicl::Table<GENIEReWeightParamConfig> const &cfg) {
  SystMetaData Othermd;

  ADD_PARAM_TO_SYST(CCQEPauliSupViaKF, Othermd);
  ADD_PARAM_TO_SYST(CCQEMomDistroFGtoSF, Othermd);

  size_t CCQEMomDistroFGtoSFIndex =
      GetParamIndex(Othermd, "CCQEMomDistroFGtoSF");
  if (IndexIsHandled(CCQEMomDistroFGtoSFIndex)) {
    Othermd[CCQEMomDistroFGtoSFIndex].paramValidityRange[0] = 0;
    Othermd[CCQEMomDistroFGtoSFIndex].paramValidityRange[1] = 1;
  }

  for (auto &hdr : Othermd.headers) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }

  return Othermd;
}

} // namespace nusyst
