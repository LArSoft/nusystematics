#include "nusystematics/systproviders/GENIEReWeightParamConfig.hh"

#include "systematicstools/interface/ISystProvider_tool.hh"

#include "systematicstools/utility/string_parsers.hh"
#include "systematicstools/utility/ResponselessParamUtility.hh"

#include <iomanip>
#include <iostream>
#include <memory>

using namespace systtools;

namespace nusyst {

SystMetaData ConfigureQEParameterHeaders(fhicl::ParameterSet const &cfg,
                                         paramId_t firstParamId,
                                         fhicl::ParameterSet &tool_options) {
  SystMetaData QEmd;

  bool MaCCQEIsShapeOnly = cfg.get<bool>("MaCCQEIsShapeOnly", false);

  tool_options.put("MaCCQEIsShapeOnly", MaCCQEIsShapeOnly);

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  // Axial FFs
  bool DipoleNormCCQEIsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "NormCCQE");
  bool DipoleIsShapeOnly = MaCCQEIsShapeOnly || DipoleNormCCQEIsUsed;
  bool DipoleMaCCQEIsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "MaCCQE");

  bool IsDipoleReWeight =
      DipoleIsShapeOnly || DipoleNormCCQEIsUsed || DipoleMaCCQEIsUsed;

  bool ZNormIsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "ZNormCCQE");
  bool ZExpA1IsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "ZExpA1CCQE");
  bool ZExpA2IsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "ZExpA2CCQE");
  bool ZExpA3IsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "ZExpA3CCQE");
  bool ZExpA4IsUsed =
      FHiCLSimpleToolConfigurationParameterExists(cfg, "ZExpA4CCQE");

  bool IsZExpReWeight = ZNormIsUsed || ZExpA1IsUsed || ZExpA2IsUsed ||
                        ZExpA3IsUsed || ZExpA4IsUsed;

  if (IsDipoleReWeight && IsZExpReWeight) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIReWeight_tool, both dipole and "
           "Z-expansion axial form factor dials are specified.";
  }

  if (DipoleNormCCQEIsUsed && !MaCCQEIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormCCQE was "
           "requested but MaCCQE was not specified to be shape-only.";
  }

  if (IsDipoleReWeight) {
    if (DipoleNormCCQEIsUsed) {
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(cfg, "NormCCQE", param,
                                                 firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
    if (DipoleMaCCQEIsUsed) {
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(cfg, "MaCCQE", param,
                                                 firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
  } else if (IsZExpReWeight) {
    // ZNorm enters in linearly to the weight calculation so doesn't need to
    // output its response via the a meta-parameter.
    if (ZNormIsUsed) {
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(cfg, "ZNormCCQE", param,
                                                 firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
    if (ZExpA1IsUsed || ZExpA2IsUsed || ZExpA3IsUsed || ZExpA4IsUsed) {
      SystParamHeader ZExp;
      if (!ignore_parameter_dependence) {
        ZExp.prettyName = "ZExpAVariationResponse";
        ZExp.systParamId = firstParamId++;
      }

      if (ZExpA1IsUsed) {
        systtools::SystParamHeader param;
        ParseFHiCLSimpleToolConfigurationParameter(cfg, "ZExpA1", param,
                                                   firstParamId);
        param.systParamId = firstParamId++;
        if (!ignore_parameter_dependence) {
          param.isResponselessParam = true;
          param.responseParamId = ZExp.systParamId;
          if (param.isSplineable) {
            throw invalid_ToolConfigurationFHiCL()
                << "[ERROR]: Attempted to build spline from "
                   "parameter ZExpA1 , which enters into an intrinsically "
                   "multi-parameter response calculation. Either run in random "
                   "throw mode or set \"ignore_parameter_dependence\" in the "
                   "GENIEReWeight_tool configuration.";
          }
        }
        QEmd.push_back(std::move(param));
      }
      if (ZExpA2IsUsed) {
        systtools::SystParamHeader param;
        ParseFHiCLSimpleToolConfigurationParameter(cfg, "ZExpA2", param,
                                                   firstParamId);
        param.systParamId = firstParamId++;
        if (!ignore_parameter_dependence) {
          param.isResponselessParam = true;
          param.responseParamId = ZExp.systParamId;
          if (param.isSplineable) {
            throw invalid_ToolConfigurationFHiCL()
                << "[ERROR]: Attempted to build spline from "
                   "parameter ZExpA2 , which enters into an intrinsically "
                   "multi-parameter response calculation. Either run in random "
                   "throw mode or set \"ignore_parameter_dependence\" in the "
                   "GENIEReWeight_tool configuration.";
          }
        }
        QEmd.push_back(std::move(param));
      }
      if (ZExpA3IsUsed) {
        systtools::SystParamHeader param;
        ParseFHiCLSimpleToolConfigurationParameter(cfg, "ZExpA3", param,
                                                   firstParamId);
        param.systParamId = firstParamId++;
        if (!ignore_parameter_dependence) {
          param.isResponselessParam = true;
          param.responseParamId = ZExp.systParamId;
          if (param.isSplineable) {
            throw invalid_ToolConfigurationFHiCL()
                << "[ERROR]: Attempted to build spline from "
                   "parameter ZExpA3 , which enters into an intrinsically "
                   "multi-parameter response calculation. Either run in random "
                   "throw mode or set \"ignore_parameter_dependence\" in the "
                   "GENIEReWeight_tool configuration.";
          }
        }
        QEmd.push_back(std::move(param));
      }
      if (ZExpA4IsUsed) {
        systtools::SystParamHeader param;
        ParseFHiCLSimpleToolConfigurationParameter(cfg, "ZExpA4", param,
                                                   firstParamId);
        param.systParamId = firstParamId++;
        if (!ignore_parameter_dependence) {
          param.isResponselessParam = true;
          param.responseParamId = ZExp.systParamId;
          if (param.isSplineable) {
            throw invalid_ToolConfigurationFHiCL()
                << "[ERROR]: Attempted to build spline from "
                   "parameter ZExpA4 , which enters into an intrinsically "
                   "multi-parameter response calculation. Either run in random "
                   "throw mode or set \"ignore_parameter_dependence\" in the "
                   "GENIEReWeight_tool configuration.";
          }
        }
        QEmd.push_back(std::move(param));
      }
      if (!ignore_parameter_dependence) {
        QEmd.push_back(std::move(ZExp));
      }
    }
  }

  bool VecFFCCQEIsBBA = cfg.get<bool>("VecFFCCQEIsBBA", false);

  if (!VecFFCCQEIsBBA) {
    SystParamHeader vecFFQE;
    vecFFQE.prettyName = "VecFFCCQEshape";
    vecFFQE.systParamId = firstParamId++;
    vecFFQE.isCorrection = true;
    vecFFQE.centralParamValue = 1;
    QEmd.push_back(std::move(vecFFQE));
  }

  bool AxFFCCQEDipoleToZExp =
      cfg.get<bool>("AxFFCCQEDipoleToZExp", false) || IsZExpReWeight;

  if (AxFFCCQEDipoleToZExp) {
    SystParamHeader axFFQE;
    axFFQE.prettyName = "AxFFCCQEshape";
    axFFQE.systParamId = firstParamId++;
    axFFQE.isCorrection = true;
    axFFQE.centralParamValue = 1;
    QEmd.push_back(std::move(axFFQE));
  }

  if (HasParam(QEmd, "ZExpAVariationResponse")) {
    FinalizeAndValidateDependentParameters(
        QEmd, "ZExpAVariationResponse",
        {"ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE"});
  }

  return QEmd;
}

#ifndef GRWTEST

SystMetaData ConfigureNCELParameterHeaders(fhicl::ParameterSet const &cfg,
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
    if (EtaNCELIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(EtaNCEL, NCELmd, NCELResp.systParamId);
    }
    NCELmd.push_back(std::move(NCELResp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : NCELmd) {
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

SystMetaData ConfigureRESParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId) {
  SystMetaData RESmd;

  //************* CCRES
  bool NormCCRESIsUsed = PARAM_IS_USED_BY_CFG(NormCCRES);
  bool CCRESIsShapeOnly = cfg().CCRESIsShapeOnly() || NormCCRESIsUsed;
  bool MaCCRESIsUsed = PARAM_IS_USED_BY_CFG(MaCCRES);
  bool MvCCRESIsUsed = PARAM_IS_USED_BY_CFG(MvCCRES);

  if (NormCCRESIsUsed) {
    ADD_PARAM_TO_SYST(NormCCRES, RESmd);
  }

  if (MaCCRESIsUsed || MvCCRESIsUsed) {
    SystParamHeader CCRESResp;
    CCRESResp.prettyName = "CCRESVariationResponse";
    CCRESResp.systParamId = firstParamId++;

    if (MaCCRESIsUsed) {
      ADD_PARAM_TO_SYST(MaCCRES, RESmd);
      if (CCRESIsShapeOnly) {
        GetParam(RESmd, "MaCCRES").opts.push_back("shape");
      }
    }
    if (MvCCRESIsUsed) {
      ADD_PARAM_TO_SYST(MvCCRES, RESmd);
      if (CCRESIsShapeOnly) {
        GetParam(RESmd, "MvCCRES").opts.push_back("shape");
      }
    }
    RESmd.push_back(std::move(CCRESResp));
  }

  //************* NCRES
  bool NormNCRESIsUsed = PARAM_IS_USED_BY_CFG(NormNCRES);
  bool NCRESIsShapeOnly = cfg().NCRESIsShapeOnly() || NormNCRESIsUsed;
  bool MaNCRESIsUsed = PARAM_IS_USED_BY_CFG(MaNCRES);
  bool MvNCRESIsUsed = PARAM_IS_USED_BY_CFG(MvNCRES);

  if (NormNCRESIsUsed) {
    ADD_PARAM_TO_SYST(NormNCRES, RESmd);
  }

  if (MaNCRESIsUsed || MvNCRESIsUsed) {
    SystParamHeader NCRESResp;
    NCRESResp.prettyName = "NCRESVariationResponse";
    NCRESResp.systParamId = firstParamId++;

    if (MaNCRESIsUsed) {
      ADD_PARAM_TO_SYST(MaNCRES, RESmd);
      if (NCRESIsShapeOnly) {
        GetParam(RESmd, "MaNCRES").opts.push_back("shape");
      }
    }
    if (MvNCRESIsUsed) {
      ADD_PARAM_TO_SYST(MvNCRES, RESmd);
      if (NCRESIsShapeOnly) {
        GetParam(RESmd, "MvNCRES").opts.push_back("shape");
      }
    }
    RESmd.push_back(std::move(NCRESResp));
  }

  // These are all independent and based upon the channel that was generated
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvpCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvpCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvpNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvpNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvnCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvnCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvnNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvnNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarpCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarpCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarpNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarpNC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarnCC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarnCC2pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarnNC1pi, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(NonRESBGvbarnNC2pi, RESmd);

  CHECK_USED_ADD_PARAM_TO_SYST(RDecBR1gamma, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(RDecBR1eta, RESmd);
  CHECK_USED_ADD_PARAM_TO_SYST(Theta_Delta2Npi, RESmd);

  std::unique_ptr<CLHEP::MTwistEngine> RNgine =
      std::make_unique<CLHEP::MTwistEngine>(0);
  std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
      std::make_unique<CLHEP::RandGaussQ>(*RNgine);

  for (auto &hdr : RESmd) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }

  if (HasParam(RESmd, "CCRESVariationResponse")) {
    uint64_t NCCRESVariations = 0;

    auto const &CCRESParamNames = {"MaCCRES", "MvCCRES"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(CCRESParamNames, RESmd,
                                               NCCRESVariations);
    std::vector<double> CCRESdummyParamVars;
    for (size_t i = 0; i < NCCRESVariations; ++i) {
      CCRESdummyParamVars.push_back(i);
    }

    GetParam(RESmd, "CCRESVariationResponse").paramVariations =
        std::move(CCRESdummyParamVars);
  }

  if (HasParam(RESmd, "NCRESVariationResponse")) {

    uint64_t NNCRESVariations = 0;

    auto const &NCRESParamNames = {"MaNCRES", "MvNCRES"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(NCRESParamNames, RESmd,
                                               NNCRESVariations);
    std::vector<double> NCRESdummyParamVars;
    for (size_t i = 0; i < NNCRESVariations; ++i) {
      NCRESdummyParamVars.push_back(i);
    }

    GetParam(RESmd, "NCRESVariationResponse").paramVariations =
        std::move(NCRESdummyParamVars);
  }

  return RESmd;
}

SystMetaData ConfigureCOHParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId) {
  SystMetaData COHmd;

  bool MaCOHpiIsUsed = PARAM_IS_USED_BY_CFG(MaCOHpi);
  bool R0COHpiIsUsed = PARAM_IS_USED_BY_CFG(R0COHpi);

  if (MaCOHpiIsUsed || R0COHpiIsUsed) {
    SystParamHeader COHResp;
    COHResp.prettyName = "COHVariationResponse";
    COHResp.systParamId = firstParamId++;

    if (MaCOHpiIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(MaCOHpi, COHmd, COHResp.systParamId);
    }
    if (R0COHpiIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(R0COHpi, COHmd, COHResp.systParamId);
    }
    COHmd.push_back(std::move(COHResp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : COHmd) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &COHParamNames = {"MaCOHpi", "R0COHpi"};
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

SystMetaData ConfigureDISParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId) {
  SystMetaData DISmd;

  bool DISIsShapeOnly = cfg().DISIsShapeOnly();

  bool AhtBYIsUsed = PARAM_IS_USED_BY_CFG(AhtBY);
  bool BhtBYIsUsed = PARAM_IS_USED_BY_CFG(BhtBY);
  bool CV1uBYIsUsed = PARAM_IS_USED_BY_CFG(CV1uBY);
  bool CV2uBYIsUsed = PARAM_IS_USED_BY_CFG(CV2uBY);

  std::unique_ptr<CLHEP::MTwistEngine> RNgine =
      std::make_unique<CLHEP::MTwistEngine>(0);
  std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
      std::make_unique<CLHEP::RandGaussQ>(*RNgine);

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
    DISmd.push_back(std::move(DISResp));

    for (auto &hdr : DISmd) {
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

  bool AGKY_xF1piIsUsed = PARAM_IS_USED_BY_CFG(AGKY_xF1pi);
  bool AGKY_pT1piIsUsed = PARAM_IS_USED_BY_CFG(AGKY_pT1pi);
  if (AGKY_xF1piIsUsed || AGKY_pT1piIsUsed) {
    SystParamHeader AGKYResp;
    AGKYResp.prettyName = "AGKYVariationResponse";
    AGKYResp.systParamId = firstParamId++;
    if (AGKY_xF1piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(AGKY_xF1pi, DISmd, AGKYResp.systParamId);
    }
    if (AGKY_pT1piIsUsed) {
      ADD_PARAM_TO_SYST_RESPONSELESS(AGKY_pT1pi, DISmd, AGKYResp.systParamId);
    }
    DISmd.push_back(std::move(AGKYResp));

    for (auto &hdr : DISmd) {
      MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
    }

    uint64_t NVariations = 0;

    auto const &AGKYParamNames = {"AGKY_xF1pi", "AGKY_pT1pi"};
    GET_NVARIATIONS_OF_RESPONSELESS_PARAMETERS(AGKYParamNames, DISmd,
                                               NVariations);
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(DISmd, "AGKYVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  CHECK_USED_ADD_PARAM_TO_SYST(FormZone, DISmd);

  for (auto &hdr : DISmd) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }

  return DISmd;
}

SystMetaData ConfigureFSIParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId) {
  SystMetaData FSImd;

  bool MFP_piIsUsed = PARAM_IS_USED_BY_CFG(MFP_pi);
  bool FrCEx_piIsUsed = PARAM_IS_USED_BY_CFG(FrCEx_pi);
  bool FrElas_piIsUsed = PARAM_IS_USED_BY_CFG(FrElas_pi);
  bool FrInel_piIsUsed = PARAM_IS_USED_BY_CFG(FrInel_pi);
  bool FrAbs_piIsUsed = PARAM_IS_USED_BY_CFG(FrAbs_pi);
  bool FrPiProd_piIsUsed = PARAM_IS_USED_BY_CFG(FrPiProd_pi);

  if (MFP_piIsUsed || FrCEx_piIsUsed || FrElas_piIsUsed || FrInel_piIsUsed ||
      FrAbs_piIsUsed || FrPiProd_piIsUsed) {
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
    FSImd.push_back(std::move(FSI_pi_Resp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : FSImd) {
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

  if (MFP_NIsUsed || FrCEx_NIsUsed || FrElas_NIsUsed || FrInel_NIsUsed ||
      FrAbs_NIsUsed || FrPiProd_NIsUsed) {
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
    FSImd.push_back(std::move(FSI_N_Resp));

    std::unique_ptr<CLHEP::MTwistEngine> RNgine =
        std::make_unique<CLHEP::MTwistEngine>(0);
    std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
        std::make_unique<CLHEP::RandGaussQ>(*RNgine);

    for (auto &hdr : FSImd) {
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

SystMetaData ConfigureOtherParameterHeaders(fhicl::ParameterSet const &cfg,
                                            paramId_t firstParamId) {
  SystMetaData Othermd;

  CHECK_USED_ADD_PARAM_TO_SYST(CCQEPauliSupViaKF, Othermd);
  CHECK_USED_ADD_PARAM_TO_SYST(CCQEMomDistroFGtoSF, Othermd);

  size_t CCQEMomDistroFGtoSFIndex =
      GetParamIndex(Othermd, "CCQEMomDistroFGtoSF");
  if (IndexIsHandled(Othermd, CCQEMomDistroFGtoSFIndex)) {
    Othermd[CCQEMomDistroFGtoSFIndex].paramValidityRange[0] = 0;
    Othermd[CCQEMomDistroFGtoSFIndex].paramValidityRange[1] = 1;
  }

  std::unique_ptr<CLHEP::MTwistEngine> RNgine =
      std::make_unique<CLHEP::MTwistEngine>(0);
  std::unique_ptr<CLHEP::RandGaussQ> RNJesus =
      std::make_unique<CLHEP::RandGaussQ>(*RNgine);

  for (auto &hdr : Othermd) {
    MakeThrowsIfNeeded(hdr, RNJesus, cfg().numberOfThrows());
  }

  return Othermd;
}

#endif

} // namespace nusyst
