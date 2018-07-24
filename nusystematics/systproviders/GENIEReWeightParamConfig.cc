#include "nusystematics/systproviders/GENIEReWeightParamConfig.hh"

#include "systematicstools/interface/ISystProvider_tool.hh"

#include "systematicstools/utility/ResponselessParamUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "ReWeight/GSyst.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace systtools;
using namespace genie::rew;

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
  bool DipoleNormCCQEIsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormCCQE));
  bool DipoleIsShapeOnly = MaCCQEIsShapeOnly || DipoleNormCCQEIsUsed;
  bool DipoleMaCCQEIsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_MaCCQE));

  bool IsDipoleReWeight =
      DipoleIsShapeOnly || DipoleNormCCQEIsUsed || DipoleMaCCQEIsUsed;

  bool ZNormIsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZNormCCQE));
  bool ZExpA1IsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA1CCQE));
  bool ZExpA2IsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA2CCQE));
  bool ZExpA3IsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA3CCQE));
  bool ZExpA4IsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_ZExpA4CCQE));

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
      ParseFHiCLSimpleToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_NormCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
    if (DipoleMaCCQEIsUsed) {
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_MaCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      param.prettyName = GSyst::AsString(
          DipoleIsShapeOnly ? kXSecTwkDial_MaCCQEshape : kXSecTwkDial_MaCCQE);
      QEmd.push_back(std::move(param));
    }
  } else if (IsZExpReWeight) {
    // ZNorm enters in linearly to the weight calculation so doesn't need to
    // output its response via the a meta-parameter.
    if (ZNormIsUsed) {
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(
          cfg, GSyst::AsString(kXSecTwkDial_ZNormCCQE), param, firstParamId);
      param.systParamId = firstParamId++;
      QEmd.push_back(std::move(param));
    }
    SystParamHeader ZExp;
    if (!ignore_parameter_dependence) {
      ZExp.prettyName = "ZExpAVariationResponse";
      ZExp.systParamId = firstParamId++;
    }
    for (GSyst_t gdial : {kXSecTwkDial_ZExpA1CCQE, kXSecTwkDial_ZExpA2CCQE,
                          kXSecTwkDial_ZExpA3CCQE, kXSecTwkDial_ZExpA4CCQE}) {
      if (!FHiCLSimpleToolConfigurationParameterExists(
              cfg, GSyst::AsString(gdial))) {
        continue;
      }
      systtools::SystParamHeader param;
      ParseFHiCLSimpleToolConfigurationParameter(cfg, GSyst::AsString(gdial),
                                                 param, firstParamId);
      param.systParamId = firstParamId++;
      if (!ignore_parameter_dependence) {
        param.isResponselessParam = true;
        param.responseParamId = ZExp.systParamId;
        if (param.isSplineable) {
          throw invalid_ToolConfigurationFHiCL()
              << "[ERROR]: Attempted to build spline from "
                 "parameter "
              << param.prettyName
              << ", which enters into an intrinsically "
                 "multi-parameter response calculation. Either run in "
                 "random "
                 "throw mode or set \"ignore_parameter_dependence\" in the "
                 "GENIEReWeight_tool configuration.";
        }
      }
      QEmd.push_back(std::move(param));
    }
    if (!ignore_parameter_dependence) {
      QEmd.push_back(std::move(ZExp));
      FinalizeAndValidateDependentParameters(
          QEmd, "ZExpAVariationResponse",
          {GSyst::AsString(kXSecTwkDial_ZExpA1CCQE),
           GSyst::AsString(kXSecTwkDial_ZExpA2CCQE),
           GSyst::AsString(kXSecTwkDial_ZExpA3CCQE),
           GSyst::AsString(kXSecTwkDial_ZExpA4CCQE)});
    }
  }

  bool VecFFCCQEIsBBA = cfg.get<bool>("VecFFCCQEIsBBA", false);

  if (!VecFFCCQEIsBBA) {
    SystParamHeader vecFFQE;
    vecFFQE.prettyName = GSyst::AsString(kXSecTwkDial_VecFFCCQEshape);
    vecFFQE.systParamId = firstParamId++;
    vecFFQE.isCorrection = true;
    vecFFQE.centralParamValue = 1;
    QEmd.push_back(std::move(vecFFQE));
  }

  bool AxFFCCQEDipoleToZExp =
      cfg.get<bool>("AxFFCCQEDipoleToZExp", false) || IsZExpReWeight;

  if (AxFFCCQEDipoleToZExp) {
    SystParamHeader axFFQE;
    axFFQE.prettyName = GSyst::AsString(kXSecTwkDial_AxFFCCQEshape);
    axFFQE.systParamId = firstParamId++;
    axFFQE.isCorrection = true;
    axFFQE.centralParamValue = 1;
    QEmd.push_back(std::move(axFFQE));
  }

  return QEmd;
} // namespace nusyst

SystMetaData ConfigureNCELParameterHeaders(fhicl::ParameterSet const &cfg,
                                           paramId_t firstParamId,
                                           fhicl::ParameterSet &tool_options) {
  SystMetaData NCELmd;

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystParamHeader NCELResp;
  if (!ignore_parameter_dependence) {
    NCELResp.prettyName = "NCELVariationResponse";
    NCELResp.systParamId = firstParamId++;
  }

  for (GSyst_t const &dial : {kXSecTwkDial_MaNCEL, kXSecTwkDial_EtaNCEL}) {
    std::string const &pname = GSyst::AsString(dial);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = NCELResp.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in random "
               "throw mode or set \"ignore_parameter_dependence\" in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    NCELmd.push_back(std::move(param));
  }

  if (!ignore_parameter_dependence && NCELmd.size()) {
    NCELmd.push_back(std::move(NCELResp));
    FinalizeAndValidateDependentParameters(
        NCELmd, "NCELVariationResponse",
        {GSyst::AsString(kXSecTwkDial_MaNCEL),
         GSyst::AsString(kXSecTwkDial_EtaNCEL)});
  }

  return NCELmd;
}

SystMetaData ConfigureRESParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {
  SystMetaData RESmd;

  //************* CCRES
  bool CCRESIsShapeOnly = cfg.get<bool>("CCRESIsShapeOnly", false);
  tool_options.put("CCRESIsShapeOnly", CCRESIsShapeOnly);

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);
  bool NormCCRESIsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormCCRES));

  if (NormCCRESIsUsed && !CCRESIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormCCRES was "
           "requested but CCRES was not specified to be shape-only.";
  }

  if (NormCCRESIsUsed) {
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(
        cfg, GSyst::AsString(kXSecTwkDial_NormCCRES), param, firstParamId);
    param.systParamId = firstParamId++;
    RESmd.push_back(std::move(param));
  }

  SystParamHeader CCRESresp;
  if (!ignore_parameter_dependence) {
    CCRESresp.prettyName = "CCResVariationResponse";
    CCRESresp.systParamId = firstParamId++;
  }

  size_t NCCResParams = 0;
  for (std::pair<GSyst_t, GSyst_t> const &name_dial :
       std::vector<std::pair<GSyst_t, GSyst_t>>{
           {{kXSecTwkDial_MaCCRES, kXSecTwkDial_MaCCRESshape},
            {kXSecTwkDial_MvCCRES, kXSecTwkDial_MvCCRESshape}}}) {
    std::string const &pname = GSyst::AsString(name_dial.first);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    param.prettyName = CCRESIsShapeOnly ? GSyst::AsString(name_dial.second)
                                        : GSyst::AsString(name_dial.first);

    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = CCRESresp.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << pname
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in random "
               "throw mode or set \"ignore_parameter_dependence\" in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    RESmd.push_back(std::move(param));
    NCCResParams++;
  }

  if (!ignore_parameter_dependence && NCCResParams) {
    RESmd.push_back(std::move(CCRESresp));
    FinalizeAndValidateDependentParameters(
        RESmd, "CCResVariationResponse",
        {GSyst::AsString(CCRESIsShapeOnly ? kXSecTwkDial_MaCCRESshape
                                          : kXSecTwkDial_MaCCRES),
         GSyst::AsString(CCRESIsShapeOnly ? kXSecTwkDial_MvCCRESshape
                                          : kXSecTwkDial_MvCCRES)});
  }

  //************* NCRES
  bool NCRESIsShapeOnly = cfg.get<bool>("NCRESIsShapeOnly", false);
  tool_options.put("NCRESIsShapeOnly", NCRESIsShapeOnly);

  bool NormNCRESIsUsed = FHiCLSimpleToolConfigurationParameterExists(
      cfg, GSyst::AsString(kXSecTwkDial_NormNCRES));

  if (NormNCRESIsUsed && !NCRESIsShapeOnly) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring GENIEReWeight_tool, NormNCRES was "
           "requested but NCRES was not specified to be shape-only.";
  }

  if (NormNCRESIsUsed) {
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(
        cfg, GSyst::AsString(kXSecTwkDial_NormNCRES), param, firstParamId);
    param.systParamId = firstParamId++;
    RESmd.push_back(std::move(param));
  }

  SystParamHeader NCRESresp;
  if (!ignore_parameter_dependence) {
    NCRESresp.prettyName = "NCResVariationResponse";
    NCRESresp.systParamId = firstParamId++;
  }

  size_t NNCResParams = 0;
  for (std::pair<GSyst_t, GSyst_t> const &name_dial :
       std::vector<std::pair<GSyst_t, GSyst_t>>{
           {{kXSecTwkDial_MaNCRES, kXSecTwkDial_MaNCRESshape},
            {kXSecTwkDial_MvNCRES, kXSecTwkDial_MvNCRESshape}}}) {
    std::string const &pname = GSyst::AsString(name_dial.first);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    param.prettyName = NCRESIsShapeOnly ? GSyst::AsString(name_dial.second)
                                        : GSyst::AsString(name_dial.first);
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = NCRESresp.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in random "
               "throw mode or set \"ignore_parameter_dependence\" in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    RESmd.push_back(std::move(param));
    NNCResParams++;
  }

  if (!ignore_parameter_dependence && NNCResParams) {
    RESmd.push_back(std::move(NCRESresp));
    FinalizeAndValidateDependentParameters(
        RESmd, "NCResVariationResponse",
        {GSyst::AsString(NCRESIsShapeOnly ? kXSecTwkDial_MaNCRESshape
                                          : kXSecTwkDial_MaNCRES),
         GSyst::AsString(NCRESIsShapeOnly ? kXSecTwkDial_MvNCRESshape
                                          : kXSecTwkDial_MvNCRES)});
  }

  // These are all independent and based upon the channel that was generated
  for (GSyst_t const &gdial : std::vector<GSyst_t>{
           {kXSecTwkDial_RvpCC1pi, kXSecTwkDial_RvpCC2pi, kXSecTwkDial_RvpNC1pi,
            kXSecTwkDial_RvpNC2pi, kXSecTwkDial_RvnCC1pi, kXSecTwkDial_RvnCC2pi,
            kXSecTwkDial_RvnNC1pi, kXSecTwkDial_RvnNC2pi,
            kXSecTwkDial_RvbarpCC1pi, kXSecTwkDial_RvbarpCC2pi,
            kXSecTwkDial_RvbarpNC1pi, kXSecTwkDial_RvbarpNC2pi,
            kXSecTwkDial_RvbarnCC1pi, kXSecTwkDial_RvbarnCC2pi,
            kXSecTwkDial_RvbarnNC1pi, kXSecTwkDial_RvbarnNC2pi,

            kRDcyTwkDial_BR1gamma, kRDcyTwkDial_BR1eta,
            kRDcyTwkDial_Theta_Delta2Npi}}) {
    std::string const &pname = GSyst::AsString(gdial);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    RESmd.push_back(std::move(param));
  }

  return RESmd;
}

SystMetaData ConfigureCOHParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {
  SystMetaData COHmd;

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystParamHeader COHResp;
  if (!ignore_parameter_dependence) {
    COHResp.prettyName = "COHVariationResponse";
    COHResp.systParamId = firstParamId++;
  }

  for (GSyst_t const &dial : {kXSecTwkDial_MaCOHpi, kXSecTwkDial_R0COHpi}) {
    std::string const &pname = GSyst::AsString(dial);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = COHResp.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in random "
               "throw mode or set \"ignore_parameter_dependence\" in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    COHmd.push_back(std::move(param));
  }

  if (!ignore_parameter_dependence && COHmd.size()) {
    COHmd.push_back(std::move(COHResp));
    FinalizeAndValidateDependentParameters(
        COHmd, "COHVariationResponse",
        {GSyst::AsString(kXSecTwkDial_MaCOHpi),
         GSyst::AsString(kXSecTwkDial_R0COHpi)});
  }

  return COHmd;
}

SystMetaData ConfigureDISParameterHeaders(fhicl::ParameterSet const &cfg,
                                          paramId_t firstParamId,
                                          fhicl::ParameterSet &tool_options) {
  SystMetaData DISmd;

  bool DISBYIsShapeOnly = cfg.get<bool>("DISBYIsShapeOnly", false);
  tool_options.put("DISBYIsShapeOnly", DISBYIsShapeOnly);

  bool ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystParamHeader DISBYResponse;
  if (!ignore_parameter_dependence) {
    DISBYResponse.prettyName = "DISBYVariationResponse";
    DISBYResponse.systParamId = firstParamId++;
  }

  size_t NDISBYParams = 0;
  for (std::pair<GSyst_t, GSyst_t> const &name_dial :
       std::vector<std::pair<GSyst_t, GSyst_t>>{
           {{kXSecTwkDial_AhtBY, kXSecTwkDial_AhtBYshape},
            {kXSecTwkDial_BhtBY, kXSecTwkDial_BhtBYshape},
            {kXSecTwkDial_CV1uBY, kXSecTwkDial_CV1uBYshape},
            {kXSecTwkDial_CV2uBY, kXSecTwkDial_CV2uBYshape}}}) {
    std::string const &pname = GSyst::AsString(name_dial.first);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    param.prettyName = DISBYIsShapeOnly ? GSyst::AsString(name_dial.second)
                                        : GSyst::AsString(name_dial.first);
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = DISBYResponse.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in "
               "random "
               "throw mode or set \"ignore_parameter_dependence\" in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    DISmd.push_back(std::move(param));
    NDISBYParams++;
  }
  if (!ignore_parameter_dependence && NDISBYParams) {
    DISmd.push_back(std::move(DISBYResponse));
    FinalizeAndValidateDependentParameters(
        DISmd, "DISBYVariationResponse",
        {GSyst::AsString(DISBYIsShapeOnly ? kXSecTwkDial_AhtBYshape
                                          : kXSecTwkDial_AhtBY),
         GSyst::AsString(DISBYIsShapeOnly ? kXSecTwkDial_BhtBYshape
                                          : kXSecTwkDial_BhtBY),
         GSyst::AsString(DISBYIsShapeOnly ? kXSecTwkDial_CV1uBYshape
                                          : kXSecTwkDial_CV1uBY),
         GSyst::AsString(DISBYIsShapeOnly ? kXSecTwkDial_CV2uBYshape
                                          : kXSecTwkDial_CV2uBY)});
  }

  SystParamHeader AGKYResponse;
  if (!ignore_parameter_dependence) {
    AGKYResponse.prettyName = "AGKYVariationResponse";
    AGKYResponse.systParamId = firstParamId++;
  }

  size_t NAGKYParams = 0;
  for (GSyst_t const &gdial :
       {kHadrAGKYTwkDial_xF1pi, kHadrAGKYTwkDial_pT1pi}) {
    std::string const &pname = GSyst::AsString(gdial);
    if (!FHiCLSimpleToolConfigurationParameterExists(cfg, pname)) {
      continue;
    }
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, param, firstParamId);
    param.systParamId = firstParamId++;
    if (!ignore_parameter_dependence) {
      param.isResponselessParam = true;
      param.responseParamId = AGKYResponse.systParamId;
      if (param.isSplineable) {
        throw invalid_ToolConfigurationFHiCL()
            << "[ERROR]: Attempted to build spline from "
               "parameter "
            << param.prettyName
            << ", which enters into an intrinsically "
               "multi-parameter response calculation. Either run in "
               "random "
               "throw mode or set \"ignore_parameter_dependence\" "
               "in the "
               "GENIEReWeight_tool configuration.";
      }
    }
    DISmd.push_back(std::move(param));
    NAGKYParams++;
  }
  if (!ignore_parameter_dependence && NAGKYParams) {
    DISmd.push_back(std::move(AGKYResponse));
    FinalizeAndValidateDependentParameters(
        DISmd, "AGKYVariationResponse",
        {GSyst::AsString(kHadrAGKYTwkDial_xF1pi),
         GSyst::AsString(kHadrAGKYTwkDial_pT1pi)});
  }

  std::string const &fz_pname = GSyst::AsString(kHadrNuclTwkDial_FormZone);
  if (FHiCLSimpleToolConfigurationParameterExists(cfg, fz_pname)) {
    systtools::SystParamHeader param;
    ParseFHiCLSimpleToolConfigurationParameter(cfg, fz_pname, param,
                                               firstParamId);
    param.systParamId = firstParamId++;
    DISmd.push_back(std::move(param));
  }

  return DISmd;
}

#ifndef GRWTEST

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
