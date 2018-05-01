#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN

#include "larsyst/interface/SystMetaData.hh"
#include "larsyst/interface/types.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

#include <string>

#define NOMINAL_HELP_TEXT(DIAL_NAME)                                              \
#DIAL_NAME                                                                      \
      " Nominal value, set to 0 to use in its nominal (GENIE-defined)"            \
      "state, If unitsAreNatural == true, then this value will always be passed " \
      "to the weight engine."

#define TWEAK_DEFINITION_HELP_TEXT(DIAL_NAME)                                  \
  "Defines variations in " #DIAL_NAME ". If not specified, but "               \
  "NominalValue is, then a single response at the new nominal value will be "  \
  "generated. If specified as \"{sigma_low, sigma_up}\", then normally "       \
  "distributed throws will be made with those as the 2*half gaussian widths "  \
  "above and below the central value. If the list is specified as "            \
  "\"(-1,2,3,...)\", then spline knots will be made at specified values, "     \
  "respects unitsAreNatural. If specified as \"[1,2,3]\", then discrete "      \
  "non-splineable tweaks will be made at those values, respects "              \
  "unitsAreNatural."

#define DIAL_NAME_HELPER(DIAL_NAME, b) DIAL_NAME##b
#define TWEAKABLE_PARAMETER_DEFINITION(DIAL_NAME)                              \
  fhicl::Atom<double> DIAL_NAME_HELPER(DIAL_NAME, NominalValue){               \
      fhicl::Name(#DIAL_NAME "Nominal"),                                       \
      fhicl::Comment(NOMINAL_HELP_TEXT(DIAL_NAME)), 0xdeadb33f};               \
  fhicl::Atom<std::string> DIAL_NAME_HELPER(DIAL_NAME, TweakDefinition) {      \
    fhicl::Name(#DIAL_NAME "TweakDefinition"),                                 \
        fhicl::Comment(TWEAK_DEFINITION_HELP_TEXT(DIAL_NAME)), ""              \
  }

namespace nusyst {

struct GENIEReWeightParamConfig {
  //********* START Engine-level options ***********
  fhicl::Atom<bool> unitsAreNatural{
      fhicl::Name("unitsAreNatural"),
      fhicl::Comment(
          "Whether parameter nominal values should be interpreted as in "
          "GENIE-defined sigma units, or natural units. Default == false"),
      []() { return false; }, false};
  fhicl::Atom<uint64_t> numberOfThrows{
      fhicl::Name("numberOfThrows"),
      fhicl::Comment("Number of throws to make for parameters "), 0};
  //********* END Engine-level options ***********

  //********* START QE options ***********
  fhicl::Atom<bool> MAQEIsShapeOnly{
      fhicl::Name("MAQEIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in MAQE only affect the shape of the event rate "
          "distribution and not the total normalization. Default == false"),
      false};

  TWEAKABLE_PARAMETER_DEFINITION(NormCCQE);
  TWEAKABLE_PARAMETER_DEFINITION(MaCCQE);

  fhicl::Atom<bool> VecFFCCQEIsBBA{
      fhicl::Name("VecFFCCQEIsBBA"),
      fhicl::Comment("Use BBA vector form factors, otherwise "
                     "use dipole. {Default == true}"),
      true};

  TWEAKABLE_PARAMETER_DEFINITION(MaNCEL);
  TWEAKABLE_PARAMETER_DEFINITION(EtaNCEL);

  //********* START 1p1h Z-exp options ***********

  fhicl::Atom<int> AxFFCCQEDipoleToZExp{
      fhicl::Name("AxFFCCQEDipoleToZExp"),
      fhicl::Comment(
          "Whether to apply dipole to z-expansion reweighting. {Default == "
          "true if any z-expansion dials are used, false otherwise}"),
      -1};

  TWEAKABLE_PARAMETER_DEFINITION(ZNormCCQE);
  TWEAKABLE_PARAMETER_DEFINITION(ZExpA1CCQE);
  TWEAKABLE_PARAMETER_DEFINITION(ZExpA2CCQE);
  TWEAKABLE_PARAMETER_DEFINITION(ZExpA3CCQE);
  TWEAKABLE_PARAMETER_DEFINITION(ZExpA4CCQE);

  //********* END 1p1h Z-exp options ***********

  //********* END QE options ***********

  //********* START RES options ***********
  fhicl::Atom<bool> CCRESIsShapeOnly{
      fhicl::Name("CCRESIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in CCRES parameters only affect the shape of "
          "the event rate distribution and not the total normalization. "
          "Default == false"),
      false};

  TWEAKABLE_PARAMETER_DEFINITION(NormCCRES);
  TWEAKABLE_PARAMETER_DEFINITION(MaCCRES);
  TWEAKABLE_PARAMETER_DEFINITION(MvCCRES);

  fhicl::Atom<bool> NCRESIsShapeOnly{
      fhicl::Name("NCRESIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in NCRES parameters only affect the shape of "
          "the event rate distribution and not the total normalization. "
          "Default == false"),
      false};

  TWEAKABLE_PARAMETER_DEFINITION(NormNCRES);
  TWEAKABLE_PARAMETER_DEFINITION(MaNCRES);
  TWEAKABLE_PARAMETER_DEFINITION(MvNCRES);

  //********* START Background pi-production ***********

  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvpCC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvpCC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvpNC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvpNC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvnCC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvnCC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvnNC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvnNC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarpCC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarpCC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarpNC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarpNC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarnCC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarnCC2pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarnNC1pi);
  TWEAKABLE_PARAMETER_DEFINITION(NonRESBGvbarnNC2pi);

  //********* END Background pi-production ***********

  // Resonance decay kinematics options
  TWEAKABLE_PARAMETER_DEFINITION(RDecBR1gamma);
  TWEAKABLE_PARAMETER_DEFINITION(RDecBR1eta);
  TWEAKABLE_PARAMETER_DEFINITION(Theta_Delta2Npi);

  //********* END RES options ***********

  //********* START COH options ***********
  TWEAKABLE_PARAMETER_DEFINITION(MaCOHpi);
  TWEAKABLE_PARAMETER_DEFINITION(R0COHpi);
  //********* END COH options ***********

  //********* START DIS options ***********
  fhicl::Atom<bool> DISIsShapeOnly{
      fhicl::Name("DISIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in Bodek-Yang DIS parameters only affect the "
          "shape of the event rate distribution and not the total "
          "normalization. Default == false"),
      false};

  TWEAKABLE_PARAMETER_DEFINITION(AhtBY);
  TWEAKABLE_PARAMETER_DEFINITION(BhtBY);
  TWEAKABLE_PARAMETER_DEFINITION(CV1uBY);
  TWEAKABLE_PARAMETER_DEFINITION(CV2uBY);

  TWEAKABLE_PARAMETER_DEFINITION(AGKY_xF1pi);
  TWEAKABLE_PARAMETER_DEFINITION(AGKY_pT1pi);

  TWEAKABLE_PARAMETER_DEFINITION(FormZone);

  //********* END DIS options ***********

  //********* START FSI options ***********

  TWEAKABLE_PARAMETER_DEFINITION(MFP_pi);
  TWEAKABLE_PARAMETER_DEFINITION(FrCEx_pi);
  TWEAKABLE_PARAMETER_DEFINITION(FrElas_pi);
  TWEAKABLE_PARAMETER_DEFINITION(FrInel_pi);
  TWEAKABLE_PARAMETER_DEFINITION(FrAbs_pi);
  TWEAKABLE_PARAMETER_DEFINITION(FrPiProd_pi);

  TWEAKABLE_PARAMETER_DEFINITION(MFP_N);
  TWEAKABLE_PARAMETER_DEFINITION(FrCEx_N);
  TWEAKABLE_PARAMETER_DEFINITION(FrElas_N);
  TWEAKABLE_PARAMETER_DEFINITION(FrInel_N);
  TWEAKABLE_PARAMETER_DEFINITION(FrAbs_N);
  TWEAKABLE_PARAMETER_DEFINITION(FrPiProd_N);

  //********* END FSI options ***********

  //********* START Other options ***********
  // QE nuclear model options
  TWEAKABLE_PARAMETER_DEFINITION(CCQEPauliSupViaKF);
  TWEAKABLE_PARAMETER_DEFINITION(CCQEMomDistroFGtoSF);

  //********* END Other options ***********
};

larsyst::SystMetaData
ConfigureQEParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                            larsyst::paramId_t);

larsyst::SystMetaData
ConfigureNCELParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                              larsyst::paramId_t);

larsyst::SystMetaData
ConfigureRESParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                             larsyst::paramId_t);

larsyst::SystMetaData
ConfigureCOHParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                             larsyst::paramId_t);

larsyst::SystMetaData
ConfigureDISParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                             larsyst::paramId_t);

larsyst::SystMetaData
ConfigureFSIParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                             larsyst::paramId_t);

larsyst::SystMetaData
ConfigureOtherParameterHeaders(fhicl::Table<GENIEReWeightParamConfig> const &,
                               larsyst::paramId_t);
} // namespace nusyst

#undef NOMINAL_HELP_TEXT
#undef TWEAK_DEFINITION_HELP_TEXT
#undef DIAL_NAME_HELPER
#undef TWEAKABLE_PARAMETER_DEFINITION

#endif
