#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHTCONFIG_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHTCONFIG_SEEN

#include "larsyst/interface/SystMetaData.hh"
#include "larsyst/interface/types.hh"
#include "larsyst/utility/string_parsers.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEaxial.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

struct GENIEReWeightConfig {
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

  fhicl::Atom<double> NormCCQENominalValue{
      fhicl::Name("NormCCQENominal"),
      fhicl::Comment(
          "NormCCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> NormCCQETweakDefinition{
      fhicl::Name("NormCCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in NormCCQE. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> MaCCQENominalValue{
      fhicl::Name("MaCCQENominal"),
      fhicl::Comment(
          "MaCCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> MaCCQETweakDefinition{
      fhicl::Name("MaCCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in MaCCQE. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<bool> VecFFCCQEIsBBA{
      fhicl::Name("VecFFCCQEIsBBA"),
      fhicl::Comment("Use BBA vector form factors, otherwise "
                     "use dipole. {Default == true}"),
      true};

  fhicl::Atom<double> MaNCELNominalValue{
      fhicl::Name("MaNCELNominal"),
      fhicl::Comment(
          "MaNCEL Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> MaNCELTweakDefinition{
      fhicl::Name("MaNCELTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MaNCEL. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> EtaNCELNominalValue{
      fhicl::Name("EtaNCELNominal"),
      fhicl::Comment(
          "EtaNCEL Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> EtaNCELTweakDefinition{
      fhicl::Name("EtaNCELTweakDefinition"),
      fhicl::Comment(
          "Defines variations in EtaNCEL. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  //********* START 1p1h Z-exp options ***********
  fhicl::Atom<double> ZNormCCQENominalValue{
      fhicl::Name("ZNormCCQENominal"),
      fhicl::Comment(
          "ZNormCCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> ZNormCCQETweakDefinition{
      fhicl::Name("ZNormCCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in ZNormCCQE. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> ZExpA1CCQENominalValue{
      fhicl::Name("ZExpA1CCQENominal"),
      fhicl::Comment(
          "ZExpA1CCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> ZExpA1CCQETweakDefinition{
      fhicl::Name("ZExpA1CCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in ZExpA1CCQE. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> ZExpA2CCQENominalValue{
      fhicl::Name("ZExpA2CCQENominal"),
      fhicl::Comment(
          "ZExpA2CCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> ZExpA2CCQETweakDefinition{
      fhicl::Name("ZExpA2CCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in ZExpA2CCQE. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> ZExpA3CCQENominalValue{
      fhicl::Name("ZExpA3CCQENominal"),
      fhicl::Comment(
          "ZExpA3CCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> ZExpA3CCQETweakDefinition{
      fhicl::Name("ZExpA3CCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in ZExpA3CCQE. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> ZExpA4CCQENominalValue{
      fhicl::Name("ZExpA4CCQENominal"),
      fhicl::Comment(
          "ZExpA4CCQE Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};

  fhicl::Atom<std::string> ZExpA4CCQETweakDefinition{
      fhicl::Name("ZExpA4CCQETweakDefinition"),
      fhicl::Comment(
          "Defines variations in ZExpA4CCQE. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};
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

  fhicl::Atom<bool> NCRESIsShapeOnly{
      fhicl::Name("NCRESIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in NCRES parameters only affect the shape of "
          "the event rate distribution and not the total normalization. "
          "Default == false"),
      false};

  fhicl::Atom<double> NormNCRESNominalValue{
      fhicl::Name("NormNCRESNominal"),
      fhicl::Comment(
          "NormNCRES Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> NormNCRESTweakDefinition{
      fhicl::Name("NormNCRESTweakDefinition"),
      fhicl::Comment(
          "Defines variations in NormNCRES. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> MaNCRESNominalValue{
      fhicl::Name("MaNCRESNominal"),
      fhicl::Comment(
          "MaNCRES Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> MaNCRESTweakDefinition{
      fhicl::Name("MaNCRESTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MaNCRES. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> MvNCRESNominalValue{
      fhicl::Name("MvNCRESNominal"),
      fhicl::Comment(
          "MvNCRES Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> MvNCRESTweakDefinition{
      fhicl::Name("MvNCRESTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MvNCRES. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvpCC1piNominalValue{
      fhicl::Name("RvpCC1piNominal"),
      fhicl::Comment(
          "RvpCC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvpCC1piTweakDefinition{
      fhicl::Name("RvpCC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvpCC1pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvpCC2piNominalValue{
      fhicl::Name("RvpCC2piNominal"),
      fhicl::Comment(
          "RvpCC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvpCC2piTweakDefinition{
      fhicl::Name("RvpCC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvpCC2pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvpNC1piNominalValue{
      fhicl::Name("RvpNC1piNominal"),
      fhicl::Comment(
          "RvpNC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvpNC1piTweakDefinition{
      fhicl::Name("RvpNC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvpNC1pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvpNC2piNominalValue{
      fhicl::Name("RvpNC2piNominal"),
      fhicl::Comment(
          "RvpNC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvpNC2piTweakDefinition{
      fhicl::Name("RvpNC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvpNC2pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvnCC1piNominalValue{
      fhicl::Name("RvnCC1piNominal"),
      fhicl::Comment(
          "RvnCC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvnCC1piTweakDefinition{
      fhicl::Name("RvnCC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvnCC1pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvnCC2piNominalValue{
      fhicl::Name("RvnCC2piNominal"),
      fhicl::Comment(
          "RvnCC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvnCC2piTweakDefinition{
      fhicl::Name("RvnCC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvnCC2pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvnNC1piNominalValue{
      fhicl::Name("RvnNC1piNominal"),
      fhicl::Comment(
          "RvnNC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvnNC1piTweakDefinition{
      fhicl::Name("RvnNC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvnNC1pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvnNC2piNominalValue{
      fhicl::Name("RvnNC2piNominal"),
      fhicl::Comment(
          "RvnNC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvnNC2piTweakDefinition{
      fhicl::Name("RvnNC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvnNC2pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarpCC1piNominalValue{
      fhicl::Name("RvbarpCC1piNominal"),
      fhicl::Comment(
          "RvbarpCC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarpCC1piTweakDefinition{
      fhicl::Name("RvbarpCC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarpCC1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarpCC2piNominalValue{
      fhicl::Name("RvbarpCC2piNominal"),
      fhicl::Comment(
          "RvbarpCC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarpCC2piTweakDefinition{
      fhicl::Name("RvbarpCC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarpCC2pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarpNC1piNominalValue{
      fhicl::Name("RvbarpNC1piNominal"),
      fhicl::Comment(
          "RvbarpNC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarpNC1piTweakDefinition{
      fhicl::Name("RvbarpNC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarpNC1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarpNC2piNominalValue{
      fhicl::Name("RvbarpNC2piNominal"),
      fhicl::Comment(
          "RvbarpNC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarpNC2piTweakDefinition{
      fhicl::Name("RvbarpNC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarpNC2pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarnCC1piNominalValue{
      fhicl::Name("RvbarnCC1piNominal"),
      fhicl::Comment(
          "RvbarnCC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarnCC1piTweakDefinition{
      fhicl::Name("RvbarnCC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarnCC1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarnCC2piNominalValue{
      fhicl::Name("RvbarnCC2piNominal"),
      fhicl::Comment(
          "RvbarnCC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarnCC2piTweakDefinition{
      fhicl::Name("RvbarnCC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarnCC2pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarnNC1piNominalValue{
      fhicl::Name("RvbarnNC1piNominal"),
      fhicl::Comment(
          "RvbarnNC1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarnNC1piTweakDefinition{
      fhicl::Name("RvbarnNC1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarnNC1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RvbarnNC2piNominalValue{
      fhicl::Name("RvbarnNC2piNominal"),
      fhicl::Comment(
          "RvbarnNC2pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> RvbarnNC2piTweakDefinition{
      fhicl::Name("RvbarnNC2piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RvbarnNC2pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  //********* END RES options ***********

  //********* START COH options ***********
  fhicl::Atom<double> MaCOHpiNominalValue{
      fhicl::Name("MaCOHpiNominal"),
      fhicl::Comment(
          "MaCOHpi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> MaCOHpiTweakDefinition{
      fhicl::Name("MaCOHpiTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MaCOHpi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> R0COHpiNominalValue{
      fhicl::Name("R0COHpiNominal"),
      fhicl::Comment(
          "R0COHpi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> R0COHpiTweakDefinition{
      fhicl::Name("R0COHpiTweakDefinition"),
      fhicl::Comment(
          "Defines variations in R0COHpi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};
  //********* END COH options ***********

  //********* START DIS options ***********
  fhicl::Atom<bool> DISIsShapeOnly{
      fhicl::Name("DISIsShapeOnly"),
      fhicl::Comment(
          "Whether variations in Bodek-Yang DIS parameters only affect the "
          "shape of the event rate distribution and not the total "
          "normalization. Default == false"),
      false};

  fhicl::Atom<double> NormDISCCNominalValue{
      fhicl::Name("NormDISCCNominal"),
      fhicl::Comment(
          "NormDISCC Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> NormDISCCTweakDefinition{
      fhicl::Name("NormDISCCTweakDefinition"),
      fhicl::Comment(
          "Defines variations in NormDISCC. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> AhtBYNominalValue{
      fhicl::Name("AhtBYNominal"),
      fhicl::Comment(
          "AhtBY Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> AhtBYTweakDefinition{
      fhicl::Name("AhtBYTweakDefinition"),
      fhicl::Comment(
          "Defines variations in AhtBY. If not specified, but NominalValue is, "
          "then a single response at the new nominal value will be generated. "
          "If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> BhtBYNominalValue{
      fhicl::Name("BhtBYNominal"),
      fhicl::Comment(
          "BhtBY Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> BhtBYTweakDefinition{
      fhicl::Name("BhtBYTweakDefinition"),
      fhicl::Comment(
          "Defines variations in BhtBY. If not specified, but NominalValue is, "
          "then a single response at the new nominal value will be generated. "
          "If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> CV1uBYNominalValue{
      fhicl::Name("CV1uBYNominal"),
      fhicl::Comment(
          "CV1uBY Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> CV1uBYTweakDefinition{
      fhicl::Name("CV1uBYTweakDefinition"),
      fhicl::Comment(
          "Defines variations in CV1uBY. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> CV2uBYNominalValue{
      fhicl::Name("CV2uBYNominal"),
      fhicl::Comment(
          "CV2uBY Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine."),
      0xdeadb33f};
  fhicl::Atom<std::string> CV2uBYTweakDefinition{
      fhicl::Name("CV2uBYTweakDefinition"),
      fhicl::Comment(
          "Defines variations in CV2uBY. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> RnubarnuCCNominalValue{
      fhicl::Name("RnubarnuCCNominal"),
      fhicl::Comment(
          "RnubarnuCC Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> RnubarnuCCTweakDefinition{
      fhicl::Name("RnubarnuCCTweakDefinition"),
      fhicl::Comment(
          "Defines variations in RnubarnuCC. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> DISNuclModNominalValue{
      fhicl::Name("DISNuclModNominal"),
      fhicl::Comment("DISNuclMod Nominal value, set to 0 to use in its nominal "
                     "(GENIE-defined) state. It appears that this may not be "
                     "currently working."),
      0xdeadb33f};

  fhicl::Atom<double> AGKY_xF1piNominalValue{
      fhicl::Name("AGKY_xF1piNominal"),
      fhicl::Comment(
          "AGKY_xF1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> AGKY_xF1piTweakDefinition{
      fhicl::Name("AGKY_xF1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in AGKY_xF1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> AGKY_pT1piNominalValue{
      fhicl::Name("AGKY_pT1piNominal"),
      fhicl::Comment(
          "AGKY_pT1pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> AGKY_pT1piTweakDefinition{
      fhicl::Name("AGKY_pT1piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in AGKY_pT1pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  //********* END DIS options ***********

  //********* START FSI options ***********
  fhicl::Atom<double> FormZoneNominalValue{
      fhicl::Name("FormZoneNominal"),
      fhicl::Comment(
          "FormZone Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FormZoneTweakDefinition{
      fhicl::Name("FormZoneTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FormZone. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> MFP_piNominalValue{
      fhicl::Name("MFP_piNominal"),
      fhicl::Comment(
          "MFP_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> MFP_piTweakDefinition{
      fhicl::Name("MFP_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MFP_pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrCEx_piNominalValue{
      fhicl::Name("FrCEx_piNominal"),
      fhicl::Comment(
          "FrCEx_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrCEx_piTweakDefinition{
      fhicl::Name("FrCEx_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrCEx_pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrElas_piNominalValue{
      fhicl::Name("FrElas_piNominal"),
      fhicl::Comment(
          "FrElas_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrElas_piTweakDefinition{
      fhicl::Name("FrElas_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrElas_pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrInel_piNominalValue{
      fhicl::Name("FrInel_piNominal"),
      fhicl::Comment(
          "FrInel_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrInel_piTweakDefinition{
      fhicl::Name("FrInel_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrInel_pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrAbs_piNominalValue{
      fhicl::Name("FrAbs_piNominal"),
      fhicl::Comment(
          "FrAbs_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrAbs_piTweakDefinition{
      fhicl::Name("FrAbs_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrAbs_pi. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrPiProd_piNominalValue{
      fhicl::Name("FrPiProd_piNominal"),
      fhicl::Comment(
          "FrPiProd_pi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrPiProd_piTweakDefinition{
      fhicl::Name("FrPiProd_piTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrPiProd_pi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> MFP_NNominalValue{
      fhicl::Name("MFP_NNominal"),
      fhicl::Comment(
          "MFP_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> MFP_NTweakDefinition{
      fhicl::Name("MFP_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in MFP_N. If not specified, but NominalValue is, "
          "then a single response at the new nominal value will be generated. "
          "If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrCEx_NNominalValue{
      fhicl::Name("FrCEx_NNominal"),
      fhicl::Comment(
          "FrCEx_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrCEx_NTweakDefinition{
      fhicl::Name("FrCEx_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrCEx_N. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrElas_NNominalValue{
      fhicl::Name("FrElas_NNominal"),
      fhicl::Comment(
          "FrElas_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrElas_NTweakDefinition{
      fhicl::Name("FrElas_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrElas_N. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrInel_NNominalValue{
      fhicl::Name("FrInel_NNominal"),
      fhicl::Comment(
          "FrInel_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrInel_NTweakDefinition{
      fhicl::Name("FrInel_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrInel_N. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrAbs_NNominalValue{
      fhicl::Name("FrAbs_NNominal"),
      fhicl::Comment(
          "FrAbs_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrAbs_NTweakDefinition{
      fhicl::Name("FrAbs_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrAbs_N. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> FrPiProd_NNominalValue{
      fhicl::Name("FrPiProd_NNominal"),
      fhicl::Comment(
          "FrPiProd_N Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> FrPiProd_NTweakDefinition{
      fhicl::Name("FrPiProd_NTweakDefinition"),
      fhicl::Comment(
          "Defines variations in FrPiProd_N. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};
  //********* END FSI options ***********

  //********* START Other options ***********

  fhicl::Atom<double> CCQEPauliSupViaKFNominalValue{
      fhicl::Name("CCQEPauliSupViaKFNominal"),
      fhicl::Comment(
          "CCQEPauliSupViaKF Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> CCQEPauliSupViaKFTweakDefinition{
      fhicl::Name("CCQEPauliSupViaKFTweakDefinition"),
      fhicl::Comment(
          "Defines variations in CCQEPauliSupViaKF. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<bool> CCQEMomDistroFGtoSFNominalValue{
      fhicl::Name("CCQEMomDistroFGtoSFNominal"),
      fhicl::Comment("Whether to apply CCQEMomDistroFGtoSF reweighting."),
      false};

  fhicl::Atom<double> BR1gammaNominalValue{
      fhicl::Name("BR1gammaNominal"),
      fhicl::Comment(
          "BR1gamma Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> BR1gammaTweakDefinition{
      fhicl::Name("BR1gammaTweakDefinition"),
      fhicl::Comment(
          "Defines variations in BR1gamma. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> BR1etaNominalValue{
      fhicl::Name("BR1etaNominal"),
      fhicl::Comment(
          "BR1eta Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> BR1etaTweakDefinition{
      fhicl::Name("BR1etaTweakDefinition"),
      fhicl::Comment(
          "Defines variations in BR1eta. If not specified, but NominalValue "
          "is, then a single response at the new nominal value will be "
          "generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};

  fhicl::Atom<double> Theta_Delta2NpiNominalValue{
      fhicl::Name("Theta_Delta2NpiNominal"),
      fhicl::Comment(
          "Theta_Delta2Npi Nominal value, set to 0 to use in its nominal "
          "(GENIE-defined) state, If unitsAreNatural == true, then this "
          "value will always be passed to the weight engine. It appears "
          "that this may not be currently working."),
      0xdeadb33f};
  fhicl::Atom<std::string> Theta_Delta2NpiTweakDefinition{
      fhicl::Name("Theta_Delta2NpiTweakDefinition"),
      fhicl::Comment(
          "Defines variations in Theta_Delta2Npi. If not specified, but "
          "NominalValue is, then a single response at the new nominal value "
          "will be generated. If specified as \"{sigma_low, "
          "sigma_up}\", then normally distributed throws will be "
          "made with those as the 2*half gaussian widths above and below the "
          "central value. If the list is specified as \"(-1,2,3,...)\", then "
          "spline knots will be made at specified values, respects "
          "unitsAreNatural. If specified as \"[1,2,3]\", then discrete "
          "non-splineable tweaks will be made at those values, respects "
          "unitsAreNatural."),
      ""};
  //********* END Other options ***********
};

void MakeThrowsIfNeeded(larsyst::SystParamHeader &sph,
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

larsyst::SystParamHeader
BuildHeaderFromNomAndTweakDefintion(double nominal = 0xdeadb33f,
                                    std::string tweakDefinition = "") {
  larsyst::SystParamHeader sph;
  sph.centralParamValue = (nominal == 0xdeadb33f) ? 0 : nominal;
  larsyst::trim(tweakDefinition);

  if (tweakDefinition.size()) {
    char fchar = tweakDefinition.front();
    tweakDefinition = tweakDefinition.substr(1, tweakDefinition.length() - 2);
    larsyst::trim(tweakDefinition);
    if (fchar == '(') { // Spline knots

      sph.paramVariations = larsyst::BuildDoubleList(tweakDefinition);
      sph.isSplineable = true;
    } else if (fchar == '[') { // Discrete tweaks

      sph.paramVariations = larsyst::ParseToVect<double>(tweakDefinition, ",");
    } else if (fchar == '{') { // OneSigmaShifts

      std::vector<double> sigShifts =
          larsyst::ParseToVect<double>(tweakDefinition, ",");
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

size_t GetParamIndex(larsyst::SystMetaData const &md, std::string const &name) {
  for (size_t it = 0; it < md.headers.size(); ++it) {
    if (md.headers[it].prettyName == name) {
      return it;
    }
  }
  return larsyst::kParamUnhandled<size_t>;
}
bool HasParam(larsyst::SystMetaData const &md, std::string const &name) {
  return (GetParamIndex(md, name) != larsyst::kParamUnhandled<size_t>);
}

larsyst::SystParamHeader const &GetParam(larsyst::SystMetaData const &md,
                                         std::string const &name) {
  size_t idx = GetParamIndex(md, name);
  if (idx != larsyst::kParamUnhandled<size_t>) {
    return md.headers[idx];
  }
  std::cout << "[ERROR]: Tried to get parameter named " << std::quoted(name)
            << " from a SystMetaData instance, but it doesn't exist."
            << std::endl;
  throw;
}
larsyst::SystParamHeader &GetParam(larsyst::SystMetaData &md,
                                   std::string const &name) {
  for (auto &sph : md.headers) {
    if (sph.prettyName == name) {
      return sph;
    }
  }
  throw;
}

larsyst::SystMetaData
ConfigureQEParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg,
                            larsyst::paramId_t firstParamId) {
  larsyst::SystMetaData QEmd;

  // Axial FFs
  bool DipoleNormCCQEIsUsed = (cfg().NormCCQENominalValue() != 0xdeadb33f) ||
                              cfg().NormCCQETweakDefinition().size();
  bool DipoleIsShapeOnly = cfg().MAQEIsShapeOnly() || DipoleNormCCQEIsUsed;
  bool DipoleMaIsUsed = (cfg().MaCCQENominalValue() != 0xdeadb33f) ||
                        cfg().MaCCQETweakDefinition().size();

  bool IsDipoleReWeight =
      DipoleIsShapeOnly || DipoleNormCCQEIsUsed || DipoleMaIsUsed;

  bool ZNormIsUsed = (cfg().ZNormCCQENominalValue() != 0xdeadb33f) ||
                     cfg().ZNormCCQETweakDefinition().size();
  bool ZExpA1IsUsed = (cfg().ZExpA1CCQENominalValue() != 0xdeadb33f) ||
                      cfg().ZExpA1CCQETweakDefinition().size();
  bool ZExpA2IsUsed = (cfg().ZExpA2CCQENominalValue() != 0xdeadb33f) ||
                      cfg().ZExpA2CCQETweakDefinition().size();
  bool ZExpA3IsUsed = (cfg().ZExpA3CCQENominalValue() != 0xdeadb33f) ||
                      cfg().ZExpA3CCQETweakDefinition().size();
  bool ZExpA4IsUsed = (cfg().ZExpA4CCQENominalValue() != 0xdeadb33f) ||
                      cfg().ZExpA4CCQETweakDefinition().size();

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
      larsyst::SystParamHeader qenorm = BuildHeaderFromNomAndTweakDefintion(
          cfg().NormCCQENominalValue(), cfg().NormCCQETweakDefinition());
      qenorm.prettyName = "NormCCQE";
      qenorm.unitsAreNatural = false;
      qenorm.systParamId = firstParamId++;
      QEmd.headers.push_back(std::move(qenorm));
    }
    if (DipoleMaIsUsed) {
      larsyst::SystParamHeader maqe = BuildHeaderFromNomAndTweakDefintion(
          cfg().MaCCQENominalValue(), cfg().MaCCQETweakDefinition());
      maqe.prettyName = "MaCCQE";
      maqe.unitsAreNatural = false;
      maqe.systParamId = firstParamId++;
      if (DipoleIsShapeOnly) {
        maqe.opts.push_back("shape");
      }
      QEmd.headers.push_back(std::move(maqe));
    }
  } else if (IsZExpReWeight) {
    larsyst::SystParamHeader zexpaxFF;
    zexpaxFF.prettyName = "ZExpAxFF";
    zexpaxFF.systParamId = firstParamId++;
    zexpaxFF.isCorrection = true;
    zexpaxFF.centralParamValue = 1;
    QEmd.headers.push_back(std::move(zexpaxFF));

    if (ZNormIsUsed) {
      larsyst::SystParamHeader znorm = BuildHeaderFromNomAndTweakDefintion(
          cfg().ZNormCCQENominalValue(), cfg().ZNormCCQETweakDefinition());
      znorm.prettyName = "ZNormCCQE";
      znorm.unitsAreNatural = false;
      znorm.systParamId = firstParamId++;
      if (znorm.isSplineable) {
        std::cout << "[ERROR]: Attempted to build spline from Z-Expansion "
                     "ZNorm parameter. Multi-dimensional splines are not yet "
                     "supported, this parameter must currently be used as a "
                     "multi-sim or with hand-picked offsets."
                  << std::endl;
        throw;
      }
      QEmd.headers.push_back(std::move(znorm));
    }
    if (ZExpA1IsUsed || ZExpA2IsUsed || ZExpA3IsUsed || ZExpA4IsUsed) {
      larsyst::SystParamHeader ZExp;
      ZExp.prettyName = "ZExpAVariationResponse";
      ZExp.systParamId = firstParamId++;
      QEmd.headers.push_back(std::move(ZExp));

      if (ZExpA1IsUsed) {
        larsyst::SystParamHeader ZA1 = BuildHeaderFromNomAndTweakDefintion(
            cfg().ZExpA1CCQENominalValue(), cfg().ZExpA1CCQETweakDefinition());
        ZA1.prettyName = "ZExpA1CCQE";
        ZA1.unitsAreNatural = false;

        ZA1.systParamId = firstParamId++;
        if (ZA1.isSplineable) {
          std::cout << "[ERROR]: Attempted to build spline from Z-Expansion "
                       "ZExpA1 parameter. Multi-dimensional splines are not "
                       "yet supported, this parameter must currently be used "
                       "as a multi-sim or with hand-picked offsets."
                    << std::endl;
          throw;
        }
        ZA1.isResponselessParam = true;
        ZA1.responseParamId = ZExp.systParamId;
        QEmd.headers.push_back(std::move(ZA1));
      }
      if (ZExpA2IsUsed) {
        larsyst::SystParamHeader ZA2 = BuildHeaderFromNomAndTweakDefintion(
            cfg().ZExpA2CCQENominalValue(), cfg().ZExpA2CCQETweakDefinition());
        ZA2.prettyName = "ZExpA2CCQE";
        ZA2.unitsAreNatural = false;

        ZA2.systParamId = firstParamId++;
        if (ZA2.isSplineable) {
          std::cout << "[ERROR]: Attempted to build spline from Z-Expansion "
                       "ZExpA2 parameter. Multi-dimensional splines are not "
                       "yet supported, this parameter must currently be used "
                       "as a multi-sim or with hand-picked offsets."
                    << std::endl;
          throw;
        }
        ZA2.isResponselessParam = true;
        ZA2.responseParamId = ZExp.systParamId;
        QEmd.headers.push_back(std::move(ZA2));
      }
      if (ZExpA3IsUsed) {
        larsyst::SystParamHeader ZA3 = BuildHeaderFromNomAndTweakDefintion(
            cfg().ZExpA3CCQENominalValue(), cfg().ZExpA3CCQETweakDefinition());
        ZA3.prettyName = "ZExpA3CCQE";
        ZA3.unitsAreNatural = false;

        ZA3.systParamId = firstParamId++;
        if (ZA3.isSplineable) {
          std::cout << "[ERROR]: Attempted to build spline from Z-Expansion "
                       "ZExpA3 parameter. Multi-dimensional splines are not "
                       "yet supported, this parameter must currently be used "
                       "as a multi-sim or with hand-picked offsets."
                    << std::endl;
          throw;
        }
        ZA3.isResponselessParam = true;
        ZA3.responseParamId = ZExp.systParamId;
        QEmd.headers.push_back(std::move(ZA3));
      }
      if (ZExpA4IsUsed) {
        larsyst::SystParamHeader ZA4 = BuildHeaderFromNomAndTweakDefintion(
            cfg().ZExpA4CCQENominalValue(), cfg().ZExpA4CCQETweakDefinition());
        ZA4.prettyName = "ZExpA4CCQE";
        ZA4.unitsAreNatural = false;

        ZA4.systParamId = firstParamId++;
        if (ZA4.isSplineable) {
          std::cout << "[ERROR]: Attempted to build spline from Z-Expansion "
                       "ZExpA4 parameter. Multi-dimensional splines are not "
                       "yet supported, this parameter must currently be used "
                       "as a multi-sim or with hand-picked offsets."
                    << std::endl;
          throw;
        }
        ZA4.isResponselessParam = true;
        ZA4.responseParamId = ZExp.systParamId;
        QEmd.headers.push_back(std::move(ZA4));
      }
    }
  }

  bool DipoleVecFFShapeIsUsed = !cfg().VecFFCCQEIsBBA();
  if (DipoleVecFFShapeIsUsed) {
    larsyst::SystParamHeader vecFFQE;
    vecFFQE.prettyName = "dipoleVecFFCCQE";
    vecFFQE.systParamId = firstParamId++;
    vecFFQE.isCorrection = true;
    vecFFQE.centralParamValue = 1;
    QEmd.headers.push_back(std::move(vecFFQE));
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
    for (auto &name :
         {"ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE"}) {
      if (!HasParam(QEmd, name)) {
        continue;
      }
      if (!NVariations) {
        NVariations = GetParam(QEmd, name).paramVariations.size();
        continue;
      }
      if (NVariations != GetParam(QEmd, name).paramVariations.size()) {
        std::cout
            << "[ERROR]: ZExpansion configuration error. Each form factor "
               "parameter specified must have the same number of variations, "
            << name << " specifies "
            << GetParam(QEmd, name).paramVariations.size()
            << ", but a previous parameter specified " << NVariations << "."
            << std::endl;
        throw;
      }
    }
    std::vector<double> dummyParamVars;
    for (size_t i = 0; i < NVariations; ++i) {
      dummyParamVars.push_back(i);
    }

    GetParam(QEmd, "ZExpAVariationResponse").paramVariations =
        std::move(dummyParamVars);
  }

  return QEmd;
}

larsyst::SystMetaData
ConfigureRESParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg) {
  larsyst::SystMetaData RESmd;
  return RESmd;
}

larsyst::SystMetaData
ConfigureCOHParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg) {
  larsyst::SystMetaData COHmd;
  return COHmd;
}

larsyst::SystMetaData
ConfigureDISParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg) {
  larsyst::SystMetaData DISmd;
  return DISmd;
}

larsyst::SystMetaData
ConfigureFSIParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg) {
  larsyst::SystMetaData FSImd;
  return FSImd;
}

larsyst::SystMetaData
ConfigureOtherParameterHeaders(fhicl::Table<GENIEReWeightConfig> const &cfg) {
  larsyst::SystMetaData Othermd;
  return Othermd;
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> ConfigureQEWeightEngine(
    larsyst::SystMetaData const &QEmd,
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
      bool maqeHasShapeOpt =
          (HasParam(QEmd, "MaCCQE") &&
           (std::find(QEmd.headers[maqeidx].opts.begin(),
                      QEmd.headers[maqeidx].opts.end(),
                      "all") != QEmd.headers[maqeidx].opts.end()));
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

      size_t zaxidx = GetParamIndex(QEmd, "ZExpAxFF");
      GReWeightEngine->AdoptWghtCalc("xsec_ccqe_ZExp",
                                     new genie::rew::GReWeightNuXSecCCQEaxial);
      param_map[zaxidx].insert(
          {genie::rew::kXSecTwkDial_AxFFCCQEshape, zaxidx});

      if (HasParam(QEmd, "ZNormCCQE")) {
        size_t znormidx = GetParamIndex(QEmd, "ZNormCCQE");
        param_map[znormidx].insert(
            {genie::rew::kXSecTwkDial_ZNormCCQE, znormidx});
      }

      // Turns of each zexpansion dial effect a response via a single parameter.
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

  if (HasParam(QEmd, "dipoleVecFFCCQE")) {
    GReWeightEngine->AdoptWghtCalc("xsec_ccqe_vecFF",
                                   new genie::rew::GReWeightNuXSecCCQEvec);
    size_t vecffidx = GetParamIndex(QEmd, "dipoleVecFFCCQE");
    param_map[vecffidx].insert(
        {genie::rew::kXSecTwkDial_VecFFCCQEshape, vecffidx});
  }
  return param_map;
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureRESWeightEngine(
    larsyst::SystMetaData const &RESmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureCOHWeightEngine(
    larsyst::SystMetaData const &COHmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureDISWeightEngine(
    larsyst::SystMetaData const &DISmd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureFSIWeightEngine(
    larsyst::SystMetaData const &FSImd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureOtherWeightEngine(
    larsyst::SystMetaData const &Othermd,
    std::unique_ptr<genie::rew::GReWeight> &GReWeightEngine) {
  return {};
}

#endif
