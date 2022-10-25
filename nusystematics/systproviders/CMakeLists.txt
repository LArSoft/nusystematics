
cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
SET(GENIE_EV_LIBS
  ${GENIE_LIB_LIST}
  log4cpp::log4cpp
  XML2
  LHAPDF
  ROOT::Core
  ROOT::EG
  ROOT::TreePlayer
  ROOT::FFTW
  ROOT::EGPythia6
  ROOT::Gui
  PYTHIA6
)

art_make_library(LIBRARY_NAME nusystematics_GENIEReWeight_config
  SOURCE
    GENIEReWeightParamConfig.cc
    GENIEReWeightEngineConfig.cc
  LIBRARIES
    systematicstools_utility
    systematicstools_interface
    nugen::EventGeneratorBase_GENIE
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    ${GENIE_EV_LIBS}
    ifdh_art::IFDH_service
  )

simple_plugin(GENIEReWeight tool
    nusystematics_GENIEReWeight_config
    systematicstools_utility
    systematicstools_interface
    nugen::EventGeneratorBase_GENIE
    nusimdata::SimulationBase
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Provenance
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    ${GENIE_EV_LIBS}
 )

SET(TOOLS MKSinglePiTemplate BeRPAWeight MiscInteractionSysts NOvAStyleNonResPionNorm MINERvAE2p2h EbLepMomShift MINERvAq0q3Weighting FSILikeEAvailSmearing)
foreach(tool ${TOOLS})
  simple_plugin(${tool} tool
    systematicstools_utility
    systematicstools_interface
    nugen::EventGeneratorBase_GENIE
    nusimdata::SimulationBase
      art::Framework_Core
      art::Framework_Principal
      art::Framework_Services_Registry
      art::Framework_Services_Optional_RandomNumberGenerator_service
      art::Persistency_Provenance
      canvas::canvas
      messagefacility::MF_MessageLogger
      fhiclcpp::fhiclcpp
      cetlib::cetlib
      cetlib_except::cetlib_except
      CLHEP::CLHEP
      ${GENIE_EV_LIBS}
   )
endforeach()
