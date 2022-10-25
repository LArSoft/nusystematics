
SET(GENIE_EV_LIBS
  ${GENIE_LIB_LIST}
  log4cpp::log4cpp
  LibXml2::LibXml2
  LHAPDF::LHAPDF
  ROOT::Core
  ROOT::EG
  ROOT::TreePlayer
  ROOT::FFTW
  ROOT::EGPythia6
  ROOT::Gui
  Pythia6::Pythia6
)

art_make_library(LIBRARY_NAME nusystematics_GENIEReWeight_config
  SOURCE
    GENIEReWeightParamConfig.cc
    GENIEReWeightEngineConfig.cc
  LIBRARIES
    systematicstools::utility
    systematicstools::interface
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

cet_build_plugin(GENIEReWeight art::tool
  LIBRARIES PRIVATE
    nusystematics_GENIEReWeight_config
    systematicstools::utility
    systematicstools::interface
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
  cet_build_plugin(${tool} art::tool
    LIBRARIES PRIVATE
      systematicstools::utility
      systematicstools::interface
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
