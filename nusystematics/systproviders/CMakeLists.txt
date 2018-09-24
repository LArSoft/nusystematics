
SET(GENIE_EV_LIBS
  ${GENIE_LIB_LIST}
  ${LOG4CPP}
  ${XML2}
  ${LHAPDF}
  ${ROOT_BASIC_LIB_LIST}
  ${ROOT_EG}
  ${ROOT_TREEPLAYER}
  ${ROOT_FFTW}
  ${ROOT_REFLEX}
  ${ROOTSYS}/lib/libEGPythia6.so
  ${ROOT_GUI}
  ${PYTHIA6}
)

art_make_library(LIBRARY_NAME nusystematics_GENIEReWeight_config
  SOURCE
    GENIEReWeightParamConfig.cc
    GENIEReWeightEngineConfig.cc
  LIBRARIES
    systematicstools_utility
    systematicstools_interface
    nutools_EventGeneratorBase_GENIE
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib cetlib_except
    ${CLHEP}
    ${GENIE_EV_LIBS}
  )

simple_plugin(GENIEReWeight tool
    nusystematics_GENIEReWeight_config
    systematicstools_utility
    systematicstools_interface
    nutools_EventGeneratorBase_GENIE
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib cetlib_except
    ${CLHEP}
    ${GENIE_EV_LIBS}
 )

simple_plugin(MKSinglePiTemplate tool
  systematicstools_utility
  systematicstools_interface
  nutools_EventGeneratorBase_GENIE
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib cetlib_except
    ${CLHEP}
    ${GENIE_EV_LIBS}
 )

simple_plugin(MINERvAq0q3Weighting tool
  systematicstools_utility
  systematicstools_interface
  nutools_EventGeneratorBase_GENIE
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib cetlib_except
    ${CLHEP}
    ${GENIE_EV_LIBS}
 )
