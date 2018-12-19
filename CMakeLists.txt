# ======================================================================
#  nusystematics main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/nusystematics/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.3 or later
cmake_minimum_required (VERSION 3.3)

project(nusystematics)

if(NOT DEFINED USEART) #by default build in ART mode
  SET(USEART TRUE)
endif()

if(USEART) # This is the standard, CET-empowered build for building the art module and associated tools

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
                      $ENV{CETBUILDTOOLS_DIR}/Modules
                      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic )

cet_report_compiler_flags()

find_ups_product( systematicstools )
find_ups_product( nutools )
find_ups_product( art )
find_ups_product( cetbuildtools )
find_ups_product( genie )
find_ups_product( log4cpp )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# include search path
include_directories ( $ENV{GENIE_INC}/GENIE )

# get libxml2
cet_find_library( XML2 NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib NO_DEFAULT_PATH )

# BEGIN GENIE VERSION CHECK
# Determine the major version number of GENIE based on the GENIE_VERSION
# environment variable. Based on this, define the GENIE_PRE_R3 preprocessor
# macro (needed to set the right GENIE header locations, etc.) if needed.
#
# The top-level nutools CMakeLists.txt does this as well, but I don't think that
# we can access the result of that version check in the nusystematics build.

parse_ups_version( ${GENIE_VERSION} )

# the GENIE libraries are already defined, but we need the lib list
if( ${VMAJ} LESS 3 )

  add_definitions(-DGENIE_PRE_R3)

  # GENIE v2 library list
  set(GENIE_LIB_LIST  ${GALGORITHM}
                    ${GBARYONRESONANCE}
                    ${GBASE}
                    ${GBODEKYANG}
                    ${GCHARM}
                    ${GCOH}
                    ${GDFRC}
                    ${GDIS}
                    ${GCROSSSECTIONS}
                    ${GDECAY}
                    ${GELAS}
                    ${GELFF}
                    ${GHEP}
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GEVGDRIVERS}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GINTERACTION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GMESSENGER}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNTUPLE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GUTILS}
                    ${GGEO}
                    ${GFLUXDRIVERS}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )

  #MESSAGE("--NUSYSTEMATICS-- GENIE_VERSION=${VMAJ} ${VMIN} -DGENIE_PRE_R3")
else()

  # GENIE v3 library list
  set (GENIE_LIB_LIST ${GFWMSG}
                      ${GFWREG}
                      ${GFWALG}
                      ${GFWINT}
                      ${GFWGHEP}
                      ${GFWNUM}
                      ${GSL}  # FWNUM relies on GSL
                      ${GFWUTL}
                      ${GFWPARDAT}
                      ${GFWEG}
                      ${GFWNTP}
                      ${GPHXSIG}
                      ${GPHPDF}
                      ${GPHNUCLST}
                      ${GPHCMN}
                      ${GPHDCY}
                      ${GPHHADTRANSP}
                      ${GPHHADNZ}
                      ${GPHDEEX}
                      ${GPHAMNGXS}
                      ${GPHAMNGEG}
                      ${GPHCHMXS}
                      ${GPHCOHXS}
                      ${GPHCOHEG}
                      ${GPHDISXS}
                      ${GPHDISEG}
                      ${GPHDFRCXS}
                      ${GPHDFRCEG}
                      ${GPHGLWRESXS}
                      ${GPHGLWRESEG}
                      ${GPHIBDXS}
                      ${GPHIBDEG}
                      ${GPHMNUCXS}
                      ${GPHMNUCEG}
                      ${GPHMEL}
                      ${GPHNUELXS}
                      ${GPHNUELEG}
                      ${GPHQELXS}
                      ${GPHQELEG}
                      ${GPHRESXS}
                      ${GPHRESEG}
                      ${GPHSTRXS}
                      ${GPHSTREG}
                      ${GPHNDCY}
                      ${GTLGEO}
                      ${GTLFLX}
                      ${GRWFWK}
                      ${GRWIO}
                      ${GRWCLC})

  #MESSAGE("--NUSYSTEMATICS-- GENIE_VERSION=${VMAJ} ${VMIN} not defining -DGENIE_PRE_R3")
endif()

#MESSAGE("--NUSYSTEMATICS-- GENIE ${VMAJ} GENIE_LIB_LIST=${GENIE_LIB_LIST}")
# END GENIE VERSION CHECK

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(nusystematics)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

else() #

  #Changes default install path to be a subdirectory of the build dir.
  #Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
  if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
    "/usr/local")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
  elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
  endif()
  message(STATUS "Installing to: ${CMAKE_INSTALL_PREFIX}")

  add_subdirectory(nusystematics/artless)

endif()


