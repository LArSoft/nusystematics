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
cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(nusystematics VERSION 01.00.15 LANGUAGES CXX)

if(NOT DEFINED USEART) #by default build in ART mode
  SET(USEART TRUE)
endif()

if(USEART) # This is the standard, CET-empowered build for building the art module and associated tools


# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic )

cet_report_compiler_flags()

unset(GSL_FOUND CACHE) # Workaround infelicity in CMake's FindPkgConfig.cmake.
find_ups_product(gsl)
find_ups_product( systematicstools )
find_ups_product( nusimdata )
find_ups_product( nugen )
find_ups_product( art )
find_ups_product(art_root_io)
#find_ups_product( genie )
find_ups_product( log4cpp )
find_ups_product( ifdh_art )
find_ups_product( lhapdf )
find_package(GENIE REQUIRED)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# include search path
include_directories ( $ENV{GENIE_INC}/GENIE )

cet_find_library( PYTHIA6 NAMES Pythia6 PATHS ENV PYLIB NO_DEFAULT_PATH )
cet_find_library( XML2 NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib NO_DEFAULT_PATH )
cet_find_library( LHAPDF NAMES LHAPDF PATHS ENV LHAPDF_LIB NO_DEFAULT_PATH )
cet_find_library( LOG4CPP NAMES log4cpp PATHS ENV LOG4CPP_LIB NO_DEFAULT_PATH )


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
