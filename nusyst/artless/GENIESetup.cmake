# Adapted for nusyst by L. Pickering, originally from the NUISANCE package

# Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

################################################################################
#    This file is part of NUISANCE.
#
#    NUISANCE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NUISANCE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

#################################  GENIE  ######################################
if(GENIE STREQUAL "")
  message(FATAL_ERROR "Variable GENIE is not defined. "
    "The location of a pre-built GENIE install must be defined either as"
    " $ cmake -DGENIE=/path/to/GENIE or as and environment vairable"
    " $ export GENIE=/path/to/GENIE")
endif()

find_program(GENIECFG genie-config)
if(GENIECFG STREQUAL "GENIECFG-NOTFOUND")
  message(FATAL_ERROR "Failed to find the genie-config program. Is the GENIE environment set up?")
endif()

execute_process (COMMAND ${GENIECFG}
  --libs OUTPUT_VARIABLE GENIE_LD_FLAGS_STR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND ${GENIECFG}
  --topsrcdir OUTPUT_VARIABLE GENIE_INCLUDES_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REGEX MATCH "-L\([^ ]+\) \(.*\)$" PARSE_GENIE_LIBS_MATCH ${GENIE_LD_FLAGS_STR})

message(STATUS "genie-config --libs: ${GENIE_LD_FLAGS_STR}")

if(NOT PARSE_GENIE_LIBS_MATCH)
  message(FATAL_ERROR "Expected to be able to parse the result of genie-config --libs to a lib directory and a list of libraries to include, but got: \"${GENIE_LD_FLAGS_STR}\"")
endif()

set(GENIE_LIB_DIR ${CMAKE_MATCH_1})
set(GENIE_LIBS_RAW ${CMAKE_MATCH_2})
string(REPLACE "-l" "" GENIE_LIBS_STRIPED "${GENIE_LIBS_RAW}")

message(STATUS "GENIE version : ${GENIE_VERSION}")
message(STATUS "GENIE libdir  : ${GENIE_LIB_DIR}")
message(STATUS "GENIE libs    : ${GENIE_LIBS_STRIPED}")

string(REGEX MATCH "ReWeight" WASMATCHED ${GENIE_LIBS_STRIPED})
if(NOT WASMATCHED)
  set(GENIE_LIBS_STRIPED "GReWeight ${GENIE_LIBS_STRIPED}")
  message(STATUS "Force added ReWeight library: ${GENIE_LIBS_STRIPED}")
endif()

string(REPLACE " " ";" GENIE_LIBS_LIST "${GENIE_LIBS_STRIPED}")
message(STATUS "genie-config --libs -- MATCH1: ${CMAKE_MATCH_1}")
message(STATUS "genie-config --libs -- MATCH2: ${CMAKE_MATCH_2}")
message(STATUS "genie-config --libs -- libs stripped: ${GENIE_LIBS_STRIPED}")
message(STATUS "genie-config --libs -- libs list: ${GENIE_LIBS_LIST}")

##############################  VARIABLES  #####################################
function(CheckAndSetDefaultEnv VARNAME DEFAULT CACHETYPE DOCSTRING ENVNAME)
  if(NOT DEFINED ${VARNAME})
    if(DEFINED ENV{${ENVNAME}} AND NOT $ENV{${ENVNAME}} STREQUAL "")
      set(${VARNAME} $ENV{${ENVNAME}} CACHE ${CACHETYPE} ${DOCSTRING})
    else()
      set(${VARNAME} ${DEFAULT} CACHE ${CACHETYPE} ${DOCSTRING})
    endif()
  else()
    set(${VARNAME} ${${VARNAME}} CACHE ${CACHETYPE} ${DOCSTRING})
    unset(${VARNAME})
  endif()
endfunction()

CheckAndSetDefaultEnv(GENIE "" PATH "Path to GENIE source tree root directory. Overrides environment variable \$GENIE <>" GENIE)
CheckAndSetDefaultEnv(LHAPDF_LIB "" PATH "Path to pre-built LHAPDF libraries. Overrides environment variable \$LHAPDF_LIB. <>" LHAPDF_LIB)
CheckAndSetDefaultEnv(LHAPDF_INC "" PATH "Path to installed LHAPDF headers. Overrides environment variable \$LHAPDF_INC. <>" LHAPDF_INC)
CheckAndSetDefaultEnv(LHAPATH "" PATH "Path to LHA PDF inputs. Overrides environment variable \$LHAPATH. <>" LHAPATH)
CheckAndSetDefaultEnv(LIBXML2_LIB "" PATH "Path to pre-built LIBXML2 libraries. Overrides environment variable \$LIBXML2_LIB. <>" LIBXML2_LIB)
CheckAndSetDefaultEnv(LIBXML2_INC "" PATH "Path to installed LIBXML2 headers. Overrides environment variable \$LIBXML2_INC. <>" LIBXML2_INC)
CheckAndSetDefaultEnv(LOG4CPP_LIB "" PATH "Path to pre-built LOG4CPP libraries. Overrides environment variable \$LOG4CPP_LIB. <>" LOG4CPP_LIB)
CheckAndSetDefaultEnv(LOG4CPP_INC "" PATH "Path to installed LOG4CPP headers. Overrides environment variable \$LOG4CPP_INC. <>" LOG4CPP_INC)
CheckAndSetDefaultEnv(PYTHIA6 "" PATH "Path to directory containing Pythia6 library. Overrides environment variable \$PYTHIA6. <>" PYTHIA6)

################################  LHAPDF  ######################################
if(LHAPDF_LIB STREQUAL "")
  message(FATAL_ERROR "Variable LHAPDF_LIB is not defined. The location of a pre-built lhapdf install must be defined either as $ cmake -DLHAPDF_LIB=/path/to/LHAPDF_libraries or as and environment vairable $ export LHAPDF_LIB=/path/to/LHAPDF_libraries")
endif()

if(LHAPDF_INC STREQUAL "")
  message(FATAL_ERROR "Variable LHAPDF_INC is not defined. The location of a pre-built lhapdf install must be defined either as $ cmake -DLHAPDF_INC=/path/to/LHAPDF_includes or as and environment vairable $ export LHAPDF_INC=/path/to/LHAPDF_includes")
endif()

if(LHAPATH STREQUAL "")
  message(FATAL_ERROR "Variable LHAPATH is not defined. The location of a the LHAPATH directory must be defined either as $ cmake -DLHAPATH=/path/to/LHAPATH or as and environment variable $ export LHAPATH=/path/to/LHAPATH")
endif()

################################  LIBXML  ######################################
if(LIBXML2_LIB STREQUAL "")
  message(FATAL_ERROR "Variable LIBXML2_LIB is not defined. The location of a pre-built libxml2 install must be defined either as $ cmake -DLIBXML2_LIB=/path/to/LIBXML2_libraries or as and environment vairable $ export LIBXML2_LIB=/path/to/LIBXML2_libraries")
endif()

if(LIBXML2_INC STREQUAL "")
  message(FATAL_ERROR "Variable LIBXML2_INC is not defined. The location of a pre-built libxml2 install must be defined either as $ cmake -DLIBXML2_INC=/path/to/LIBXML2_includes or as and environment vairable $ export LIBXML2_INC=/path/to/LIBXML2_includes")
endif()
###############################  log4cpp  ######################################
if(LOG4CPP_LIB STREQUAL "")
  find_program(LOG4CPPCFG log4cpp-config)
  if(NOT LOG4CPPCFG STREQUAL "LOG4CPPCFG-NOTFOUND")
    execute_process (COMMAND ${LOG4CPPCFG}
    --pkglibdir OUTPUT_VARIABLE LOG4CPP_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    message(FATAL_ERROR "Variable LOG4CPP_LIB is not defined. The location of a pre-built log4cpp install must be defined either as $ cmake -DLOG4CPP_LIB=/path/to/LOG4CPP_libraries or as and environment vairable $ export LOG4CPP_LIB=/path/to/LOG4CPP_libraries")
  endif()
endif()

if(LOG4CPP_INC  STREQUAL "")
  find_program(LOG4CPPCFG log4cpp-config)
  if(NOT LOG4CPPCFG STREQUAL "LOG4CPPCFG-NOTFOUND")
    execute_process (COMMAND ${LOG4CPPCFG}
    --pkgincludedir OUTPUT_VARIABLE LOG4CPP_INC OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    message(FATAL_ERROR "Variable LOG4CPP_INC is not defined. The location of a pre-built log4cpp install must be defined either as $ cmake -DGENIE_LOG4CPP_INC=/path/to/LOG4CPP_includes or as and environment vairable $ export LOG4CPP_INC=/path/to/LOG4CPP_includes")
  endif()
endif()

if(PYTHIA6 STREQUAL "")
  message(FATAL_ERROR "Variable PYTHIA6 is not defined. The location of a pre-built PYTHIA6 prefix must be defined either as $ cmake -DPYTHIA6=/path/to/PYTHIA6 or as and environment vairable $ export PYTHIA6=/path/to/PYTHIA6")
endif()
################################################################################

SET(GENIE_CXX_FLAGS -D__GENIE_ENABLED__)

SET(GENIE_INCLUDE_DIRS)
LIST(APPEND GENIE_INCLUDE_DIRS
  ${GENIE_INCLUDES_DIR}
  ${LHAPDF_INC}
  ${LIBXML2_INC}
  ${LOG4CPP_INC})

LIST(APPEND GENIE_LINK_DIRS
  ${GENIE_LIB_DIR}
  ${PYTHIA6}
  ${LHAPDF_LIB}
  ${LIBXML2_LIB}
  ${LOG4CPP_LIB})

LIST(APPEND GENIE_LIBS ${GENIE_LIBS_LIST})
LIST(APPEND GENIE_LIBS LHAPDF xml2 log4cpp)

if(NOT GENIE_LINK_DIRS STREQUAL "")
  string(REPLACE ";" " -L" STR_GENIE_LINK_DIRS "-L${GENIE_LINK_DIRS}")
endif()

if(NOT GENIE_LIBS STREQUAL "")
  string(REPLACE ";" " -l" STR_GENIE_LIBS "-l${GENIE_LIBS}")
endif()

message (STATUS "[GENIE]: Includes : ${GENIE_INCLUDE_DIRS} ")
message (STATUS "[GENIE]: Link dirs: ${STR_GENIE_LINK_DIRS} ")
message (STATUS "[GENIE]: Link libs: ${STR_GENIE_LIBS} ")

SET(GENIE_LIBS ${GENIE_LIBS};${GENIE_LIBS})
