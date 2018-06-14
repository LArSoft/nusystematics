###### FHICL set up
include(ExternalProject)

ExternalProject_Add(larsyst
  PREFIX "${PROJECT_BINARY_DIR}/larsyst"
  GIT_REPOSITORY https://github.com/luketpickering/larsyst.git
  # GIT_REPOSITORY file:///home/picker24/software/larsyst/larsyst
  GIT_TAG develop
  CMAKE_ARGS -DUSEART=0
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})

###### ROOT set up

if (NOT DEFINED ENV{ROOTSYS})
  message (FATAL_ERROR "$ROOTSYS is not defined, please set up ROOT first.")
else()
  message(STATUS "Using ROOT installed at $ENV{ROOTSYS}")
  set(CMAKE_ROOTSYS $ENV{ROOTSYS})
endif()

# Get cflags from ROOT
execute_process (COMMAND root-config
  --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE " " ";" ROOT_CXX_FLAGS_PREPARSE "${ROOT_CXX_FLAGS_RAW}")
set(ROOT_CXX_FLAGS_SCRUBBED "")
foreach(FLG ${ROOT_CXX_FLAGS_PREPARSE})
  string(REGEX MATCH "-std=" ISSTDFLAG ${FLG})
  if(NOT ISSTDFLAG) #ignore -std= flags from ROOT config
    SET(ROOT_CXX_FLAGS_SCRUBBED "${ROOT_CXX_FLAGS_SCRUBBED} ${FLG}")
  endif()
endforeach()
# Get libdir from ROOT
execute_process (COMMAND root-config
  --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get version from ROOT
execute_process (COMMAND root-config
  --version OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get features from ROOT
execute_process (COMMAND root-config
  --features OUTPUT_VARIABLE ROOT_FEATURES OUTPUT_STRIP_TRAILING_WHITESPACE)

SET(ROOT_LIBS)
LIST(APPEND ROOT_LIBS
  Core
  RIO
  XMLIO
  Net
  Hist
  Graf
  Graf3d
  Gpad
  Tree
  Rint
  Postscript
  Matrix
  Physics
  MathCore
  Thread
  EG
  Geom
  GenVector
  EGPythia6
  Pythia6)

message (STATUS "[ROOT]: root-config --version: ${ROOT_VERSION} ")
message (STATUS "[ROOT]: root-config --cflags : ${ROOT_CXX_FLAGS_SCRUBBED} ")
message (STATUS "[ROOT]: libs in use          : ${ROOT_LIBS} ")

###### GENIE setup
  include(${CMAKE_CURRENT_SOURCE_DIR}/GENIESetup.cmake)
###### Compiler set up
#-Wl,--no-as-needed used to ensuring dynamic linking of GENIE algorithm libraries so that TClass::GetClass has access to the dictionaries at runtime
set(EXTRA_LINK_FLAGS "-Wl,--no-as-needed -L${CMAKE_INSTALL_PREFIX}/lib ${STR_GENIE_LINK_DIRS} -L${ROOT_LIBDIR}")
if(CMAKE_LINK_FLAGS STREQUAL "")
  SET(CMAKE_LINK_FLAGS "${EXTRA_LINK_FLAGS}")
else()
  SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${EXTRA_LINK_FLAGS}")
endif()

set(EXTRA_CXX_FLAGS "-DNO_ART -std=c++1y ${ROOT_CXX_FLAGS_SCRUBBED} ${GENIE_CXX_FLAGS}")
if(CMAKE_CXX_FLAGS STREQUAL "")
  SET(CMAKE_CXX_FLAGS ${EXTRA_CXX_FLAGS})
else()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()

include_directories(${GENIE_INCLUDES_DIR})

include_directories(${CMAKE_INSTALL_PREFIX}/include)
include_directories(${CMAKE_SOURCE_DIR})


get_property(INCDIRS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${INCDIRS})
  message(STATUS "INCLUDE: ${dir}")
endforeach()

####### systproviders library
SET(SP_IMPLFILES
  ${CMAKE_SOURCE_DIR}/nusyst/systproviders/GENIEReWeight_tool.cc
  ${CMAKE_SOURCE_DIR}/nusyst/systproviders/GENIEReWeightEngineConfig.cc)

SET(SP_HDRFILES
  ${CMAKE_SOURCE_DIR}/nusyst/systproviders/GENIEReWeight_tool.hh)

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
  add_library(nusyst_systproviders STATIC ${SP_IMPLFILES})
else(CMAKE_BUILD_TYPE MATCHES RELEASE)
  add_library(nusyst_systproviders SHARED ${SP_IMPLFILES})
endif()

INSTALL(TARGETS nusyst_systproviders DESTINATION lib)
INSTALL(FILES ${SP_HDRFILES} DESTINATION include/nusyst/systproviders)

add_executable(ReWeightGHep ReWeightGHep.cc)
add_dependencies(ReWeightGHep larsyst)
target_link_libraries(ReWeightGHep nusyst_systproviders)

set_target_properties(ReWeightGHep PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
target_link_libraries(ReWeightGHep
  larsyst_utility
  larsyst_interpreters)
target_link_libraries(ReWeightGHep ${GENIE_LIBS})
target_link_libraries(ReWeightGHep ${ROOT_LIBS})

INSTALL(TARGETS ReWeightGHep DESTINATION bin)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})
