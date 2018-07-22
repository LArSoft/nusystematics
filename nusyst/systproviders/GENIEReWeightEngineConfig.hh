#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN

#include "larsyst/interface/SystMetaData.hh"

#include "nusyst/systproviders/GENIEResponseParameterAssociation.hh"

#include "fhiclcpp/ParameterSet.h"

// GENIE
#include "ReWeight/GReWeight.h"

#include <map>
#include <memory>

namespace nusyst {

std::vector<GENIEResponseParameter>
ConfigureQEWeightEngine(larsyst::SystMetaData const &,
                        fhicl::ParameterSet const &tool_options);

#ifndef GRWTEST

std::vector<GENIEResponseParameter>
ConfigureNCELWeightEngine(larsyst::SystMetaData const &,
                          fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureRESWeightEngine(larsyst::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureCOHWeightEngine(larsyst::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureDISWeightEngine(larsyst::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureFSIWeightEngine(larsyst::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureOtherWeightEngine(larsyst::SystMetaData const &,
                           fhicl::ParameterSet const &tool_options);

#endif

} // namespace nusyst

#endif
