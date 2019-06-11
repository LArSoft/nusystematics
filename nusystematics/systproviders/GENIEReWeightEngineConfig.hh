#ifndef nusystematics_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN
#define nusystematics_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN

#include "systematicstools/interface/SystMetaData.hh"

#include "nusystematics/systproviders/GENIEResponseParameterAssociation.hh"

#include "fhiclcpp/ParameterSet.h"

// GENIE
#include "RwFramework/GReWeight.h"

#include <map>
#include <memory>

namespace nusyst {

std::vector<GENIEResponseParameter>
ConfigureQEWeightEngine(systtools::SystMetaData const &,
                        fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureNCELWeightEngine(systtools::SystMetaData const &,
                          fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureRESWeightEngine(systtools::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureCOHWeightEngine(systtools::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureDISWeightEngine(systtools::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureFSIWeightEngine(systtools::SystMetaData const &,
                         fhicl::ParameterSet const &tool_options);

std::vector<GENIEResponseParameter>
ConfigureOtherWeightEngine(systtools::SystMetaData const &,
                           fhicl::ParameterSet const &tool_options);

} // namespace nusyst

#endif
