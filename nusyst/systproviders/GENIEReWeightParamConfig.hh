#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN

#include "larsyst/interface/SystMetaData.hh"

#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nusyst {

larsyst::SystMetaData
ConfigureQEParameterHeaders(fhicl::ParameterSet const &, larsyst::paramId_t,
                            fhicl::ParameterSet &tool_options);

#ifndef GRWTEST

larsyst::SystMetaData ConfigureNCELParameterHeaders(fhicl::ParameterSet const &,
                                                    larsyst::paramId_t);

larsyst::SystMetaData ConfigureRESParameterHeaders(fhicl::ParameterSet const &,
                                                   larsyst::paramId_t);

larsyst::SystMetaData ConfigureCOHParameterHeaders(fhicl::ParameterSet const &,
                                                   larsyst::paramId_t);

larsyst::SystMetaData ConfigureDISParameterHeaders(fhicl::ParameterSet const &,
                                                   larsyst::paramId_t);

larsyst::SystMetaData ConfigureFSIParameterHeaders(fhicl::ParameterSet const &,
                                                   larsyst::paramId_t);
larsyst::SystMetaData
ConfigureOtherParameterHeaders(fhicl::ParameterSet const &, larsyst::paramId_t);

#endif

} // namespace nusyst

#endif
