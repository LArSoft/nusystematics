#ifndef nusystematics_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN
#define nusystematics_SYSTPROVIDERS_GENIEREWEIGHTPARAMCONFIG_SEEN

#include "systematicstools/interface/SystMetaData.hh"

#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nusyst {

systtools::SystMetaData
ConfigureQEParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                            fhicl::ParameterSet &tool_options);

#ifndef GRWTEST

systtools::SystMetaData ConfigureNCELParameterHeaders(fhicl::ParameterSet const &,
                                                    systtools::paramId_t);

systtools::SystMetaData ConfigureRESParameterHeaders(fhicl::ParameterSet const &,
                                                   systtools::paramId_t);

systtools::SystMetaData ConfigureCOHParameterHeaders(fhicl::ParameterSet const &,
                                                   systtools::paramId_t);

systtools::SystMetaData ConfigureDISParameterHeaders(fhicl::ParameterSet const &,
                                                   systtools::paramId_t);

systtools::SystMetaData ConfigureFSIParameterHeaders(fhicl::ParameterSet const &,
                                                   systtools::paramId_t);
systtools::SystMetaData
ConfigureOtherParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t);

#endif

} // namespace nusyst

#endif
