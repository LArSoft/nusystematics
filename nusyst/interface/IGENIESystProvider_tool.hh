#ifndef NUSYST_INTERFACE_IGENIESYSTPROVIDER_TOOL_SEEN
#define NUSYST_INTERFACE_IGENIESYSTPROVIDER_TOOL_SEEN

#include "larsyst/interface/ISystProvider_tool.hh"

#include "fhiclcpp/ParameterSet.h"

// GENIE
#include "EVGCore/EventRecord.h"

namespace nusyst {
class IGENIESystProvider_tool : public larsyst::ISystProvider_tool {
public:
  IGENIESystProvider_tool(fhicl::ParameterSet const &ps)
      : ISystProvider_tool(ps) {}
  virtual larsyst::event_unit_response_t
  GetEventResponse(genie::EventRecord &) = 0;
};
} // namespace nusyst

#endif
