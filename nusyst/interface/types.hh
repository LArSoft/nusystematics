#ifndef NUSYST_INTERFACE_TYPES_SEEN
#define NUSYST_INTERFACE_TYPES_SEEN

#include <map>
#include <memory>

namespace nusyst {

class IGENIESystProvider_tool;

typedef std::map<std::string, std::unique_ptr<IGENIESystProvider_tool>>
    genie_provider_map_t;

} // namespace nusyst

#endif
