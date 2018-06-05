#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHTENGINECONFIG_SEEN

#include "larsyst/interface/SystMetaData.hh"

//GENIE
#include "ReWeight/GReWeight.h"

#include <map>
#include <memory>

namespace nusyst {

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureQEWeightEngine(larsyst::SystMetaData const &,
                        std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureNCELWeightEngine(larsyst::SystMetaData const &,
                         std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureRESWeightEngine(larsyst::SystMetaData const &,
                         std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureCOHWeightEngine(larsyst::SystMetaData const &,
                         std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureDISWeightEngine(larsyst::SystMetaData const &,
                         std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureFSIWeightEngine(larsyst::SystMetaData const &,
                         std::unique_ptr<genie::rew::GReWeight> &);

std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
ConfigureOtherWeightEngine(larsyst::SystMetaData const &,
                           std::unique_ptr<genie::rew::GReWeight> &);
} // namespace nusyst

#endif
