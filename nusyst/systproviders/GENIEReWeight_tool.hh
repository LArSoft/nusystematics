#ifndef NUSYST_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN
#define NUSYST_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN

#include "nusyst/interface/IGENIESystProvider_tool.hh"

// GENIE
#include "EVGCore/EventRecord.h"
#include "ReWeight/GReWeight.h"

#include <memory>
#include <set>

class GENIEReWeight : public nusyst::IGENIESystProvider_tool {
public:
  explicit GENIEReWeight(fhicl::ParameterSet const &);

  bool Configure();

#ifndef NO_ART
  std::unique_ptr<larsyst::EventResponse> GetEventResponse(art::Event &);
  larsyst::SystMetaData ConfigureFromFHICL(fhicl::ParameterSet const &,
                                           paramId_t);
#endif

  larsyst::event_unit_response_t GetEventResponse(genie::EventRecord &);

  std::string AsString();

private:
  std::unique_ptr<genie::rew::GReWeight> GReWeightEngine;

  std::map<size_t, std::map<genie::rew::GSyst_t, size_t>>
      ResponseToGENIEParameters;

  void extend_ResponseToGENIEParameters(
      std::map<size_t, std::map<genie::rew::GSyst_t, size_t>> &&);
  std::set<genie::rew::GSyst_t> GENIEEngineDials;

#ifndef NO_ART
  std::string fGENIEModuleLabel = "generator";
#endif
};

#endif
