#ifndef NUSYST_SYSTPROVIDERS_MKSINGLEPIENUQ0Q3_TOOL_SEEN
#define NUSYST_SYSTPROVIDERS_MKSINGLEPIENUQ0Q3_TOOL_SEEN

#include "nusyst/interface/IGENIESystProvider_tool.hh"

#include "nusyst/responsecalculators/MKSinglePiEnuq0q3_ReWeight.hh"

// GENIE
#include "EVGCore/EventRecord.h"
#include "Interaction/SppChannel.h"

#include <memory>
#include <string>

class MKSinglePiEnuq0q3 : public nusyst::IGENIESystProvider_tool {

  std::unique_ptr<nusyst::MKSinglePiEnuq0q3_ReWeight> templateReweighter;

  larsyst::paramId_t ResponseParameterId;
  std::map<genie::SppChannel_t, larsyst::paramId_t> ChannelParameterMapping;

public:
  explicit MKSinglePiEnuq0q3(fhicl::ParameterSet const &);

  bool Configure();

#ifndef NO_ART
  std::unique_ptr<larsyst::EventResponse> GetEventResponse(art::Event &);
  larsyst::SystMetaData ConfigureFromFHICL(fhicl::ParameterSet const &,
                                           larsyst::paramId_t);
#endif

  larsyst::event_unit_response_t GetEventResponse(genie::EventRecord &);

  std::string AsString();
};

#endif
