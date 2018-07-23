#ifndef nusystematics_SYSTPROVIDERS_MKSINGLEPIENUQ0Q3_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_MKSINGLEPIENUQ0Q3_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/MKSinglePiEnuq0q3_ReWeight.hh"

// GENIE
#include "EVGCore/EventRecord.h"
#include "Interaction/SppChannel.h"

#include <memory>
#include <string>

class MKSinglePiEnuq0q3 : public nusyst::IGENIESystProvider_tool {

  std::unique_ptr<nusyst::MKSinglePiEnuq0q3_ReWeight> templateReweighter;

  systtools::paramId_t ResponseParameterId;
  std::map<genie::SppChannel_t, systtools::paramId_t> ChannelParameterMapping;

public:
  explicit MKSinglePiEnuq0q3(fhicl::ParameterSet const &);

#ifndef NO_ART
  std::unique_ptr<systtools::EventResponse> GetEventResponse(art::Event &);
#endif

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord &);

  std::string AsString();

private:
  fhicl::ParameterSet tool_options;
};

#endif
