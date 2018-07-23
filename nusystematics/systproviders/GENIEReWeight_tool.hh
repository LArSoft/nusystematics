#ifndef nusystematics_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/GENIEResponseParameterAssociation.hh"

// GENIE
#include "EVGCore/EventRecord.h"
#include "ReWeight/GReWeight.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <set>

// HERG: HIRD OF RAMPAGING GENIES, HIRD: HERG OF INFINITELY REPEATING DEPTHS

class GENIEReWeight : public nusyst::IGENIESystProvider_tool {
public:
  explicit GENIEReWeight(fhicl::ParameterSet const &);

#ifndef NO_ART
  std::unique_ptr<systtools::EventResponse> GetEventResponse(art::Event &);
#endif

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord &);

  std::string AsString();

  ~GENIEReWeight();

private:
  std::vector<nusyst::GENIEResponseParameter> ResponseToGENIEParameters;

  void extend_ResponseToGENIEParameters(
      std::vector<nusyst::GENIEResponseParameter> &&);

  fhicl::ParameterSet tool_options;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu;
  double Enu, Q2, q0, q3, W;
  std::vector<double> weights;

#ifndef NO_ART
  std::string fGENIEModuleLabel = "generator";
#endif
};

#endif
