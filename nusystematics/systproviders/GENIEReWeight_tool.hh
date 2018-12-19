#ifndef nusystematics_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_GENIEREWEIGHT_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/GENIEResponseParameterAssociation.hh"

// GENIE includes
#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "ReWeight/GReWeight.h"
#else
  // Use these for GENIE v3
  #include "RwFramework/GReWeight.h"
#endif

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <set>

// HERG: HIRD OF RAMPAGING GENIES, HIRD: HERG OF INFINITELY REPEATING DEPTH

class GENIEReWeight : public nusyst::IGENIESystProvider_tool {
public:
  NEW_SYSTTOOLS_EXCEPT(invalid_engine_state);

  explicit GENIEReWeight(fhicl::ParameterSet const &);

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  double GetEventWeightResponse(genie::EventRecord const &,
                                systtools::param_value_list_t const &);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &,
                                                    systtools::paramId_t);

  std::string AsString();

  ~GENIEReWeight();

private:
  /// Set when GetEventWeightResponse has been used as it reconfigures weight
  /// engines such that GetEventResponse will not perform as expected.
  bool fHaveReconfiguredOneOfTheHERG;

  systtools::ParamResponses
  GetEventGENIEParameterResponse(genie::EventRecord const &, size_t idx);

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
};

#endif
