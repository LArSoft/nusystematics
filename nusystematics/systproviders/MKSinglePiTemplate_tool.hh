#ifndef nusystematics_SYSTPROVIDERS_MKSinglePiTemplate_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_MKSinglePiTemplate_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/MKSinglePiTemplate_ReWeight.hh"

// GENIE includes
#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "Interaction/SppChannel.h"
#else
  // Use these for GENIE v3
  #include "Framework/Interaction/SppChannel.h"
#endif

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class MKSinglePiTemplate : public nusyst::IGENIESystProvider_tool {

  size_t ResponseParameterIdx;

  struct TemplateHelper {
    std::unique_ptr<nusyst::MKSinglePiTemplate_ReWeight> Template;
    bool ZeroIsValid;
  };

  std::map<genie::SppChannel_t, TemplateHelper>
      ChannelParameterMapping;

public:
  explicit MKSinglePiTemplate(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~MKSinglePiTemplate();

private:
  /// Whether input templates should be interpreted as in Enuq0q3 or EnuQ2W
  bool use_Q2W_templates;
  bool Q2_or_q0_is_x;

  fhicl::ParameterSet tool_options;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  bool SuppressNeutrinoBkgSPP;
  bool SuppressAntiNeutrinoBkgSPP;

  int NEUTMode, Pdgnu, pdgfslep, pdghmfspi, SppChannel, IsDIS;
  double Enu, momfslep, cthetafslep, momhmfspi, cthetahmfspi, Q2, q0, q3,
      Enu_nuc_rest_frame, q0_nuc_rest_frame, q3_nuc_rest_frame, W, weight;
};

#endif
