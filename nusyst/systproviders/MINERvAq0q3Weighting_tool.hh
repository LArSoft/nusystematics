#ifndef NUSYST_SYSTPROVIDERS_MINERVAQ0Q3_TOOL_SEEN
#define NUSYST_SYSTPROVIDERS_MINERVAQ0Q3_TOOL_SEEN

#include "nusyst/interface/IGENIESystProvider_tool.hh"

#include "nusyst/responsecalculators/MINERvARPAq0q3_ReWeight.hh"
#include "nusyst/responsecalculators/MINERvAq0q3Weighting_data.hh"

#include "nusyst/utility/GENIEUtils.hh"

// GENIE
#include "EVGCore/EventRecord.h"

#include <memory>
#include <string>

///
///
/// Can apply the MINERvA RPA tune and two discrete 'tweaks', the 2p2h gaussian
/// enhancement tune and 3 discrete 'tweaks', can also apply the 2p2h gaussian
/// with arbitrary parameters.
class MINERvAq0q3Weighting : public nusyst::IGENIESystProvider_tool {

  enum class param_t {
    kGaussNorm = 0,
    kGaussMeanQ0,
    kGaussMeanQ3,
    kGaussSigmaQ0,
    kGaussSigmaQ3,
    kGaussCorrelation,
    kGaussResponse,
    kMINERvARPA,
    kMINERvA2p2h
  };

  struct param_t_name {
    // The parameter prettyName expected to be found in the input metadata
    std::string name;
    // The local enum id of the parameter
    param_t lid;
  };

  std::array<double, 6> Gauss2DParams = nusyst::Gauss2DParams_CV;
  std::unique_ptr<nusyst::MINERvARPAq0q3_ReWeight> RPATemplateReweighter;
  std::map<param_t, larsyst::paramId_t> ConfiguredParameters;

public:
  explicit MINERvAq0q3Weighting(fhicl::ParameterSet const &);

  bool Configure();

#ifndef NO_ART
  std::unique_ptr<larsyst::EventResponse> GetEventResponse(art::Event &);
  larsyst::SystMetaData ConfigureFromFHICL(fhicl::ParameterSet const &,
                                           larsyst::paramId_t);
#endif

  double GetMINERvARPATuneWeight(double val, double q0, double q3);
  double GetMINERvA2p2hTuneWeight(double val, double q0, double q3,
                                  nusyst::QELikeTarget_t QELTarget);

  larsyst::event_unit_response_t GetEventResponse(genie::EventRecord &);

  std::string AsString();
};

#endif
