#ifndef nusystematics_RESPONSE_CALCULATORS_MKSinglePiTemplate_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_MKSinglePiTemplate_REWEIGHT_HH_SEEN

#include "nusystematics/responsecalculators/EnuBinnedTemplateResponseCalculator.hh"
#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

namespace nusyst {

class TemplateResponseQ0Q3 : public TemplateResponse2DDiscrete {
public:
  std::string GetCalculatorName() const { return "TemplateResponseQ0Q3"; }
};

typedef EnuBinnedTemplateResponseCalculator<TemplateResponseQ0Q3>
    MKSinglePiTemplate_ReWeight;
} // namespace nusyst

#endif
