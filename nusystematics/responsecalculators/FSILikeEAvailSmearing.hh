#ifndef nusystematics_RESPONSE_CALCULATORS_FSILIKEEAVAILSMEARING_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_FSILIKEEAVAILSMEARING_HH_SEEN

#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

namespace nusyst {

class TemplateResponseQ0Q3EAvail : public TemplateResponse3DDiscrete {
public:
  std::string GetCalculatorName() const { return "TemplateResponseQ0Q3EAvail"; }
};

typedef TemplateResponseQ0Q3EAvail FSILikeEAvailSmearing_ReWeight;
} // namespace nusyst

#endif
