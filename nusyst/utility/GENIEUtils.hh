#ifndef NUSYST_UTILITY_GENIEUTILS_SEEN
#define NUSYST_UTILITY_GENIEUTILS_SEEN

#include "GHEP/GHepUtils.h"
#include "Interaction/SppChannel.h"

namespace nusyst {
/// Gets the GENIE SPP channel enum for a supplied GHepEvent
///
/// N.B. going via NEUT mode as genie::SppChannel::FromInteraction doesn't work
/// for RES events.
inline genie::SppChannel_t SPPChannelFromGHep(genie::EventRecord const &ev) {
  int NEUTCh = genie::utils::ghep::NeutReactionCode(&ev);

  switch (NEUTCh) {
  case 11: {
    return genie::kSpp_vp_cc_10100;
  }
  case 12: {
    return genie::kSpp_vn_cc_10010;
  }
  case 13: {
    return genie::kSpp_vn_cc_01100;
  }
  case -11: {
    return genie::kSpp_vbn_cc_01001;
  }
  case -12: {
    return genie::kSpp_vbp_cc_01010;
  }
  case -13: {
    return genie::kSpp_vbp_cc_10001;
  }
  default: { return genie::kSppNull; }
  }
}

enum class QELikeTarget_t { kNN = 0, knp, kQE, kInvalidTopology };

inline QELikeTarget_t GetQELikeTarget(genie::EventRecord const &ev) {

  if (ev.Summary()->ProcInfo().IsQuasiElastic()) {
    return QELikeTarget_t::kQE;
  }
  if (ev.Summary()->ProcInfo().IsMEC()) {
    return QELikeTarget_t::kNN;
  }
  return QELikeTarget_t::kInvalidTopology;
}

} // namespace nusyst

#endif
