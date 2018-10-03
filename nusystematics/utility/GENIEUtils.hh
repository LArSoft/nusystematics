#ifndef nusystematics_UTILITY_GENIEUTILS_SEEN
#define nusystematics_UTILITY_GENIEUTILS_SEEN

#include "EVGCore/EventRecord.h"

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
  case 11: { // kCCSPP_PPip
    return genie::kSpp_vp_cc_10100;
  }
  case 12: { // kCCSPP_PPi0
    return genie::kSpp_vn_cc_10010;
  }
  case 13: { // kCCSPP_NPip
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
NEW_SYSTTOOLS_EXCEPT(indeterminable_QELikeTarget);

inline QELikeTarget_t GetQELikeTarget(genie::EventRecord const &ev) {

  if (ev.Summary()->ProcInfo().IsQuasiElastic()) {
    return QELikeTarget_t::kQE;
  }
  if (ev.Summary()->ProcInfo().IsMEC()) {
    genie::Target const &tgt = ev.Summary()->InitState().Tgt();
    size_t nuc_pdg = tgt.HitNucPdg();

    switch (nuc_pdg - 2000000200) {
    case 0:
    case 2: {
      return QELikeTarget_t::kNN;
    }
    case 1: {
      return QELikeTarget_t::knp;
    }
    default: {
      throw indeterminable_QELikeTarget()
          << "[ERROR]: Failed to determine 2p2h topology from interaction: "
          << ev.Summary()->AsString()
          << ", expected the target nucleon PDG: " << nuc_pdg
          << " - 2000000200 = " << (nuc_pdg - 2000000200)
          << " to be 0, 1, or 2.";
    }
    }
  }

  return QELikeTarget_t::kInvalidTopology;
}

} // namespace nusyst

#endif
