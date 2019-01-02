#ifndef nusystematics_UTILITY_GENIEUTILS_SEEN
#define nusystematics_UTILITY_GENIEUTILS_SEEN

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/utility/simbUtility.hh"

// GENIE includes
#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "EVGCore/EventRecord.h"
  #include "GHEP/GHepParticle.h"
  #include "GHEP/GHepUtils.h"
  #include "Interaction/SppChannel.h"
#else
  // Use these for GENIE v3
  #include "Framework/EventGen/EventRecord.h"
  #include "Framework/GHEP/GHepParticle.h"
  #include "Framework/GHEP/GHepUtils.h"
  #include "Framework/Interaction/SppChannel.h"
#endif

#include <sstream>

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

  if (ev.Summary()->ProcInfo().IsQuasiElastic() &&
      !ev.Summary()->ExclTag().IsCharmEvent()) {
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

/// Encoded as IsNeutrino * 1 + IsCC * 10 + TargetIsProton * 100 + NPi * 1000 +
/// NPiplus * 10000 + NPiminus * 100000 + NPi0 * 1000000
///
///\n.b. is TargetIsProton == 2 then it can be either proton or neutron target
typedef size_t NRPiChan_t;

inline bool IsNeutrinoNRPiChan(NRPiChan_t ch) { return ch % 10; }
inline bool IsCCNRPiChan(NRPiChan_t ch) { return (ch / 10) % 10; }
inline bool IsProtonTargetNRPiChan(NRPiChan_t ch) {
  return ((ch / 100) % 10) & 2;
}
inline bool IsNeutronTargetNRPiChan(NRPiChan_t ch) {
  return ((ch / 100) % 10) & 1;
}
inline size_t GetNRPiChanNPi(NRPiChan_t ch) { return (ch / 1000) % 10; }
inline size_t GetNRPiChanNPip(NRPiChan_t ch) { return (ch / 10000) % 10; }
inline size_t GetNRPiChanNPim(NRPiChan_t ch) { return (ch / 100000) % 10; }
inline size_t GetNRPiChanNPi0(NRPiChan_t ch) { return (ch / 1000000) % 10; }

inline NRPiChan_t BuildNRPiChannel(bool IsNeutrino, bool IsCC,
                                   int TargetNucleon, size_t NPi,
                                   size_t NPiplus = 0, size_t NPiminus = 0,
                                   size_t NPi0 = 0) {
  return IsNeutrino * 1 + IsCC * 10 + TargetNucleon * 100 + NPi * 1000 +
         NPiplus * 10000 + NPiminus * 100000 + NPi0 * 1000000;
}

inline NRPiChan_t GetNRPiChannel(genie::EventRecord const &ev) {

  if (!ev.Summary()->ProcInfo().IsDeepInelastic()) {
    return 0;
  }

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  if (!tgt.HitNucIsSet()) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to get hit nucleon kinematics as it was not "
           "included in this GHep event. This is a fatal error.";
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }
  int NPi0 = 0, NPip = 0, NPim = 0;

  bool nuclear_target = tgt.IsNucleus();

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    genie::GHepStatus_t ghep_ist = (genie::GHepStatus_t)p->Status();
    int ghep_pdgc = p->Pdg();
    int ghep_fm = p->FirstMother();
    int ghep_fmpdgc = (ghep_fm == -1) ? 0 : ev.Particle(ghep_fm)->Pdg();

    // For nuclear targets use hadrons marked as 'hadron in the nucleus'
    // which are the ones passed in the intranuclear rescattering
    // For free nucleon targets use particles marked as 'final state'
    // but make an exception for decayed pi0's,eta's (count them and not their
    // daughters)

    bool decayed =
        (ghep_ist == genie::kIStDecayedState &&
         (ghep_pdgc == genie::kPdgPi0 || ghep_pdgc == genie::kPdgEta));
    bool parent_included =
        (ghep_fmpdgc == genie::kPdgPi0 || ghep_fmpdgc == genie::kPdgEta);

    bool count_it =
        (nuclear_target && ghep_ist == genie::kIStHadronInTheNucleus) ||
        (!nuclear_target && decayed) ||
        (!nuclear_target && ghep_ist == genie::kIStStableFinalState &&
         !parent_included);

    if (!count_it) {
      continue;
    }
    if (ghep_pdgc == genie::kPdgPiP) {
      NPip++;
    } else if (ghep_pdgc == genie::kPdgPiM) {
      NPim++;
    } else if (ghep_pdgc == genie::kPdgPi0) {
      NPi0++;
    }
  }

  return BuildNRPiChannel(ISLep->Pdg() > 0, ev.Summary()->ProcInfo().IsWeakCC(),
                          (tgt.HitNucPdg() == genie::kPdgProton) ? 2 : 1,
                          NPi0 + NPip + NPim, NPip, NPim, NPi0);
}

inline std::string GetNRPiChannelName(NRPiChan_t ch) {
  std::stringstream ss("");

  ss << "NR_" << (IsNeutrinoNRPiChan(ch) ? "nu" : "nubar") << "_"
     << (IsNeutronTargetNRPiChan(ch) ? "n" : "")
     << (IsProtonTargetNRPiChan(ch) ? "p" : "") << "_"
     << (IsCCNRPiChan(ch) ? "CC" : "NC") << "_" << GetNRPiChanNPi(ch) << "Pi";

  return ss.str();
}

/// Check if an event channel is equivalent to an overall channel, n and p
/// target events are equivalent to np channels and NPions > NMax pions are
/// equivalent to NPions == NMax pions
inline bool ChannelsAreEquivalent(NRPiChan_t ch, NRPiChan_t event_ch,
                                  size_t NMaxPions) {
  if (IsNeutrinoNRPiChan(ch) != IsNeutrinoNRPiChan(event_ch)) {
    return false;
  }
  if (IsCCNRPiChan(ch) != IsCCNRPiChan(event_ch)) {
    return false;
  }
  if (!((IsProtonTargetNRPiChan(event_ch) && IsProtonTargetNRPiChan(ch)) ||
        (IsNeutronTargetNRPiChan(event_ch) && IsNeutronTargetNRPiChan(ch)))) {
    return false;
  }
  size_t NPiChannel = (GetNRPiChanNPi(event_ch) > NMaxPions)
                          ? NMaxPions
                          : GetNRPiChanNPi(event_ch);
  if (NPiChannel != GetNRPiChanNPi(ch)) {
    return false;
  }
  return true;
}

inline double GetErecoil_MINERvA_LowRecoil(genie::EventRecord const &ev) {
  // Get total energy of hadronic system.
  double Erecoil = 0.0;

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    if (p->Status() != genie::kIStStableFinalState) {
      continue;
    }
    switch (p->Pdg()) {
    case 2212:
    case 211:
    case -211: {
      Erecoil += p->KinE();
      break;
    }
    case 111:
    case 11:
    case -11:
    case -22: {
      Erecoil += p->E();
      break;
    }
    default: {}
    }
  }

  // For nue CC scattering, we would have counted the E of the charged lepton,
  // subtract it off here
  if (ev.Summary()->ProcInfo().IsWeakCC() && (abs(ev.Probe()->Pdg()) == 12)) {
    Erecoil -= ev.FinalStatePrimaryLepton()->P4()->E();
  }

  return Erecoil;
}

inline simb_mode_copy GetSimbMode(genie::EventRecord const &ev) {

  simb_mode_copy mode = simb_mode_copy::kUnknownInteraction;

  if (ev.Summary()->ProcInfo().IsQuasiElastic()) {
    mode = simb_mode_copy::kQE;
  } else if (ev.Summary()->ProcInfo().IsMEC()) {
    mode = simb_mode_copy::kMEC;
  } else if (ev.Summary()->ProcInfo().IsResonant()) {
    mode = simb_mode_copy::kRes;
  } else if (ev.Summary()->ProcInfo().IsDeepInelastic()) {
    mode = simb_mode_copy::kDIS;
  } else if (ev.Summary()->ProcInfo().IsCoherent()) {
    mode = simb_mode_copy::kCoh;
  } else if (ev.Summary()->ProcInfo().IsCoherentElas()) {
    mode = simb_mode_copy::kCohElastic;
  } else if (ev.Summary()->ProcInfo().IsElectronScattering()) {
    mode = simb_mode_copy::kElectronScattering;
  } else if (ev.Summary()->ProcInfo().IsIMDAnnihilation()) {
    mode = simb_mode_copy::kIMDAnnihilation;
  } else if (ev.Summary()->ProcInfo().IsInverseBetaDecay()) {
    mode = simb_mode_copy::kInverseBetaDecay;
  } else if (ev.Summary()->ProcInfo().IsGlashowResonance()) {
    mode = simb_mode_copy::kGlashowResonance;
  } else if (ev.Summary()->ProcInfo().IsAMNuGamma()) {
    mode = simb_mode_copy::kAMNuGamma;
  } else if (ev.Summary()->ProcInfo().IsDiffractive()) {
    mode = simb_mode_copy::kDiffractive;
  } else if (ev.Summary()->ProcInfo().IsEM()) {
    mode = simb_mode_copy::kEM;
  } else if (ev.Summary()->ProcInfo().IsWeakMix()) {
    mode = simb_mode_copy::kWeakMix;
  }

  return mode;
}

inline std::string DumpGENIEEv(genie::EventRecord const &ev) {
  std::stringstream ss("");
  ss << ev.Summary()->AsString() << std::endl;

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();

  bool nuclear_target = tgt.IsNucleus();

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;
  size_t p_it = 0;

  while ((p = dynamic_cast<genie::GHepParticle *>(event_iter.Next()))) {
    genie::GHepStatus_t ghep_ist = (genie::GHepStatus_t)p->Status();
    int ghep_pdgc = p->Pdg();
    int ghep_fm = p->FirstMother();
    int ghep_fmpdgc = (ghep_fm == -1) ? 0 : ev.Particle(ghep_fm)->Pdg();

    bool decayed =
        (ghep_ist == genie::kIStDecayedState &&
         (ghep_pdgc == genie::kPdgPi0 || ghep_pdgc == genie::kPdgEta));
    bool parent_included =
        (ghep_fmpdgc == genie::kPdgPi0 || ghep_fmpdgc == genie::kPdgEta);

    bool count_it =
        (nuclear_target && ghep_ist == genie::kIStHadronInTheNucleus) ||
        (!nuclear_target && decayed) ||
        (!nuclear_target && ghep_ist == genie::kIStStableFinalState &&
         !parent_included);

    ss << "Part: " << p_it++ << ", pdg = " << ghep_pdgc << ", decayed ? "
       << decayed << ", parent_included ? " << parent_included
       << ", counting ? " << count_it << std::endl;
  }
  ss << std::endl;
  return ss.str();
}

} // namespace nusyst

#endif
