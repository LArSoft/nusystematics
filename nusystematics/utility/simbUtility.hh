#ifndef nusystematics_UTILITY_simbUtility_HH_SEEN
#define nusystematics_UTILITY_simbUtility_HH_SEEN

namespace nusyst {
  ///N.B. This is bad. It is a hard copy of simb::int_type from LArSoft/nusimdata/SimulationBase/MCNeutrino.h It means that standalone nusystematics has access to this enum as well
enum class simb_mode_copy {
  kUnknownInteraction = -1,
  kQE = 0,
  kRes = 1,
  kDIS = 2,
  kCoh = 3,
  kCohElastic = 4,
  kElectronScattering = 5,
  kIMDAnnihilation = 6,
  kInverseBetaDecay = 7,
  kGlashowResonance = 8,
  kAMNuGamma = 9,
  kMEC = 10,
  kDiffractive = 11,
  kEM = 12,
  kWeakMix = 13,
};
}

#endif
