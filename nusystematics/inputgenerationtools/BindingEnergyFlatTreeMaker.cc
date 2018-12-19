#include "systematicstools/utility/string_parsers.hh"

// GENIE includes

#ifdef GENIE_PRE_R3
  // Use these for GENIE v2
  #include "EVGCore/EventRecord.h"
  #include "GHEP/GHepParticle.h"
  #include "GHEP/GHepUtils.h"
  #include "Messenger/Messenger.h"
  #include "Ntuple/NtpMCEventRecord.h"
#else
  // Use these for GENIE v3
  #include "Framework/EventGen/EventRecord.h"
  #include "Framework/GHEP/GHepParticle.h"
  #include "Framework/GHEP/GHepUtils.h"
  #include "Framework/Messenger/Messenger.h"
  #include "Framework/Ntuple/NtpMCEventRecord.h"
#endif

#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace systtools;

namespace cliopts {
std::string genie_input = "";
std::string outputfile = "";
size_t NMax = std::numeric_limits<size_t>::max();
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-o <out.root>    : File to write validation canvases to.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::genie_input = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);
  if (!cliopts::genie_input.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }
  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(cliopts::genie_input.c_str())) {
    std::cout << "[ERROR]: Failed to find any TTrees named "
              << std::quoted("gtree") << ", from TChain::Add descriptor: "
              << std::quoted(cliopts::genie_input) << "." << std::endl;
    return 3;
  }

  size_t NEvs = gevs->GetEntries();

  if (!NEvs) {
    std::cout << "[ERROR]: Input TChain contained no entries." << std::endl;
    return 4;
  }

  genie::NtpMCEventRecord *GenieNtpl = nullptr;

  if (gevs->SetBranchAddress("gmcrec", &GenieNtpl) != TTree::kMatch) {
    std::cout << "[ERROR]: Failed to set branch address on ghep tree."
              << std::endl;
    return 5;
  }

  TFile *outf = new TFile(cliopts::outputfile.c_str(), "RECREATE");

  TTree *outt = new TTree("EBFlatTree", "");

  double Enu;
  double q0, q3;
  double cthetamu;
  double pmu;

  outt->Branch("Enu", &Enu, "Enu/D");
  outt->Branch("q0", &q0, "q0/D");
  outt->Branch("q3", &q3, "q3/D");
  outt->Branch("cthetamu", &cthetamu, "cthetamu/D");
  outt->Branch("pmu", &pmu, "pmu/D");

  size_t NToRead = std::min(NEvs, cliopts::NMax);
  size_t NToShout = NToRead / 20;
  NToShout = NToShout ? NToShout : 1;
  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;

    if (!(ev_it % NToShout)) {
      std::cout << (ev_it ? "\r" : "") << "Event #" << ev_it << "/" << NToRead
                << ", Interaction: " << GenieGHep.Summary()->AsString()
                << std::flush;
    }

    if (!(GenieGHep.Summary()->ProcInfo().IsQuasiElastic() &&
          GenieGHep.Summary()->ProcInfo().IsWeakCC())) {
      continue;
    }

    genie::GHepParticle *FSLep = GenieGHep.FinalStatePrimaryLepton();
    genie::GHepParticle *ISLep = GenieGHep.Probe();
    TLorentzVector FSLepP4 = *FSLep->P4();
    TLorentzVector ISLepP4 = *ISLep->P4();

    TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
    Enu = ISLepP4.E();
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();
    cthetamu = FSLepP4.Vect().CosTheta();
    pmu = FSLepP4.Vect().Mag();
    outt->Fill();

  }
  std::cout << std::endl;

  outf->Write();
}
