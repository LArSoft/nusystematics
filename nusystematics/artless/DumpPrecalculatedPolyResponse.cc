#include "nusystematics/artless/response_helper.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "systematicstools/interpreters/PrecalculatedResponseReader.hh"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepUtils.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"

// Included by fhiclcpp-simple will not be available in art-ful
#include "string_parsers/from_string.hxx"
#include "string_parsers/to_string.hxx"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace fhicl;
using namespace systtools;
using namespace nusyst;

namespace cliopts {
std::string fclname = "";
std::string inputfile = "";
std::string outputfile = "";
size_t NMax = std::numeric_limits<size_t>::max();

} // namespace cliopts

constexpr size_t Order = 5;
constexpr size_t NCoeffs = Order + 1;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help         : Show this message.\n"
               "\t-c <config.fcl>   : fhicl file to read.\n"
               "\t-i <ghep.root>    : GENIE event file to read.\n"
               "\t-o <output.rooot> : Response file to write.\n"
               "\t-n <NMax>         : Only calculate splines for the first "
               "NMax events."
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  if (argc == 1) {
    SayUsage(argv);
    exit(1);
  }
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::inputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-n") {
      cliopts::NMax = string_parsers::str2T<size_t>(argv[++opt]);
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

  response_helper nrh(cliopts::fclname);
  std::cout << "[INFO]: Loaded parameters: " << std::endl
            << nrh.GetHeaderInfo() << std::endl;

  TFile *f = TFile::Open(cliopts::inputfile.c_str());
  if (!f || !f->IsOpen()) {
    std::cout << "[ERROR]: Failed to open " << cliopts::inputfile
              << " for reading." << std::endl;
    return 2;
  }
  TTree *gevs = dynamic_cast<TTree *>(f->Get("gtree"));
  if (!gevs) {
    std::cout << "[ERROR]: Failed to read TTree, " << std::quoted("gtree")
              << ", from " << argv[2] << "." << std::endl;
    return 3;
  }

  genie::NtpMCEventRecord *GenieNtpl = nullptr;

  if (gevs->SetBranchAddress("gmcrec", &GenieNtpl) != TTree::kMatch) {
    std::cout << "[ERROR]: Failed to set branch address on ghep tree."
              << std::endl;
    return 4;
  }

  size_t NEvs = gevs->GetEntries();

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  TFile *of = new TFile(cliopts::outputfile.c_str(), "RECREATE");
  TTree *ot = new TTree("resp_tree", "");

  param_list_t params = nrh.GetParameters();
  std::unique_ptr<PrecalculatedResponseReader<6>> prr =
      PrecalculatedResponseReader<6>::MakeTreeWriter(nrh.GetHeaders(), ot);

  size_t NToRead = std::min(NEvs, cliopts::NMax);
  size_t NToShout = NToRead / 100;
  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    if (NToShout && !(ev_it % NToShout)) {
      std::cout << "Event #" << ev_it
                << ", Interaction: " << GenieNtpl->event->Summary()->AsString()
                << std::endl;
    }

    prr->AddEventResponses(nrh.GetEventResponses(*GenieNtpl->event));
  }
  of->Write();
  of->Close();
}
