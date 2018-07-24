#include "nusystematics/artless/response_helper.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"

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

int main(int argc, char const *argv[]) {

  if ((argc == 2) &&
      ((std::string("--help") == argv[1]) || (std::string("-?") == argv[1]))) {
    std::cout << "[USAGE]: " << argv[0]
              << " <syst configuration fhicl> <GHEP file>" << std::endl;
    return 0;
  }

  if (argc < 3) {
    std::cout << "[ERROR]: Expected to be passed 2 parameters." << std::endl;
    std::cout << "[USAGE]: " << argv[0]
              << " <syst configuration fhicl> <GHEP file>" << std::endl;
    return 1;
  }
  size_t NMax = std::numeric_limits<size_t>::max();
  if (argc >= 4) {
    NMax = string_parsers::str2T<size_t>(argv[3]);
  }

  response_helper nrh(argv[1]);
  std::cout << "[INFO]: Loaded parameters: " << std::endl
            << nrh.GetHeaderInfo() << std::endl;

  TFile *f = TFile::Open(argv[2]);
  if (!f || !f->IsOpen()) {
    std::cout << "[ERROR]: Failed to open " << argv[2] << " for reading."
              << std::endl;
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

  param_list_t params = nrh.GetParameters();
  size_t NToRead = std::min(NEvs, NMax);
  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord *GenieGHep = GenieNtpl->event;
    std::cout << "Event #" << ev_it
              << ", Interaction: " << GenieGHep->Summary()->AsString()
              << std::endl;

    event_unit_response_t resp = nrh.GetEventResponses(*GenieGHep);
    std::cout << "[INFO]: Response =  " << std::endl
              << nrh.GetEventResponseInfo(resp) << std::endl;
  }
}
