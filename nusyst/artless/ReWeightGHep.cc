#include "configure_nusyst_providers.hh"

#include "nusyst/systproviders/GENIEReWeight_tool.hh"

#include "larsyst/interpreters/ParamHeaderHelper.hh"
#include "larsyst/interpreters/load_parameter_headers.hh"

#include "fhiclcpp/make_ParameterSet.h"
#include "string_parsers/to_string.hxx"

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCEventRecord.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace fhicl;
using namespace larsyst;

int main(int argc, char const *argv[]) {

  if ((argc == 2) &&
      ((std::string("--help") == argv[1]) || (std::string("-?") == argv[1]))) {
    std::cout << "[USAGE]: " << argv[0]
              << " <syst configuration fhicl> <GHEP file>" << std::endl;
    return 0;
  }

  if (argc != 3) {
    std::cout << "[ERROR]: Expected to be passed 2 parameters." << std::endl;
    std::cout << "[USAGE]: " << argv[0]
              << " <syst configuration fhicl> <GHEP file>" << std::endl;
    return 1;
  }

  ParameterSet ps = make_ParameterSet(argv[1]).get<ParameterSet>(
      "generated_systematic_provider_configuration");

  provider_map_t syst_providers = load_syst_provider_configuration(ps);

  param_header_map_t configuredParameterHeaders =
      load_syst_provider_headers(ps);
  if (!configuredParameterHeaders.size()) {
    std::cout << "[ERROR]: Expected to find some parameter headers in "
              << std::quoted(argv[1]) << std::endl;
    throw;
  }

  ParamHeaderHelper phh(configuredParameterHeaders);

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

  for (size_t ev_it = 0; ev_it < NEvs; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord *GenieGHep = GenieNtpl->event;
    std::cout << "Event " << ev_it << std::endl;

    for (auto &sp : syst_providers) {
      GENIEReWeight *grw = dynamic_cast<GENIEReWeight *>(sp.second.get());
      if (!grw) {
        std::cout
            << "[ERROR]: Only know how to handle GENIEReWeight tools, not "
            << std::quoted(sp.second->GetToolType()) << std::endl;
        return 5;
      }
      larsyst::event_unit_response_t response =
          grw->GetEventResponse(*GenieGHep);
      std::cout << "  Provider: " << grw->GetToolType() << ":"
                << grw->GetInstanceName() << std::endl;
      for (auto pr_pair : response) {
        std::cout << "    Parameter: " << pr_pair.first << " = "
                  << std::quoted(phh.GetHeader(pr_pair.first).prettyName)
                  << " response = "
                  << string_parsers::T2Str<std::vector<double>>(pr_pair.second)
                  << std::endl;
        std::cout << "      Spline response at p = 0.345 = "
                  << phh.GetParameterResponse(pr_pair.first, 0.345,
                                              pr_pair.second)
                  << std::endl;
      }
    }
  }
}
