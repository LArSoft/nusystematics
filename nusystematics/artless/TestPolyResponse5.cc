#include "nusystematics/artless/response_helper.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"

// Included by fhiclcpp-simple will not be available in art-ful
#include "string_parsers/from_string.hxx"
#include "string_parsers/to_string.hxx"

#include "TFile.h"
#include "TTree.h"

#include <chrono>
#include <iostream>
#include <random>

using namespace fhicl;
using namespace systtools;
using namespace nusyst;

namespace cliopts {
std::string fclname = "";
std::string inputfile = "";
std::string outputfile = "";
size_t NMax = std::numeric_limits<size_t>::max();

} // namespace cliopts

constexpr size_t NTests = 10;

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
  response_helper nrh_other(cliopts::fclname);
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
  TTree *ot = new TTree("resp_valid_tree", "");
  size_t NHeaders = nrh.GetHeaders().size();
  Int_t NIds;
  Int_t *id = new Int_t[NHeaders];
  Double_t *vals = new Double_t[NHeaders * NTests];
  Double_t *responses4 = new Double_t[NHeaders * NTests];
  Double_t *responses5 = new Double_t[NHeaders * NTests];
  Double_t *responses6 = new Double_t[NHeaders * NTests];
  Double_t *spline3 = new Double_t[NHeaders * NTests];
  Double_t *calced = new Double_t[NHeaders * NTests];

  ot->Branch("nids", &NIds, "nids/I");
  ot->Branch("id", id, "id[nids]/I");
  ot->Branch("vals", vals,
             (std::string("vals[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());
  ot->Branch("responses4", responses4,
             (std::string("responses4[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());
  ot->Branch("responses5", responses5,
             (std::string("responses5[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());
  ot->Branch("responses6", responses6,
             (std::string("responses6[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());
  ot->Branch("spline3", spline3,
             (std::string("spline3[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());
  ot->Branch("calced", calced,
             (std::string("calced[nids][") +
              string_parsers::T2Str<size_t>(NTests) + "]/D")
                 .c_str());

  param_list_t params = nrh.GetParameters();
  size_t NToRead = std::min(NEvs, cliopts::NMax);
  size_t NToShout = NToRead / 100;

  std::mt19937_64 generator(
      std::chrono::steady_clock::now().time_since_epoch().count());
  std::normal_distribution<double> distribution(0, 1.5);
  std::function<double()> RNJesus = std::bind(distribution, generator);

  size_t NDumped_005 = 0;
  size_t NDumped_05 = 0;
  size_t NDumped_95 = 0;

  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    if (NToShout && !(ev_it % NToShout)) {
      std::cout << "Event #" << ev_it
                << ", Interaction: " << GenieNtpl->event->Summary()->AsString()
                << std::endl;
    }

    event_unit_response_t resp = nrh.GetEventResponses(*GenieNtpl->event);
    ScrubUnityEventResponses(resp);

    NIds = 0;
    for (auto const &pr : resp) {
      id[NIds] = pr.pid;

      PolyResponse<4> const &poly4 = nrh.GetPolyResponse<4>(pr.pid, resp);
      PolyResponse<5> const &poly5 = nrh.GetPolyResponse<5>(pr.pid, resp);
      PolyResponse<6> const &poly6 = nrh.GetPolyResponse<6>(pr.pid, resp);
      TSpline3 sp = nrh.GetSpline(pr.pid, resp);

      for (size_t i = 0; i < NTests; ++i) {
        vals[NIds * NTests + i] = RNJesus();
        while (fabs(vals[NIds * NTests + i]) > 2) {
          vals[NIds * NTests + i] = RNJesus();
        }
        responses4[NIds * NTests + i] = poly4.eval(vals[NIds * NTests + i]);
        responses5[NIds * NTests + i] = poly5.eval(vals[NIds * NTests + i]);
        responses6[NIds * NTests + i] = poly6.eval(vals[NIds * NTests + i]);
        spline3[NIds * NTests + i] = sp.Eval(vals[NIds * NTests + i]);
        calced[NIds * NTests + i] = nrh_other.GetEventWeightResponse(
            *GenieNtpl->event, {{pr.pid, vals[NIds * NTests + i]}});

        size_t NDumpSteps = 20;
        double poly5_diff =
            fabs((calced[NIds * NTests + i] - responses5[NIds * NTests + i]) /
                 calced[NIds * NTests + i]);

        if (!i && (poly5_diff > 0.05) && (poly5_diff < 0.1) &&
            (NDumped_005 < 10)) {

          TGraph g_thrown, g_calced, g_responses5;
          g_thrown.SetPoint(0, vals[NIds * NTests + i],
                            calced[NIds * NTests + i]);
          std::cout << "At " << vals[NIds * NTests + i] << " = "
                    << calced[NIds * NTests + i] << " | "
                    << responses5[NIds * NTests + i] << std::endl;

          SystParamHeader const &hdr = nrh.GetHeader(pr.pid);

          for (size_t step = 0; step < hdr.paramVariations.size(); ++step) {
            g_calced.SetPoint(step, hdr.paramVariations[step],
                              pr.responses[step]);
          }

          double step_w =
              (hdr.paramVariations.back() - hdr.paramVariations.front()) /
              double(NDumpSteps);

          for (size_t step = 0; step < NDumpSteps; ++step) {
            g_responses5.SetPoint(
                step, hdr.paramVariations.front() + step * step_w,
                poly5.eval(hdr.paramVariations.front() + step * step_w));
          }

          g_calced.Write(
              (std::string("g_calced_005_") + std::to_string(NDumped_005))
                  .c_str());
          g_responses5.Write(
              (std::string("g_responses5_005_") + std::to_string(NDumped_005))
                  .c_str());
          g_thrown.Write(
              (std::string("g_thrown_005_") + std::to_string(NDumped_005))
                  .c_str());

          NDumped_005++;
        }

        if (!i && (poly5_diff > 0.5) && (poly5_diff < 0.8) &&
            (NDumped_05 < 10)) {

          TGraph g_thrown, g_calced, g_responses5;
          g_thrown.SetPoint(0, vals[NIds * NTests + i],
                            calced[NIds * NTests + i]);
          std::cout << "At " << vals[NIds * NTests + i] << " = "
                    << calced[NIds * NTests + i] << " | "
                    << responses5[NIds * NTests + i] << std::endl;

          SystParamHeader const &hdr = nrh.GetHeader(pr.pid);

          for (size_t step = 0; step < hdr.paramVariations.size(); ++step) {
            g_calced.SetPoint(step, hdr.paramVariations[step],
                              pr.responses[step]);
          }

          double step_w =
              (hdr.paramVariations.back() - hdr.paramVariations.front()) /
              double(NDumpSteps);

          for (size_t step = 0; step < NDumpSteps; ++step) {
            g_responses5.SetPoint(
                step, hdr.paramVariations.front() + step * step_w,
                poly5.eval(hdr.paramVariations.front() + step * step_w));
          }

          g_calced.Write(
              (std::string("g_calced_05_") + std::to_string(NDumped_05))
                  .c_str());
          g_responses5.Write(
              (std::string("g_responses5_05_") + std::to_string(NDumped_05))
                  .c_str());
          g_thrown.Write(
              (std::string("g_thrown_05_") + std::to_string(NDumped_05))
                  .c_str());

          NDumped_05++;
        }

        if (!i && (poly5_diff > 0.5) && (poly5_diff < 0.8) &&
            (NDumped_95 < 10)) {

          TGraph g_thrown, g_calced, g_responses5;
          g_thrown.SetPoint(0, vals[NIds * NTests + i],
                            calced[NIds * NTests + i]);
          std::cout << "At " << vals[NIds * NTests + i] << " = "
                    << calced[NIds * NTests + i] << " | "
                    << responses5[NIds * NTests + i] << std::endl;
                    
          SystParamHeader const &hdr = nrh.GetHeader(pr.pid);

          for (size_t step = 0; step < hdr.paramVariations.size(); ++step) {
            g_calced.SetPoint(step, hdr.paramVariations[step],
                              pr.responses[step]);
          }

          double step_w =
              (hdr.paramVariations.back() - hdr.paramVariations.front()) /
              double(NDumpSteps);

          for (size_t step = 0; step < NDumpSteps; ++step) {
            g_responses5.SetPoint(
                step, hdr.paramVariations.front() + step * step_w,
                poly5.eval(hdr.paramVariations.front() + step * step_w));
          }

          g_calced.Write(
              (std::string("g_calced_95_") + std::to_string(NDumped_95))
                  .c_str());
          g_responses5.Write(
              (std::string("g_responses5_95_") + std::to_string(NDumped_95))
                  .c_str());
          g_thrown.Write(
              (std::string("g_thrown_95_") + std::to_string(NDumped_95))
                  .c_str());

          NDumped_95++;
        }
      }
      NIds++;
    }
    ot->Fill();
  }
  of->Write();
  of->Close();
}
