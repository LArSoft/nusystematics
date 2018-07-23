#include "nusystematics/artless/response_helper.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepUtils.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"

// Included by fhiclcpp-simple will not be available in art-ful
#include "string_parsers/to_string.hxx"
#include "string_parsers/from_string.hxx"

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

  bool have_MK = nrh.HaveHeader("MKSPP_Enuq0q3_response");
  bool have_MnvRPA = nrh.HaveHeader("MINERvATune_RPA");
  bool have_Mnv2p2h = nrh.HaveHeader("MINERvATune_2p2hGaussEnhancement");

  systtools::paramId_t MKRespParam =
      have_MK ? nrh.GetHeaderId("MKSPP_Enuq0q3_response")
              : systtools::kParamUnhandled<systtools::paramId_t>;
  systtools::paramId_t MnvRPARespParam =
      have_MnvRPA ? nrh.GetHeaderId("MINERvATune_RPA")
                  : systtools::kParamUnhandled<systtools::paramId_t>;
  systtools::paramId_t Mnv2p2hRespParam =
      have_Mnv2p2h ? nrh.GetHeaderId("MINERvATune_2p2hGaussEnhancement")
                   : systtools::kParamUnhandled<systtools::paramId_t>;

  TFile *fout = TFile::Open("wout.root", "RECREATE");
  TTree *MKValidTree;
  double W, Wght;
  int NEUTCh;
  if (have_MK) {
    MKValidTree = new TTree("MKValidTree", "");
    MKValidTree->SetDirectory(fout);
    MKValidTree->Branch("W", &W, "W/D");
    MKValidTree->Branch("Wght", &Wght, "Wght/D");
    MKValidTree->Branch("NEUTCh", &NEUTCh, "NEUTCh/I");
  }

  TTree *HERGValidTree;
  int event_it, param_id;
  double var, presp;
  HERGValidTree = new TTree("HERGValidTree", "");
  HERGValidTree->SetDirectory(fout);
  HERGValidTree->Branch("event_it", &event_it, "event_it/I");
  HERGValidTree->Branch("param_id", &param_id, "param_id/I");
  HERGValidTree->Branch("var", &var, "var/D");
  HERGValidTree->Branch("presp", &presp, "presp/D");

  TTree *MINERvATUNEValidTree;
  double q0, q3, Wght_RPA;
  double Wght_2p2h[4];
  int QELikeTarget;
  if (have_MnvRPA || have_Mnv2p2h) {
    MINERvATUNEValidTree = new TTree("MINERvATUNEValidTree", "");
    MINERvATUNEValidTree->SetDirectory(fout);

    MINERvATUNEValidTree->Branch("q0", &q0, "q0/D");
    MINERvATUNEValidTree->Branch("q3", &q3, "q3/D");
    if (have_MnvRPA) {
      MINERvATUNEValidTree->Branch("Wght_RPA", &Wght_RPA, "Wght_RPA/D");
    }
    if (have_Mnv2p2h) {
      MINERvATUNEValidTree->Branch("Wght_2p2h", Wght_2p2h, "Wght_2p2h[4]/D");
      MINERvATUNEValidTree->Branch("QELikeTarget", &QELikeTarget,
                                   "QELikeTarget/I");
    }
  }

  param_list_t params = nrh.GetParameters();
  size_t NToRead = std::min(NEvs, NMax);
  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord *GenieGHep = GenieNtpl->event;
    std::cout << "Event #" << ev_it
              << ", Interaction: " << GenieGHep->Summary()->AsString()
              << std::endl;

    if (!ev_it) {
      genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
          "Messenger_whisper.xml");
    }

    event_unit_response_t resp = nrh.GetEventResponses(*GenieGHep);
    std::cout << "[INFO]: Response =  " << std::endl
              << nrh.GetEventResponseInfo(resp) << std::endl;

    if (!GenieGHep->Summary()->ProcInfo().IsWeakCC()) {
      continue;
    }

    for (paramId_t p : params) {
      event_it = ev_it;
      param_id = p;
      for (size_t var_it = 0; var_it < nrh.GetNDiscreteVariations(p);
           ++var_it) {
        var = nrh.GetHeader(p).paramVariations[var_it];
        presp = nrh.GetDiscreteResponse(p, var_it, resp);
        HERGValidTree->Fill();
      }
    }

    if (GenieGHep->Summary()->ProcInfo().IsResonant() && have_MK) {
      W = GenieGHep->Summary()->Kine().W(true);
      Wght = nrh.GetDiscreteResponse(MKRespParam, 0, resp);

      NEUTCh = genie::utils::ghep::NeutReactionCode(GenieGHep);
      MKValidTree->Fill();
    }

    if ((GenieGHep->Summary()->ProcInfo().IsQuasiElastic() ||
         GenieGHep->Summary()->ProcInfo().IsMEC()) &&
        (have_MnvRPA || have_Mnv2p2h)) {

      TLorentzVector FSLepP4 = GenieGHep->Summary()->Kine().FSLeptonP4();
      TLorentzVector ISLepP4 =
          *(GenieGHep->Summary()->InitState().GetProbeP4(genie::kRfLab));
      TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

      q0 = emTransfer.E();
      q3 = emTransfer.Vect().Mag();
      if (have_MnvRPA) {
        Wght_RPA = nrh.GetDiscreteResponse(MnvRPARespParam, 0, resp);
      }
      if (have_Mnv2p2h) {
        auto const &dr = nrh.GetDiscreteResponses(Mnv2p2hRespParam, resp);
        std::cout << "[INFO]: 2p2h response size = " << dr.size() << std::endl;
        Wght_2p2h[0] = dr[0];
        Wght_2p2h[1] = dr[1];
        Wght_2p2h[2] = dr[2];
        Wght_2p2h[3] = dr[3];
      }
      QELikeTarget = e2i(GetQELikeTarget(*GenieGHep));

      MINERvATUNEValidTree->Fill();
    }
  }
  fout->Write();
  fout->Close();
}
