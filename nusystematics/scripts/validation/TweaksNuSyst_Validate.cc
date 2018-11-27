#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <vector>

// Validation script to run on output from DumpConfiguredTweaksNuSyst
// Makes a multi-page pdf of mode contributions, variations and bad parameters
// (nan on neg weights)
//
// Run by root -l -b -q 'TweaksNuSyst_Validate.cc("YOUR_OUTPUT_FILE_HERE.root")'

std::string GoodHistogram(TH1D *Hist) {
  // First check integral
  if (TMath::IsNaN(Hist->Integral()))
    return std::string("isnan");

  // Check for negative content
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) {
    double content = Hist->GetBinContent(i + 1);
    // Check negative content
    if (content < 0)
      return std::string("negative");
  }

  return std::string("");
}

void TweaksNuSyst_Validate(std::string filename) {

  std::cout << "[INFO]: Reading file: " << filename << std::endl;

  TFile *file = new TFile(filename.c_str(), "OPEN");
  TTree *tree = (TTree *)file->Get("events")->Clone();

  // Read through the file
  //
  // Find the tweak responses (tweak_response_ParName)
  // The number of tweaks (ntweaks_ParName)
  // The CV weight (paramCVWeight_ParName)
  //
  // The number of tweaks _SHOULD_ be 7?

  int nbr = tree->GetListOfBranches()->GetEntries();

  // Array of names of tweaks
  std::vector<std::string> TweakArray;
  // Array of number of tweaks
  std::vector<std::string> NumberArray;
  // Array of CV of tweaks
  std::vector<std::string> CVArray;
  // Parameter names
  std::vector<std::string> ParamNames;

  // Get the names of all the branches
  for (int i = 0; i < nbr; ++i) {
    std::string name = std::string(tree->GetListOfBranches()->At(i)->GetName());
    if (name.find("tweak_responses_") != std::string::npos) {
      TweakArray.push_back(name);
      std::string paramname =
          name.substr(std::string("tweak_responses_").size(), name.size());
      ParamNames.push_back(paramname);
    } else if (name.find("ntweaks_") != std::string::npos)
      NumberArray.push_back(name);
    else if (name.find("paramCVWeight_") != std::string::npos)
      CVArray.push_back(name);
  }

  std::cout << "Found total " << TweakArray.size() << " tweaks" << std::endl;
  std::cout << "Found total " << NumberArray.size() << " numbers" << std::endl;
  size_t NumberOfPoints = 7;
  for (size_t i = 0; i < NumberArray.size(); ++i) {
    int max = tree->GetMaximum(NumberArray[i].c_str());
    int min = tree->GetMinimum(NumberArray[i].c_str());
    if (max > int(NumberOfPoints)) {
      std::cerr << "Eeek, minimum is not maximum for parameter " << i << " = "
                << NumberArray[i] << std::endl;
      std::cerr << "  Minimum: " << min << ", Maximum: " << max << std::endl;
      break;
    }
  }

  std::cout << "Found total " << CVArray.size() << " CVs" << std::endl;

  int nEntries = tree->GetEntries();
  (void)nEntries;

  // The 1D distributions to draw
  std::vector<std::string> DrawStr;
  DrawStr.push_back("e_nu_GeV");
  DrawStr.push_back("Q2_GeV2");
  DrawStr.push_back("W_GeV2");
  DrawStr.push_back("q0_GeV");
  DrawStr.push_back("q3_GeV");

  // Binning for above distributions
  std::vector<std::string> BinningStr;
  BinningStr.push_back("40, 0, 5");
  BinningStr.push_back("40, 0, 2");
  BinningStr.push_back("40, 0.9, 3.5");
  BinningStr.push_back("40, 0, 2");
  BinningStr.push_back("40, 0, 2");

  // The interaction modes
  std::vector<std::string> IntStr;
  IntStr.push_back("is_qe");
  IntStr.push_back("is_mec");
  IntStr.push_back("is_res");
  IntStr.push_back("is_dis");

  // The CC NC flag
  std::vector<std::string> CCStr;
  CCStr.push_back("is_cc");
  CCStr.push_back("!is_cc");

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetTopMargin(canv->GetTopMargin() * 0.8);
  canv->SetRightMargin(canv->GetRightMargin() * 0.7);
  canv->SetBottomMargin(canv->GetBottomMargin() * 0.8);
  canv->SetLeftMargin(canv->GetLeftMargin() * 1.2);
  std::string canvname =
      filename.substr(0, filename.find(".root")) + "_VALIDATION";
  canv->Print((canvname + ".pdf[").c_str());

  // Print the parameters
  // Make a TPaveText at the end
  TPaveText *pttitle = new TPaveText(0.01, 0.90, 0.99, 0.99, "NDC");
  pttitle->AddText("Parameters");
  pttitle->Draw();
  TPaveText *pt = new TPaveText(0.01, 0.01, 0.495, 0.89, "NDC");
  pt->Draw("same");
  for (size_t i = 0; i < ParamNames.size() / 2 + 1; ++i) {
    pt->AddText(ParamNames[i].c_str());
  }
  pt->SetTextAlign(11);
  TPaveText *pt2 = new TPaveText(0.505, 0.01, 0.99, 0.89, "NDC");
  pt2->Draw("same");
  for (size_t i = ParamNames.size() / 2 + 1; i < ParamNames.size(); ++i) {
    pt2->AddText(ParamNames[i].c_str());
  }
  pt2->SetTextAlign(11);
  canv->Print((canvname + ".pdf").c_str());
  delete pttitle;
  delete pt;
  delete pt2;

  // First draw some silly 1D distributions
  for (size_t i = 0; i < CCStr.size() - 1; ++i) {
    canv->Clear();
    tree->Draw(Form("%s>>CCStr_%ld", CCStr[i].c_str(), i));
    TH1D *temp = (TH1D *)gDirectory->Get(Form("CCStr_%ld", i));
    temp->Scale(1 / temp->Integral());
    temp->GetYaxis()->SetTitle("Fractional contribution");
    temp->Draw();
    temp->SetLineColor(kRed);
    temp->SetLineWidth(2);
    canv->Print((canvname + ".pdf").c_str());
  }

  for (size_t i = 0; i < IntStr.size(); ++i) {
    canv->Clear();
    tree->Draw(Form("%s>>IntStr_%ld", IntStr[i].c_str(), i));
    TH1D *temp = (TH1D *)gDirectory->Get(Form("IntStr_%ld", i));
    temp->Scale(1 / temp->Integral());
    temp->GetYaxis()->SetTitle("Fractional contribution");
    temp->SetLineColor(kRed);
    temp->SetLineWidth(2);
    temp->Draw();
    canv->Print((canvname + ".pdf").c_str());
  }

  std::vector<double> Integrals;
  std::vector<std::string> IntegralNames;
  for (size_t i = 0; i < CCStr.size(); ++i) {
    for (size_t j = 0; j < IntStr.size(); ++j) {
      canv->Clear();
      tree->Draw(Form("1>>CCStr_%ld_IntStr_%ld", i, j),
                 Form("%s*%s", CCStr[i].c_str(), IntStr[j].c_str()));
      TH1D *temp = (TH1D *)gDirectory->Get(Form("CCStr_%ld_IntStr_%ld", i, j));
      Integrals.push_back(temp->Integral());
      IntegralNames.push_back(CCStr[i] + " " + IntStr[j]);
    }
  }

  size_t ModeSize = IntStr.size() * CCStr.size();
  TH1D *ModeValidation =
      new TH1D("ModeValidation", "ModeValidation", ModeSize, 0, ModeSize);
  for (size_t i = 0; i < ModeSize; ++i) {
    ModeValidation->SetBinContent(i + 1, Integrals[i]);
    ModeValidation->GetXaxis()->SetBinLabel(i + 1, IntegralNames[i].c_str());
  }
  ModeValidation->Scale(1 / ModeValidation->Integral());
  ModeValidation->GetYaxis()->SetTitle("Fractional contribution");
  TLine *ModeLine =
      new TLine(ModeValidation->GetBinLowEdge(IntStr.size() + 1), 0,
                ModeValidation->GetBinLowEdge(IntStr.size() + 1),
                ModeValidation->GetMaximum());
  ModeLine->SetLineWidth(4);
  ModeLine->SetLineColor(kRed);
  canv->Clear();
  ModeValidation->Draw();
  ModeLine->Draw("same");
  canv->Print((canvname + ".pdf").c_str());
  delete ModeLine;
  delete ModeValidation;

  // The histograms containing the -3,-2,-1... variations
  TH1D **PlotList = new TH1D *[NumberOfPoints];
  // Premade legend
  TLegend *leg = new TLegend(0.13, 0.5, 0.5, 0.92);
  for (size_t z = 0; z < NumberOfPoints; ++z) {
    PlotList[z] = new TH1D(Form("List_%ld", z), Form("List_%ld", z), 1, 0, 1);
    PlotList[z]->SetLineWidth(2);
    int color = 0;
    std::string name = "";
    switch (z) {
    case 0:
      color = kRed - 4;
      name = "-3#sigma";
      break;
    case 1:
      color = kYellow - 3;
      name = "-2#sigma";
      break;
      break;
    case 2:
      color = kBlue - 4;
      name = "-1#sigma";
      break;
    case 3:
      color = kGray + 2;
      name = "Central";
      break;
    case 4:
      color = kBlue - 4;
      name = "+1#sigma";
      break;
    case 5:
      color = kYellow - 3;
      name = "+2#sigma";
      break;
    case 6:
      color = kRed - 4;
      name = "+3#sigma";
      break;
    }
    PlotList[z]->SetTitle(name.c_str());
    PlotList[z]->SetLineColor(color);
    if (z > 3)
      PlotList[z]->SetLineStyle(kDashed);
  }
  for (size_t z = 0; z < NumberOfPoints; ++z) {
    leg->AddEntry(PlotList[NumberOfPoints - z - 1], "", "l");
  }

  // Have a vector with no effect and bad effect parameters and selections
  std::vector<std::string> NoEffect;
  std::vector<std::string> BadEffect;

  // Output the found parameters
  std::cout << "Parameters: " << std::endl;
  for (size_t i = 0; i < ParamNames.size(); ++i)
    std::cout << "   " << ParamNames[i] << std::endl;

  int ndraws = 0;
  // Loop over the central values
  for (size_t a = 0; a < CVArray.size(); ++a) {

    // Loop over CC and not CC (i.e. NC)
    for (size_t ab = 0; ab < CCStr.size(); ++ab) {

      // Loop over the interaction string selection
      for (size_t b = 0; b < IntStr.size(); ++b) {
        std::vector<std::string> SelStr;
        SelStr.push_back(CVArray[a] + "*" + IntStr[b] + "*" + CCStr[ab]);

        // Loop over the variables we want to draw
        for (size_t i = 0; i < DrawStr.size(); ++i) {

          // Loop over the selections
          for (size_t j = 0; j < SelStr.size(); ++j) {
            // Get the reference distribution (no tune)
            tree->Draw(Form("%s>>ref_%ld_%ld_%ld_%ld(%s)", DrawStr[i].c_str(), a, b,
                            ab, i, BinningStr[i].c_str()),
                       std::string(IntStr[b] + "*" + CCStr[ab]).c_str());
            TH1D *ref =
                (TH1D *)gDirectory->Get(Form("ref_%ld_%ld_%ld_%ld", a, b, ab, i))
                    ->Clone();
            ref->GetXaxis()->SetTitle(DrawStr[i].c_str());

            // Get the maximum
            ref->SetTitle(Form("REF %s %s %ld", ParamNames[a].c_str(),
                               std::string(IntStr[b] + "*" + CCStr[ab]).c_str(),
                               j));
            // Skip empty histograms
            if (ref->Integral() == 0) {
              delete ref;
              ref = NULL;
              continue;
            }

            // Check bad reference histograms (no weighting)
            if (!GoodHistogram(ref).empty()) {
              std::cout << DrawStr[i] << " not good: " << GoodHistogram(ref)
                        << std::endl;
              break;
            }

            bool draw = true;
            // Loop over the -3, -2, -1, 0... variations
            for (size_t z = 0; z < NumberOfPoints; ++z) {
              PlotList[z] = NULL;
              if (!draw)
                break;

              std::string rawit = std::string(Form(
                  "%s*%s[%ld]", SelStr[j].c_str(), TweakArray[a].c_str(), z));

              // Make the draw string
              tree->Draw(Form("%s>>hist_%ld_%ld_%ld_%ld_%ld_%ld(%s)",
                              DrawStr[i].c_str(), i, j, a, b, ab, z,
                              BinningStr[i].c_str()),
                         rawit.c_str());
              TH1D *temp =
                  (TH1D *)gDirectory
                      ->Get(Form("hist_%ld_%ld_%ld_%ld_%ld_%ld", i, j, a, b, ab, z))
                      ->Clone();
              temp->GetXaxis()->SetTitle(ref->GetXaxis()->GetTitle());
              temp->SetTitle(Form("%s %s*%s", ParamNames[a].c_str(),
                                  CCStr[ab].c_str(), IntStr[b].c_str()));
              PlotList[z] = temp;
              if (temp->Integral() == 0) {
                draw = false;
                NoEffect.push_back(rawit);
                delete PlotList[z];
                PlotList[z] = NULL;
                continue;
              }
              draw = true;

              // Check the histograms is good
              if (!GoodHistogram(temp).empty()) {
                BadEffect.push_back(GoodHistogram(temp) + "_" + SelStr[j]);
                if (GoodHistogram(temp) == "isnan") {
                  draw = false;
                  delete ref;
                  delete PlotList[z];
                  PlotList[z] = NULL;
                  break;
                }
              }
              // std::cout << DrawStr[i] << " with " << rawit << ": " <<
              // temp->Integral() << std::endl;
            }

            if (draw) {
              canv->Clear();
              ref->SetLineWidth(2);
              ref->SetLineColor(kBlack);
              TLine *line = new TLine(ref->GetXaxis()->GetBinLowEdge(1), 1,
                                      ref->GetXaxis()->GetBinLowEdge(
                                          ref->GetXaxis()->GetNbins() + 1),
                                      1);
              line->SetLineWidth(2);
              line->SetLineStyle(kDashed);
              line->SetLineColor(kBlack);
              double maximum = 0.0;
              double minimum = 10;
              for (size_t z = 0; z < NumberOfPoints; ++z) {
                PlotList[z]->SetLineWidth(2);
                int color = 0;
                switch (z) {
                case 0:
                  color = kRed - 4;
                  break;
                case 1:
                  color = kYellow - 3;
                  break;
                case 2:
                  color = kBlue - 4;
                  break;
                case 3:
                  color = kGray + 2;
                  break;
                case 4:
                  color = kBlue - 4;
                  break;
                case 5:
                  color = kYellow - 3;
                  break;
                case 6:
                  color = kRed - 4;
                  break;
                }
                if (z > 3)
                  PlotList[z]->SetLineStyle(kDashed);
                PlotList[z]->GetYaxis()->SetTitle("Tune/Nominal Value");
                PlotList[z]->SetLineColor(color);
                PlotList[z]->Divide(ref);
                if (maximum < PlotList[z]->GetMaximum())
                  maximum = PlotList[z]->GetMaximum() * 1.2;
                if (minimum > PlotList[z]->GetMinimum())
                  minimum = PlotList[z]->GetMinimum() * 0.8;
              }

              // We now have all of the plots in PlotList
              // Check to see if they're just all the same
              if (double(PlotList[0]->Integral()) ==
                  double(PlotList[3]->Integral())) {
                NoEffect.push_back(SelStr[j]);
                delete ref;
                delete line;
                for (size_t z = 0; z < NumberOfPoints; ++z) {
                  delete PlotList[z];
                  PlotList[z] = NULL;
                }
                continue;
              }

              // Now draw
              for (size_t z = 0; z < NumberOfPoints; ++z) {
                if (z == 0) {
                  PlotList[z]->Draw();
                  PlotList[z]->GetYaxis()->SetRangeUser(minimum, maximum);
                } else {
                  PlotList[z]->Draw("same");
                }
              }
              line->Draw("same");
              leg->Draw("same");
              canv->Print((canvname + ".pdf").c_str());
              ndraws++;

              delete ref;
              delete line;
              for (size_t z = 0; z < NumberOfPoints; ++z) {
                delete PlotList[z];
                PlotList[z] = NULL;
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Draw " << ndraws << " parameters" << std::endl;
  std::cout << "NoEffect " << NoEffect.size() << " combinations" << std::endl;
  std::cout << "Bad: " << std::endl;
  std::vector<std::string> FancyFormattedBad;

  if (!BadEffect.empty()) {
    std::string first = BadEffect[0];
    std::cout << "    " << BadEffect[0] << std::endl;
    FancyFormattedBad.push_back(BadEffect[0]);
    for (size_t i = 0; i < BadEffect.size(); ++i) {
      if (BadEffect[i] != first) {
        std::cout << "    " << BadEffect[i] << std::endl;
        first = BadEffect[i];
        FancyFormattedBad.push_back(first);
      }
    }
  }

  // Print the parameters
  // Make a TPaveText at the end
  canv->Clear();
  TPaveText *pttitle_2 = new TPaveText(0.01, 0.90, 0.99, 0.99, "NDC");
  pttitle_2->AddText("Parameters");
  pttitle_2->AddText("(bad)");
  ((TText *)pttitle_2->GetListOfLines()->Last())->SetTextColor(kRed);
  pttitle_2->Draw();
  TPaveText *pt_2 = new TPaveText(0.01, 0.01, 0.495, 0.89, "NDC");
  pt_2->Draw("same");
  for (size_t i = 0; i < ParamNames.size() / 2 + 1; ++i) {
    pt_2->AddText(ParamNames[i].c_str());
    // Loop over the bad parameters
    for (size_t j = 0; j < FancyFormattedBad.size(); ++j) {
      if (FancyFormattedBad[j].find(ParamNames[i]) != std::string::npos) {
        ((TText *)pt_2->GetListOfLines()->Last())->SetTextColor(kRed);
      }
    }
  }
  pt_2->SetTextAlign(11);
  TPaveText *pt2_2 = new TPaveText(0.505, 0.01, 0.99, 0.89, "NDC");
  pt2_2->Draw("same");
  for (size_t i = ParamNames.size() / 2 + 1; i < ParamNames.size(); ++i) {
    pt2_2->AddText(ParamNames[i].c_str());
    for (size_t j = 0; j < FancyFormattedBad.size(); ++j) {
      if (FancyFormattedBad[j].find(ParamNames[i]) != std::string::npos) {
        ((TText *)pt2_2->GetListOfLines()->Last())->SetTextColor(kRed);
      }
    }
  }
  pt2_2->SetTextAlign(11);
  canv->Print((canvname + ".pdf").c_str());
  delete pttitle_2;
  delete pt_2;
  delete pt2_2;

  // Make a TPaveText at the end
  canv->Clear();
  pttitle = new TPaveText(0.01, 0.90, 0.99, 0.99, "NDC");
  pttitle->AddText("Bad parameters (nan or neg weights)");
  pttitle->Draw();
  pt = new TPaveText(0.01, 0.01, 0.99, 0.89, "NDC");
  pt->Draw("same");
  for (size_t i = 0; i < FancyFormattedBad.size(); ++i) {
    pt->AddText(FancyFormattedBad[i].c_str());
  }
  pt->SetTextAlign(11);
  canv->Print((canvname + ".pdf").c_str());
  delete pttitle;
  delete pt;

  // Loop over all the parameters
  for (size_t i = 0; i < ParamNames.size(); ++i) {
    // Make a TPaveText at the end
    canv->Clear();
    pttitle = new TPaveText(0.01, 0.90, 0.99, 0.99, "NDC");
    pttitle->AddText(Form("No Effect %s", ParamNames[i].c_str()));
    pttitle->Draw();
    pt = new TPaveText(0.01, 0.01, 0.99, 0.89, "NDC");
    pt->Draw("same");
    std::string oldstring = "";
    for (size_t j = 0; j < NoEffect.size(); ++j) {
      if (NoEffect[j].find(ParamNames[i]) != std::string::npos) {
        // Make a "pretty string"
        std::string tempstring = NoEffect[j];
        // Cut out the parameter name
        tempstring = tempstring.substr(tempstring.find(ParamNames[i]) +
                                           ParamNames[i].size(),
                                       tempstring.size());
        // Cut out the "tweak_responses" part
        tempstring = tempstring.substr(0, tempstring.find(TweakArray[i]) - 1);
        // Finally drop the "*" character
        while (tempstring.find("*") != std::string::npos) {
          tempstring.replace(tempstring.find("*"), 1, std::string(" "));
        }
        if (oldstring != tempstring)
          pt->AddText(tempstring.c_str());
        oldstring = tempstring;
      }
    }
    pt->SetTextAlign(11);
    // pt->SetTextSize(6);
    canv->Print((canvname + ".pdf").c_str());
    delete pttitle;
    delete pt;
  }

  canv->Print((canvname + ".pdf]").c_str());
}

#ifndef __CINT__
int main(int, char const *argv[]) { TweaksNuSyst_Validate(argv[1]); }
#endif
