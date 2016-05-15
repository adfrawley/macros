#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TRandom1.h>
#include <TPolyLine.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
//#include <iomanip.h>
#include <TLorentzVector.h>

//#define TEST

void look_quarkonia_histos()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);

  //TFile *fin = new TFile("ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("maps_7layer_hijing_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("maps_7layer_single_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("maps_5layer_single_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("pixels_strips_7layers_embed_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("pixel_maps_mapsouter_single_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("pixel_maps_mapsouter_embed_ups1s_qual3.00_dca2d0.10.root");
  //TFile *fin = new TFile("maps_pixel_mapsouter_single_ups1s_qual3.00_dca2d0.10.root");
  TFile *fin = new TFile("maps_pixel_mapsouter_embed_ups1s_qual3.00_dca2d0.10.root");

  if(!fin)
    {
      cout << " failed to open input file " << endl;
      exit(1);
    }


  TH1F *hrquality;
  fin->GetObject("hrquality",hrquality);
  TH1F *hrdca2d;
  fin->GetObject("hrdca2d",hrdca2d);
  if(!hrquality)
    {
      cout << "Did not find hrquality" << endl;
      exit(1);
    }
  if(!hrdca2d)
    {
      cout << "Did not find hrdca2d" << endl;
     exit(1);
    }

  TCanvas *cq = new TCanvas("cq","cq",5,5,600,600 );
  cq->Divide(1,2);
  cq->cd(1);
  hrquality->Draw();
  cq->cd(2);
  gPad->SetLogy(1);
  hrdca2d->Draw();

  // Mass histos

  TH1D *recomass_primary;
  fin->GetObject("recomass_primary",recomass_primary);
  TH1D *recomass;
  fin->GetObject("recomass",recomass);
  TH1D *g4mass_primary;
  fin->GetObject("g4mass_primary",g4mass_primary);
  TH1D *g4mass;
  fin->GetObject("g4mass",g4mass);
  
  TCanvas *cmass = new TCanvas("cmass","cmass",10,10,800,600);
  gPad->SetLogy(1);

  recomass_primary->SetLineColor(kRed);
  recomass_primary->DrawCopy();  

  recomass->SetLineColor(kBlack);
  recomass->DrawCopy("same");  


  TCanvas *cm_comp = new TCanvas("cm_comp","cm_comp",10,10,800,800);
  cm_comp->Divide(2,1);
  cm_comp->cd(1);
  gPad->SetLogy(1);
  recomass_primary->Draw();
  recomass->Draw("same");
  // we want from 7 to 11 GeV/c^2 - the whole range
  double yreco = recomass->Integral();
 double yreco_primary = recomass_primary->Integral();
 cout << "Reconstructed mass spectrum has " << yreco_primary << " entries from primary tracks and " << yreco << " entries total " << endl;

  cm_comp->cd(2);
  gPad->SetLogy(1);
  
  g4mass_primary->SetLineColor(kRed);
  g4mass_primary->Draw();
  g4mass->Draw("same");
  double yg4_primary = g4mass_primary->Integral();
  double yg4 = g4mass->Integral();
  cout << "G4 mass spectrum has " << yg4_primary << " entries from primary tracks and  " << yg4 << " entries total" << endl;

  cout << "Reconstruction efficiency is " << yreco_primary/yg4_primary << endl;

}
