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

void look_quarkonia_histos_condor()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);

  TH1D *recomass=0;
  TH1D *g4mass=0;
  TH1D *recomass_primary=0;
  TH1D *g4mass_primary=0;
  TH1F *hrquality=0;
  TH1F *hrdca2d=0;

  for(int i=0;i<50;i++)
    {
      char fname[1000];
      sprintf(fname,"batch/ups1s_qual3.00_dca2d0.10_batch%i.root",i);

      TFile *fin = new TFile(fname);
      if(!fin)
	{
	  cout << " failed to open input file " << endl;
	  continue;
	}

      TH1F *hrquality_in;
      fin->GetObject("hrquality",hrquality_in);
      if(!hrquality_in)
	{
	  cout << "Did not find hrquality" << endl;
	  continue;
	}

      TH1F *hrdca2d_in;   fin->GetObject("hrdca2d",hrdca2d_in);
      TH1D *recomass_primary_in;  fin->GetObject("recomass_primary",recomass_primary_in);
      TH1D *recomass_in; fin->GetObject("recomass",recomass_in);
      TH1D *g4mass_primary_in; fin->GetObject("g4mass_primary",g4mass_primary_in);
      TH1D *g4mass_in; fin->GetObject("g4mass",g4mass_in);


      if(!recomass)
	{
	  recomass =  recomass_in;
	  recomass_primary =  recomass_primary_in;
	  g4mass = g4mass_in;
	  g4mass_primary = g4mass_primary_in;
	  hrquality = hrquality_in;
	  hrdca2d = hrdca2d_in;
	}
      else
	{
	  recomass->Add(recomass_in);
	  recomass_primary->Add(recomass_primary_in);
	  g4mass->Add(g4mass_in);
	  g4mass_primary->Add(g4mass_primary_in);
	  hrquality->Add(hrquality_in);
	  hrdca2d->Add(hrdca2d_in);
	}

    }
      
  

  TCanvas *cq = new TCanvas("cq","cq",5,5,600,600 );
  cq->Divide(1,2);
  cq->cd(1);
  hrquality->Draw();
  cq->cd(2);
  gPad->SetLogy(1);
  hrdca2d->Draw();

  // Mass histos

  
  TCanvas *cmass = new TCanvas("cmass","cmass",10,10,600,800);
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

  char fname[1000];
  sprintf(fname,"ups1s_qual3.00_dca2d0.10.root");

  cout << "Create output file " << fname << endl;

  TFile *fout1 = new TFile(fname,"recreate");
  recomass->Write();
  recomass_primary->Write();
  g4mass->Write();
  g4mass_primary->Write();
  //hrpt->Write();
  hrquality->Write();
  hrdca2d->Write();
  fout1->Close();

  cout << "Finished write to file " << fname << endl;



}
