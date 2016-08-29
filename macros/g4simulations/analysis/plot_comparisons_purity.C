#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include  <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

void plot_comparisons_purity()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  
  static const int NPLOTS = 2;
  static const int NPLOTS_UPS = 2;

  TFile *fin[NPLOTS];

  fin[0] = new TFile("root_files/maps3+tpc60_look_purity_out.root");  
  fin[1] = new TFile("root_files/maps3+intt4+tpc60_look_purity_out.root");  
  //fin[0] = new TFile("root_files/nofit_maps3+tpc60_look_purity_out.root");  
  //fin[1] = new TFile("root_files/nofit_maps3+intt4+tpc60_look_purity_out.root");  

  TGraph *gdca[NPLOTS];
  TGraph *geff[NPLOTS];
  TGraph *grdpt[NPLOTS];

  int col[3] = {kRed, kBlue, kBlack};

  TLine *lmax[2];
  lmax[0] = new TLine(0.0, 0.97, 40.0, 0.97);
  lmax[0]->SetLineColor(col[0]);
  lmax[0]->SetLineStyle(2);
  lmax[0]->SetLineWidth(3.0);
  lmax[1] = new TLine(0.0, 0.93, 40.0, 0.93);
  lmax[1]->SetLineColor(col[1]);
  lmax[1]->SetLineStyle(2);
  lmax[1]->SetLineWidth(2.0);

  for(int i=0;i<NPLOTS;i++)
    {
      if(!fin[i])
	{
	  cout << "Did not find file " << i << "  quit!" << endl;
	  exit(1);
	}
      fin[i]->GetObject("single_track_efficiency",geff[i]);
      fin[i]->GetObject("dca2d_resolution",gdca[i]);
      fin[i]->GetObject("pt_resolution",grdpt[i]);
      
      if(!geff[i])
	{
	  cout << "Failed to find geff for " << i  << endl;
	  exit(1);
	}
      if(!gdca[i])
	{
	  cout << "Failed to find gdca for " << i  << endl;
	  exit(1);
	}
      if(!grdpt[i])
	{
	  cout << "Failed to find grdpt for " << i  << endl;
	  exit(1);
	}
      
    }
  
  
   TCanvas *ceff = new TCanvas("ceff","ceff",50,50,800,600);
  ceff->SetLeftMargin(0.12);
  
  TH1F *hd = new TH1F("hd","hd",100, 0.0, 40.0);
  hd->SetMinimum(0.0);
  hd->SetMaximum(1.0);
  hd->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hd->GetXaxis()->SetTitleOffset(1.15);
  hd->GetYaxis()->SetTitle("Single track efficiency");
  hd->GetYaxis()->SetTitleOffset(1.3);
  hd->Draw();

  for(int i=0;i<NPLOTS;i++)
    {
      geff[i]->SetMarkerColor(col[i]);

      if(i==0)
	geff[i]->Draw("p");
      else
	geff[i]->Draw("same p");
    }

  TLegend *lpd = new TLegend(0.35, 0.25, 0.80, 0.40,"","NDC");
  lpd->SetBorderSize(0);
  lpd->SetFillColor(0);
  lpd->SetFillStyle(0);
  lpd->AddEntry(geff[0], "MAPS(3)+TPC(60)", "p");
  lpd->AddEntry(geff[1], "MAPS(3)+INTT(4)+TPC(60)", "p");
  lpd->Draw();
  
  for(int i=0;i<2;i++)
    lmax[i]->Draw();

  TCanvas *cdca = new TCanvas("cdca","cdca",50,50,800,600);
  cdca->SetLeftMargin(0.15);

  TH1F *hdca = new TH1F("hdca","hdca",100, 0.0, 40.0);
  hdca->SetMinimum(0.0);
  hdca->SetMaximum(0.006);
  hdca->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hdca->GetXaxis()->SetTitleOffset(1.15);
  hdca->GetYaxis()->SetTitle("dca2d resolution");
  hdca->GetYaxis()->SetTitleOffset(1.9);
  hdca->Draw();

  for(int i=0;i<NPLOTS;i++)
    {
      gdca[i]->SetMarkerColor(col[i]);

      if(i==0)
	gdca[i]->Draw("p");
      else
	gdca[i]->Draw("same p");
    }

  TLegend *lpd1 = new TLegend(0.45, 0.55, 0.89, 0.70,"","NDC");
  lpd1->SetBorderSize(0);
  lpd1->SetFillColor(0);
  lpd1->SetFillStyle(0);
  lpd1->AddEntry(geff[0], "MAPS(3)+TPC(60)", "p");
  lpd1->AddEntry(geff[1], "MAPS(3)+INTT(4)+TPC(60)", "p");
  lpd1->Draw();
  
  TCanvas *crdpt = new TCanvas("crdpt","crdpt",50,50,800,600);
  crdpt->SetLeftMargin(0.15);

  TH1F *hrdpt = new TH1F("hrdpt","hrdpt",100, 0.0, 40.0);
  hrdpt->SetMinimum(0.0);
  hrdpt->SetMaximum(0.10);
  hrdpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hrdpt->GetXaxis()->SetTitleOffset(1.15);
  hrdpt->GetYaxis()->SetTitle("#Delta p_{T} / p_{T}");
  hrdpt->GetYaxis()->SetTitleOffset(1.8);
  hrdpt->Draw();

  for(int i=0;i<NPLOTS;i++)
    {
      grdpt[i]->SetMarkerColor(col[i]);

      if(i==0)
	grdpt[i]->Draw("p");
      else
	grdpt[i]->Draw("same p");
    }

  TLegend *lpd2 = new TLegend(0.25, 0.65, 0.70, 0.80,"","NDC");
  lpd2->SetBorderSize(0);
  lpd2->SetFillColor(0);
  lpd2->SetFillStyle(0);
  lpd2->AddEntry(geff[0], "MAPS(3)+TPC(60)", "p");
  lpd2->AddEntry(geff[1], "MAPS(3)+INTT(4)+TPC(60)", "p");
  lpd2->Draw();
  

  /*

  // now the upsilon spectra

  TFile *fin_ups[NPLOTS_UPS];

  fin_ups[0] = new TFile("look_quarkonia_histos_2mpas+4maps_out.root");  
  fin_ups[1] = new TFile("look_quarkonia_histos_2pixels+4maps_out.root");  

  TH1D *recomass[NPLOTS_UPS];

  for(int i=0;i<NPLOTS_UPS;i++)
    {
      fin_ups[i]->GetObject("recomass",recomass[i]);
    }

  TCanvas *cups = new TCanvas("cups","cups",50,50,1000,600);
  cups->SetLeftMargin(0.12);

  for(int i=0;i<NPLOTS_UPS;i++)
    {
      recomass[i]->SetMarkerColor(col[i]);
      recomass[i]->SetLineColor(col[i]);
      recomass[i]->SetMarkerStyle(20);
      recomass[i]->SetMarkerSize(1.2);

      if(i==0)
	recomass[i]->Draw("p");
      else
	recomass[i]->Draw("same p");
    }

  TLegend *lpq = new TLegend(0.15, 0.55, 0.40, 0.80,"Inner barrel","NDC");
  lpq->SetBorderSize(0);
  lpq->SetFillColor(0);
  lpq->SetFillStyle(0);
  lpq->AddEntry(recomass[0], "2 maps", "p");
  lpq->AddEntry(recomass[1], "2 pixels", "p");
  
  lpq->Draw();
  */

}
