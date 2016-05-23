#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include  <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>


void look_purity()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile("purity_tpc_out.root");  
  //TFile *fin = new TFile("nlayers48_purity_tpc_out.root");  

  TCanvas *c1 = new TCanvas("c1","c1",5,5,800,600);
  c1->Divide(2,1);

  TH2D *hpt_dca2d = 0;
  fin->GetObject("hpt_dca2d",hpt_dca2d);
  c1->cd(1);
  gPad->SetLeftMargin(0.12);
  hpt_dca2d->GetYaxis()->SetTitleOffset(1.5);
  hpt_dca2d->Draw();

  TH2D *hpt_compare = 0;
  fin->GetObject("hpt_compare",hpt_compare);
  c1->cd(2);
  gPad->SetLeftMargin(0.12);
  hpt_compare->GetYaxis()->SetTitleOffset(1.5);
  hpt_compare->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",20,20,800,600);

  // extract DCA resolution vs pT
  int NPT = 79;
  double pT[80];
  double dca2d[80];
  for(int i = 0;i<NPT;i++)
    {
      double ptlo = (double) i * 0.5 + 0.5;
      double pthi = ptlo + 0.5;

      int binlo = hpt_dca2d->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_dca2d->GetXaxis()->FindBin(pthi);

      std::cout << "ptlo " << ptlo << " binlo " << binlo << " pthi " << " binhi " << binhi << std::endl;

      TH1D *h = new TH1D("h","dca2d resolution",2000, 0.0, 2.0);
      hpt_dca2d->ProjectionY("h",binlo,binhi);
      h->GetXaxis()->SetTitle("p_{T}");
      h->GetXaxis()->SetTitle("#Delta dca2d");
      h->GetXaxis()->SetTitleOffset(1.0);
      h->DrawCopy();

      TF1 *f = new TF1("f","gaus");
      f->SetParameter(1,h->GetMean());
      f->SetParameter(2,h->GetRMS());
      h->Fit(f);

      pT[i] = (ptlo + pthi) / 2.0;

      dca2d[i] = f->GetParameter(2);
      cout << " pT " << pT[i] << " dca2d " << dca2d[i] << endl;
    }

  TCanvas *c3 = new TCanvas("c3","c3",100,100,800,600);
  TGraph *grdca2d = new TGraph(NPT,pT,dca2d);
  grdca2d->SetMarkerStyle(20);
  grdca2d->SetMarkerSize(1.2);

  TH1D *hdummy = new TH1D("hdummy","#Delta dca2d vs p_{T}",100,0.0,40.0);
  hdummy->SetMinimum(0);
  hdummy->SetMaximum(0.0050);
  hdummy->GetXaxis()->SetTitle("p_{T}");
  hdummy->GetYaxis()->SetTitle("#Delta dca2d");
  hdummy->GetYaxis()->SetTitleOffset(1.75);
  gPad->SetLeftMargin(0.15);
  hdummy->Draw();
  grdca2d->Draw("p");

  // extract pT resolution vs pT

  TCanvas *c4 = new TCanvas("c4","c4",60,60,800,600);

  // extract pT resolution vs pT
  //-----------------------------------

  double dpT[80];

  for(int i = 0;i<NPT;i++)
    {
      double ptlo = (double) i * 0.5 + 0.5;
      double pthi = ptlo + 0.5;

      int binlo = hpt_dca2d->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_dca2d->GetXaxis()->FindBin(pthi);

      TH1D *hpt = new TH1D("hpt","pT resolution ",500, -0.1, 0.1);    
      hpt_compare->ProjectionY("hpt",binlo,binhi);
      hpt->GetXaxis()->SetTitle("#Delta p_{T}/p_{T}");
      hpt->GetXaxis()->SetTitleOffset(1.0);
      hpt->DrawCopy();

      std::cout << "ptlo " << ptlo << " binlo " << binlo << " pthi " << " binhi " << binhi << " integral " << hpt->Integral() << std::endl;

      TF1 *f = new TF1("f","gaus");
      f->SetParameter(1,hpt->GetMean());
      f->SetParameter(2,hpt->GetRMS());
      hpt->Fit(f);

      pT[i] = (ptlo + pthi) / 2.0;

      dpT[i] = f->GetParameter(2);
      cout << " pT " << pT[i] << " dpT " << dpT[i] << endl;
    }

  TCanvas *c5 = new TCanvas("c5","c5",100,100,800,600);
  TGraph *grdpt = new TGraph(NPT,pT,dpT);
  grdpt->SetMarkerStyle(20);
  grdpt->SetMarkerSize(1.1);

  TH1D *hdummy2 = new TH1D("hdummy2","#Delta p_{T} vs p_{T}",100,0.0,40.0);
  hdummy2->SetMinimum(0);
  hdummy2->SetMaximum(0.12);
  hdummy2->GetXaxis()->SetTitle("p_{T}");
  hdummy2->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");
  hdummy2->GetYaxis()->SetTitleOffset(1.2);
  hdummy2->Draw();
  grdpt->Draw("p");

  TCanvas *c6 = new TCanvas("c6","c6",40,40,1200,600);
  c6->Divide(3,1);

  c6->cd(1);
  TH1D * hpurity[4];
  fin->GetObject("hpurity0",hpurity[0]);
  fin->GetObject("hpurity1",hpurity[1]);
  fin->GetObject("hpurity2",hpurity[2]);
  fin->GetObject("hpurity3",hpurity[3]);

  hpurity[0]->Draw();
  for(int i=1;i<4;i++)
    hpurity[i]->Draw("same");

  c6->cd(2);
  TH1D * hpurity_quality[4];
  fin->GetObject("hpurity_quality0",hpurity_quality[0]);
  fin->GetObject("hpurity_quality1",hpurity_quality[1]);
  fin->GetObject("hpurity_quality2",hpurity_quality[2]);
  fin->GetObject("hpurity_quality3",hpurity_quality[3]);

  hpurity_quality[0]->SetMaximum(100000);
  hpurity_quality[0]->GetXaxis()->SetTitle("quality");
  hpurity_quality[0]->Draw();
  for(int i=1;i<4;i++)
    hpurity_quality[i]->Draw("same");

  c6->cd(3);
  TH1D * hpurity_dca[4];
  fin->GetObject("hpurity_dca0",hpurity_dca[0]);
  fin->GetObject("hpurity_dca1",hpurity_dca[1]);
  fin->GetObject("hpurity_dca2",hpurity_dca[2]);
  fin->GetObject("hpurity_dca3",hpurity_dca[3]);

  //hpurity_dca[0]->SetMaximum(100000);
  hpurity_dca[0]->GetXaxis()->SetTitle("dca");
  hpurity_dca[0]->Draw();
  for(int i=1;i<4;i++)
    hpurity_dca[i]->Draw("same");

  TLegend *leg = new TLegend(0.35,0.4,0.88,0.7,"","NDC");
  leg->AddEntry(hpurity_dca[0], "0 noise hits", "l");
  leg->AddEntry(hpurity_dca[1], "1 noise hit", "l");
  leg->AddEntry(hpurity_dca[2], "2 noise hits", "l");
  leg->AddEntry(hpurity_dca[3], "3 noise hits", "l");
  c6->cd(2);
  leg->Draw();

  /*
  c6->cd(2);
  TH1D * hquality;
  fin->GetObject("hquality",hquality);
  hquality->Draw();
  */



}
