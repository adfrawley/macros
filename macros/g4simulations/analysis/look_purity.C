#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include  <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

void look_purity()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  bool pt_resolution_fit = true;

  TFile *fin = new TFile("root_files/maps3+intt4+tpc60_purity_out.root");
  //TFile *fin = new TFile("root_files/maps3+tpc60_purity_out.root");

  //TFile *fin = new TFile("root_files/purity_out.root");  
  //TFile *fin = new TFile("ladder_maps_purity_out.root");  
  //TFile *fin = new TFile("cylinder_maps_purity_out.root");  

  if(!fin)
    {
      cout << "Failed to find input file" << endl;
      exit(1);
    }

  TCanvas *c1 = new TCanvas("c1","c1",5,5,800,600);
  c1->Divide(2,1);

  //==========================
  // 2D histo for embedded pions
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
  
  //========================================
  // extract DCA resolution vs pT from 2D histo hpt_dca2d
  // by fitting 80 pT slices

  TCanvas *c2 = new TCanvas("c2","c2",20,20,800,600);

  //int NPT = 79;
  int NPT = 70;
  double pT[80];
  double dca2d[80];
  for(int i = 0;i<NPT;i++)
    //for(int i = 10;i<11;i++)
    {
      double ptlo = (double) i * 0.5 + 0.25;
      double pthi = ptlo + 0.5;

      int binlo = hpt_dca2d->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_dca2d->GetXaxis()->FindBin(pthi);

      std::cout << "ptlo " << ptlo << " binlo " << binlo << " pthi " << " binhi " << binhi << std::endl;

      TH1D *h = new TH1D("h","dca2d resolution",2000, 0.0, 2.0);
      hpt_dca2d->ProjectionY("h",binlo,binhi);
      h->GetXaxis()->SetTitle("p_{T}");
      h->GetXaxis()->SetTitle("#Delta dca2d");
      h->GetXaxis()->SetTitleOffset(1.0);
      if(i>30) h->Rebin(4);
      h->DrawCopy();

      TF1 *f = new TF1("f","gaus");
      f->SetParameter(1,h->GetMean());
      f->SetParameter(2,h->GetRMS());
      h->Fit(f);

      pT[i] = (ptlo + pthi) / 2.0;

      dca2d[i] = f->GetParameter(2);
      cout << " pT " << pT[i] << " dca2d " << dca2d[i] << endl;
    }

  //============================================
  // plot the dca2d resolution extracted above for embedded pions

  TCanvas *c3 = new TCanvas("c3","c3",100,100,800,600);
  TGraph *grdca2d = new TGraph(NPT,pT,dca2d);
  grdca2d->SetMarkerStyle(20);
  grdca2d->SetMarkerSize(1.2);
  grdca2d->SetName("dca2d_resolution");
  grdca2d->SetTitle("dca2d resolution");

  TH1D *hdummy = new TH1D("hdummy","#Delta dca2d vs p_{T}",100,0.0,40.0);
  hdummy->SetMinimum(0);
  hdummy->SetMaximum(0.0050);
  hdummy->GetXaxis()->SetTitle("p_{T}");
  hdummy->GetYaxis()->SetTitle("#Delta dca2d");
  hdummy->GetYaxis()->SetTitleOffset(1.75);
  gPad->SetLeftMargin(0.15);
  hdummy->Draw();
  grdca2d->Draw("p");

  //===================================
  // extract pT resolution vs pT from hpt_compare

  TCanvas *c4 = new TCanvas("c4","c4",60,60,800,600);

  double dpT[80];

  for(int i = 0;i<NPT;i++)
  //for(int i = 10;i<11;i++)
    {
      double ptlo = (double) i * 0.5 + 0.25;
      double pthi = ptlo + 0.5;

      int binlo = hpt_compare->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_compare->GetXaxis()->FindBin(pthi);

      //TH1D *hpt = new TH1D("hpt","pT resolution ",500, -0.1, 0.1);    
      TH1D *hpt = new TH1D("hpt","pT resolution ",200, -0.1, 0.1);    
      hpt_compare->ProjectionY("hpt",binlo,binhi);
      hpt->GetXaxis()->SetTitle("#Delta p_{T}/p_{T}");
      hpt->GetXaxis()->SetTitleOffset(1.0);
      if(i>30) hpt->Rebin(4);
      hpt->DrawCopy();

      std::cout << "ptlo " << ptlo << " binlo " << binlo << " pthi " << pthi << " binhi " << binhi << " integral " << hpt->Integral() << std::endl;

      TF1 *f = new TF1("f","gaus");
      //f->SetParameter(1,hpt->GetMean());
      //f->SetParameter(2,hpt->GetRMS() );
      f->SetParameter(1,0.008);
      f->SetParameter(2,0.002 );
      hpt->Fit(f);

      pT[i] = (ptlo + pthi) / 2.0;

      dpT[i] = f->GetParameter(2);
      cout << " pT " << pT[i] << " dpT " << dpT[i] << endl;
    }

  //==========================================
  // Plot pT resolution extracted above for embedded pions

  TCanvas *c5 = new TCanvas("c5","c5",100,100,800,600);
  TGraph *grdpt = new TGraph(NPT,pT,dpT);
  grdpt->SetMarkerStyle(20);
  grdpt->SetMarkerSize(1.1);
  grdpt->SetName("pt_resolution");
  grdpt->SetTitle("pT resolution");

  TH1D *hdummy2 = new TH1D("hdummy2","#Delta p_{T} vs p_{T}",100,0.0,40.0);
  hdummy2->SetMinimum(0);
  hdummy2->SetMaximum(0.12);
  hdummy2->GetXaxis()->SetTitle("p_{T}");
  hdummy2->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");
  hdummy2->GetYaxis()->SetTitleOffset(1.2);
  hdummy2->Draw();
  grdpt->Draw("p");

  // Parameterize pT resolution
  
  //TF1 *fpt = new TF1("fpt","sqrt([0]*[0] + [1]*[1]*x*x)", 0, 10.0);
  TF1 *fpt = new TF1("fpt","sqrt([0]*[0] + [1]*[1]*x*x)", 0, 25.0);
  fpt->SetParameter(0,0.0054);
  fpt->SetParameter(1,0.0023);
  if(pt_resolution_fit)  
    grdpt->Fit(fpt,"R");

  char lab[1000];
  sprintf(lab,"#frac{#Deltap_{T}}{p_{T}} = #sqrt{%.4f^{2} + (%.6f #times p_{T})^{2}}", fpt->GetParameter(0), fpt->GetParameter(1));
  TLatex *mres = new TLatex(0.2,0.65,lab);
  //mres->SetTextSize(0.1);
  mres->SetNDC();
  if(pt_resolution_fit)  
    mres->Draw();

  // For making a cut on the momentum difference between rgpT and rpT 
  double pT_sigmas = 4.0;
  double const_term =  fpt->GetParameter(0);
  double linear_term = fpt->GetParameter(1);

  TH1D *hpt_truth;
  fin->GetObject("hpt_truth",hpt_truth);
  if(!hpt_truth)
    {
      cout << "Failed to get hpt_truth, quit!" << endl;
      exit(1);
    }

  TH1D *hpt_matched = new TH1D("hpt_matched","hpt_matched", 500, 0.0, 40.0);
  double eff_pt[80];

  for(int i = 0;i<NPT;i++)
    {
      double ptlo = (double) i * 0.5 + 0.25;
      double pthi = ptlo + 0.5;
      double ptval = (ptlo+pthi)/2.0;

      int binlo = hpt_dca2d->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_dca2d->GetXaxis()->FindBin(pthi);

      TH1D *hpt1 = new TH1D("hpt1","pT resolution ",2000, 0.0, 2.0);    
      hpt_compare->ProjectionY("hpt1",binlo,binhi);

      double ptres = sqrt(pow(const_term,2) + pow(linear_term * ptval,2)); 
      double ptreslo = 1.0 - ptres * pT_sigmas;
      double ptreshi = 1.0 + ptres * pT_sigmas;

      int momlo = hpt1->GetXaxis()->FindBin(ptreslo);
      int momhi = hpt1->GetXaxis()->FindBin(ptreshi);

      hpt_matched->Fill(ptval, hpt1->Integral(momlo, momhi));

      int tlo = hpt_truth->FindBin(ptval-0.1);
      int thi = hpt_truth->FindBin(ptval+0.1);
      double truth_yield = hpt_truth->Integral(tlo, thi);;
      eff_pt[i] = hpt1->Integral(momlo, momhi) / truth_yield;


      cout << " ptval " << ptval 
	   << " ptreslo " << ptreslo
	   << " ptreshi " << ptreshi
	   << " total " << hpt1->Integral()
	   << " momlo " << momlo
	   << " momhi " << momhi
	   << "  cut " << hpt1->Integral(momlo,momhi)
	   << " tlo  " << tlo
	   << " thi " << thi
	   << " truth " << truth_yield
	   << " eff_pt " << eff_pt[i]
	   << endl;

      delete hpt1;

    }  
  
  cout << " create canvas c7" << endl;
  TCanvas *c7 = new TCanvas("c7","c7",60,60,1000,600);

  TH1F *hd = new TH1F("hd","hd",100, 0.0, 40.0);
  hd->SetMinimum(0.0);
  hd->SetMaximum(1.0);
  hd->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hd->GetYaxis()->SetTitle("Single track efficiency");
  hd->Draw();

  TGraph *gr_eff = new TGraph(NPT,pT,eff_pt);
  gr_eff->SetName("single_track_efficiency");
  gr_eff->SetMarkerStyle(20);
  gr_eff->SetMarkerSize(1);
  gr_eff->SetMarkerColor(kRed);

  gr_eff->Draw("p");


  //=======================
  // Get track purity for Hijing events by 
  // finding tracks within 3 sigmas of the 
  // truth pT

  TH2D *hpt_hijing_compare = 0;
  fin->GetObject("hpt_hijing_compare",hpt_hijing_compare);
  if(!hpt_hijing_compare)
    {
      cout << "Did not get hpt_hijing_compare - quit!" << endl;
      exit(1);
    }

  static const int NVARBINS = 36;
  double xbins[NVARBINS+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
			      1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
			      2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
			      4.5,5.0,5.5,6.0,7.0,8.0,10.0};
  
  TH1D *hpt_hijing_allreco = new TH1D("hpt_hijing_allreco","hpt_hijing_allreco",NVARBINS,xbins);
  TH1D *hpt_hijing_matched = new TH1D("hpt_hijing_matched","hpt_hijing_matched",NVARBINS,xbins);
  hpt_hijing_matched->GetYaxis()->SetTitle("Purity");

  //double purity_pt[NVARBINS];
  
  for (int i=0;i<NVARBINS;i++)  
    {
      double ptlo = xbins[i];
      double pthi = xbins[i+1];
      double ptval = (ptlo+pthi)/2.0;

      int binlo = hpt_hijing_compare->GetXaxis()->FindBin(ptlo);
      int binhi = hpt_hijing_compare->GetXaxis()->FindBin(pthi);

      TH1D *hpt1 = new TH1D("hpt1","pT resolution ",2000, 0.0, 2.0);    
      hpt_hijing_compare->ProjectionY("hpt1",binlo,binhi);

      double ptres = sqrt(pow(const_term,2) + pow(linear_term * ptval,2)); 
      double ptreslo = 1.0 - ptres * pT_sigmas;
      double ptreshi = 1.0 + ptres * pT_sigmas;

      int momlo = hpt1->GetXaxis()->FindBin(ptreslo);
      int momhi = hpt1->GetXaxis()->FindBin(ptreshi);

      // This to avoid TEfficiency's stupid insistence that the histograms cannot be filled with weights
      for (int i=0;i<(int)hpt1->Integral(momlo,momhi);i++)
	hpt_hijing_matched->Fill(ptval);
      for (int i=0;i<(int)hpt1->Integral();i++)
	hpt_hijing_allreco->Fill(ptval);

      //purity_pt[i] = hpt1->Integral(momlo, momhi) / hpt1->Integral();
    }

  cout << " create canvas c7" << endl;
  TCanvas *cpur = new TCanvas("cpur","cpur",60,60,1000,600);

  TEfficiency* pEff = 0;
  //if(TEfficiency::CheckConsistency(*hpt_hijing_matched,*hpt_hijing_allreco))
    {
      pEff = new TEfficiency(*hpt_hijing_matched,*hpt_hijing_allreco);
      pEff->SetTitle("Reconstruction efficiency (3.0#sigma p_{T}) ; p_{T} (GeV/c) ; reco'd tracks within 3 #sigma p_{T}");
      pEff->SetMarkerStyle(20);
      pEff->SetMarkerColor(kBlack);
      pEff->SetMarkerSize(1);
      pEff->Draw();
      gPad->Update();
      pEff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.0);
      pEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0,1.1);
      pEff->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.1);
    }

  /*
  hd->SetMinimum(0.0);
  hd->SetMaximum(1.0);
  hd->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hd->GetYaxis()->SetTitle("Tracks within 3 #sigma of p_{T}");
  hd->Draw();

  TGraph *gr_pur = new TGraph(NVARBINS,xbins,purity_pt);
  gr_pur->SetName("Track purity");
  gr_pur->SetMarkerStyle(20);
  gr_pur->SetMarkerSize(1);
  gr_pur->SetMarkerColor(kRed);

  gr_pur->Draw("p");
  */

  //============================================
  // plot purity for Hijing tracks

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

  //=========================
  // plot quality for Hijing tracks

  TCanvas *cpq = new TCanvas("cpq","cpq",40,40,1200,600);
  TH1D * hquality;
  fin->GetObject("hquality",hquality);
  hquality->Draw();

  //======================
  // Plot number of hits per track
  TCanvas *c8 = new TCanvas("c8","c8",40,40,1200,600);
  TH1D * hnhits;
  fin->GetObject("hnhits",hnhits);
  hnhits->GetXaxis()->SetTitle("hits per track");
  hnhits->Draw();
  cout << "hnhits integral = " << hnhits->Integral() << endl;

  //=========================

  TCanvas *cdca = new TCanvas("cdca","cdca",5,20,1200,800);
  cdca->Divide(3,1);

  TH1D *hdca2d[3];
  fin->GetObject("hdca2d0",hdca2d[0]);
  fin->GetObject("hdca2d1",hdca2d[1]);
  fin->GetObject("hdca2d2",hdca2d[2]);

  cdca->cd(1);
  gPad->SetLogy(1);
  //gPad->SetRightMargin(0.1);
  hdca2d[0]->GetXaxis()->SetNdivisions(505);

  hdca2d[0]->Draw();

  cdca->cd(2);
  gPad->SetLogy(1);
  hdca2d[1]->GetXaxis()->SetNdivisions(505);
  hdca2d[1]->Draw();

  cdca->cd(3);
  gPad->SetLogy(1);
  hdca2d[2]->GetXaxis()->SetNdivisions(505);
  hdca2d[2]->Draw();

  TF1 *fdca = new TF1("fdca","gaus",-0.1,0.1);
  fdca->SetLineColor(kRed);
  cdca->cd(1);
  hdca2d[0]->Fit(fdca);
  char fitr[500];
  sprintf(fitr,"#splitline{p_{T} = 0.5-1.0 GeV/c}{#sigma = %.1f #mum}",fdca->GetParameter(2)*10000);
  TLatex *l1 = new TLatex(0.2,0.927,fitr);
  l1->SetNDC();
  l1->SetTextSize(0.07);
  l1->Draw();

  cdca->cd(2);
  hdca2d[1]->Fit(fdca);
  sprintf(fitr,"#splitline{p_{T} = 1.0-2.0 GeV/c}{#sigma = %.1f #mum}",fdca->GetParameter(2)*10000);
  TLatex *l2 = new TLatex(0.2,0.927,fitr);
  l2->SetNDC();
  l2->SetTextSize(0.07);
  l2->Draw();


  cdca->cd(3);
  hdca2d[2]->Fit(fdca);
  sprintf(fitr,"#splitline{p_{T} > 2.0 GeV/c}{#sigma = %.1f #mum}",fdca->GetParameter(2)*10000);
  TLatex *l3 = new TLatex(0.2,0.927,fitr);
  l3->SetNDC();
  l3->SetTextSize(0.07);
  l3->Draw();

  // Output some graphs for comparison plots

 TFile *fout = new TFile("root_files/look_purity_out.root","recreate");
  gr_eff->Write();
  grdca2d->Write();
  grdpt->Write();

}
