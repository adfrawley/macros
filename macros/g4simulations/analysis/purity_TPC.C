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
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

using namespace std;

void purity_TPC()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  //=============================================================================
  // Important! Change the resolution so that the 3-sigma cut is correct
  //=============================================================================  

  // revised MIE reference design  
  double pT_res_constant = 0.0108;
  double pT_res_linear = 0.000271;

  //===============================================================================

  bool verbose = true;

  // Setup parameters
  //=================

  // should normally be true
  bool use_reco_pt = true;
  
  double mom_rescale = 1.0;

  static const int max_layers = 64;  // maximum number of tracking layers
  
  double ptmax = 12.2;
  
  double pT_sigmas = 3.0;  // number of sigmas in pT to use for track evaluation
  
  //based on plots of purity vs DCA sigmas and vs quality
  //double quality_cut = 3.0;
  double quality_cut = 1.0;
  double dca_sigma_cut = 3.0;
  bool use_dca_sigmas = false;
  double dca_cut = 0.1; // 1mm, used only if use_dca_sigmas is false
  //===============================
  
  // Open the evaluator output file 
  cout << "Reading ntuple " << endl;


   TChain* ntp_track = new TChain("ntp_track","reco tracks");
  TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
  TChain* ntp_vertex = new TChain("ntp_vertex","events");
  TChain *ntp_cluster = new TChain("ntp_cluster","clusters");

 // The condor jobs make 1000 files
  for(int i=0;i<1000;i++)
    {
      char name[500];
      sprintf(name,"../eval_output/g4svx_eval_%i.root",i);
      ntp_vertex->Add(name);
      ntp_track->Add(name);
      ntp_gtrack->Add(name);
    }

  // This include file contains the list of files to chain
  // #include "ntuple_files_tpc.C"
  
  // This include file contains the definitions of the ntuple variables, and the chain definitions
#include "ntuple_variables.C"
    
  // Get the purity for all tracks in ntp_track
  
  static const int npurity = 4;
  TH1D *hpurity[npurity];
  TH1D *hpurity_quality[npurity];
  TH1D *hpurity_dca[npurity];
  
  for(int i=0;i<npurity;i++)
    {
      char hname[500];
      sprintf(hname,"hpurity%i",i);
      hpurity[i] = new TH1D(hname,hname,100,0.0,ptmax);
      
      sprintf(hname,"hpurity_quality%i",i);
      hpurity_quality[i] = new TH1D(hname,hname,50,0.0,4.0);
      
      sprintf(hname,"hpurity_dca%i",i);
      hpurity_dca[i] = new TH1D(hname,hname,50,-0.04,0.04);
    }


  static const int NPTDCA = 3;

  TH1D *hdca2d[NPTDCA];
  for(int ipt=0;ipt<NPTDCA;ipt++)
    {  
      char hname[500];
      sprintf(hname,"hdca2d%i",ipt);

      hdca2d[ipt] = new TH1D(hname,hname,2000,-0.1,0.1);
      hdca2d[ipt]->GetXaxis()->SetTitle("DCA (cm)");
      hdca2d[ipt]->GetXaxis()->SetTitleSize(0.055);
      hdca2d[ipt]->GetXaxis()->SetLabelSize(0.055);

      hdca2d[ipt]->GetYaxis()->SetTitleSize(0.06);
      hdca2d[ipt]->GetYaxis()->SetLabelSize(0.055);
    }

  /*
  static const int NPTEMBED = 20;
  TH1D *hdca_embed[NPTEMBED];
  for(int ipt=0;ipt<NPTEMBED;ipt++)
    {  
      char hname[500];
      sprintf(hname,"hdca_embed%i",ipt);

      hdca_embed[ipt] = new TH1D(hname,hname,2000,-0.1,0.1);
      hdca_embed[ipt]->GetXaxis()->SetTitle("DCA (cm)");
      hdca_embed[ipt]->GetXaxis()->SetTitleSize(0.055);
      hdca_embed[ipt]->GetXaxis()->SetLabelSize(0.055);

      hdca_embed[ipt]->GetYaxis()->SetTitleSize(0.06);
      hdca_embed[ipt]->GetYaxis()->SetLabelSize(0.055);
    }
  */



  TH1D *hdca2dsigma = new TH1D("hdca2dsigma","hdca2dsigma",2000,-5.0,5.0);
  hdca2dsigma->GetXaxis()->SetTitle("dca2d / dca2dsigma");

  TH1D *hZdca = new TH1D("hZdca","hZdca",2000,-0.5,0.5);
  hZdca->GetXaxis()->SetTitle("Z DCA (cm)");

  TH1D *hquality = new TH1D("hquality","hquality",2000,0.0,5.0);
  hquality->GetXaxis()->SetTitle("quality");

  TH1D *heta = new TH1D("hea","heta",100,-1.2,1.2);

  static const int NVARBINS = 36;
  double xbins[NVARBINS+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
			    1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
			    2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
			      4.5,5.0,5.5,6.0,7.0,8.0,10.0};

  TH1D *hpt_allreco = new TH1D("hpt_allreco","hpt_allreco",NVARBINS,xbins);
  TH1D *hpt_goodreco = new TH1D("hpt_goodreco","hpt_goodreco",NVARBINS,xbins);

  TH2D *hpt_compare = new TH2D("hpt_compare","hpt_compare",500, 0.0, 40.0, 2000, 0.0, 2.0);
  hpt_compare->GetXaxis()->SetTitle("p_{T}");
  hpt_compare->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");

  TH2D *hpt_dca2d = new TH2D("hpt_dca2d","hpt_dca2d",500, 0.0, 40.0,500, -0.04, 0.04);
  hpt_dca2d->GetXaxis()->SetTitle("p_{T}");
  hpt_dca2d->GetYaxis()->SetTitle("dca2d");

  //TH2D *hpt_dcaz = new TH2D("hpt_dcaz","hpt_dcaz",500, 0.0, 40.0,500, -0.1, 0.1);

  TH1D *hptreco[NVARBINS];
  for(int ipt=0;ipt<NVARBINS;ipt++)
    {
      char hn[1000];
      sprintf(hn,"hptreco%i",ipt);
      hptreco[ipt] = new TH1D(hn, hn, NVARBINS, xbins);
    }



  //============================================================
  // Loop over events 
  //   Reject events with bad event vertex
  //   Loop over reco'd tracks
  //     Make quality cut (usually < 3)
  //       Fill dca2d histos
  //       Make loose dca2d cut (usually < 1 mm)
  //         Fill Z dca histo
  //         Drop tracks outside 1 mm in Z dca from evt vertex
  //           Add this reco'd track to hpt_allreco
  //           Make a cut on rpT < pT_sigmas from rgpT
  //             fill hpt_goodreco
  //   Make a TEfficiency of hpt_goodreco / hpt_allreco
  //============================================================
  int nr = 0;

  for(int iev=0;iev<ntp_vertex->GetEntries();iev++)
    {
      if(iev%100 == 0)
	cout << "Get event " << iev << endl;

      int recoget = ntp_vertex->GetEntry(iev);
      double evt_vertex_z = evz;

      for(int ir=nr;ir<nr+ntracks;ir++)
	{
	  int recoget = ntp_track->GetEntry(ir);

	  //cout << " rnhits = " << rnhits << " rpurity = " << rpurity << endl;
	  // want only tracks with hits in all layers
	  //if(rnhits != max_layers)
	  //continue;

	  /*
	  if(rgembed > 0)
	    {
	      // embedded pion
	      // we want to get the dca istribution vs pT for embedded pions

	      // get the pT and find the corresponding array index

	      double rgpT = sqrt(rgpx*rgpx+rgpy*rgpy);	  	  

	      int index = (int) (2.0 * rgpT + 0.01) - 1;

	      //cout << " rgpT " << rgpT  << " index " << index  << endl;

	      if(rquality < quality_cut)
		hdca_embed[index]->Fill(rdca2d);
	    }
	  */

	  /*
	  // The embedded pions screw this plot up, so ignore them
	  if(rgembed == 1)
	   continue;
	  */

	  double rgpT = sqrt(rgpx*rgpx+rgpy*rgpy);	  	  
	  double rpT = sqrt(rpx*rpx+rpy*rpy);

	  // Does the track pass the dca2d cut? 
	  // If so, add to hquality
	  if(use_dca_sigmas)
	    {
	      if(rpurity == max_layers && rgpT > 0.5 && rdca2d/rdca2dsigma < dca_sigma_cut)
		if(rgpT > 0.5 && rdca2d/rdca2dsigma < dca_sigma_cut)
		  {
		    hquality->Fill(rquality);
		  }
	    }
	  else
	    {
	      if(rpurity == max_layers && rgpT > 0.5 && fabs(rdca2d) < dca_cut)
		if(rgpT > 0.5 && fabs(rdca2d) < dca_cut)
		  {
		    hquality->Fill(rquality);
		  }
	    }
	  

	  // Make a histogram of purity vs quality or dca sigmas with no cuts
	  int ipurity = max_layers - (int) rpurity;
	  if(ipurity > 3)
	    {
	      continue;
	    }

	  hpurity_quality[ipurity]->Fill(rquality);
	  hpurity_dca[ipurity]->Fill(rdca2d);

	  // Now make some cuts on dca and quality

	  if(rquality > quality_cut)
	    {
	      continue;
	    }

	  // The track passed the quality cut, add it to the dca2d histos
	  if(rquality < quality_cut)
	    {
	      if(rgpT > 0.5 && rgpT <= 1.0)
		hdca2d[0]->Fill(rdca2d);
	      else if(rgpT > 1.0 && rgpT <= 2.0)
		hdca2d[1]->Fill(rdca2d);
	      else if(rgpT > 2.0)
		hdca2d[2]->Fill(rdca2d);

	      if(rgpT > 0.5)
		hdca2dsigma->Fill(rdca2d/rdca2dsigma);
	    }

	  if(use_dca_sigmas)
	    {
	      if( fabs(rdca2d) / rdca2dsigma > dca_sigma_cut )
		continue;
	    }
	  else
	    {
	      if(fabs(rdca2d) > dca_cut)
		continue;
	    }

	  // Histogram the track Z DCA
	  hZdca->Fill(rvz - evz);

	  // drop tracks outside 1 mm in Z dca from the event vertex    
	  if(fabs(rvz - evz) > 0.1)
	    {	      
	      continue;
	    }

	  // histogram the reco pT for the tracks that pass the cuts
	  hpt_allreco->Fill(rpT);
	  hpt_compare->Fill(rgpT,rpT/rgpT);
	  hpt_dca2d->Fill(rgpT, rdca2d);
	  //hpt_dcaz->Fill(rgpT, vz-evet_vertex_z);

	  for(int ipt=0;ipt<NVARBINS-1;ipt++)
	    {
	      if(rpT > xbins[ipt] && rpT < xbins[ipt+1])
		hptreco[ipt]->Fill(rpT);
	    }

	  // histogram the reco pT for tracks that:
	  //   - have a reco momentum that is within +/- n sigmas of the G4 pT

	  // make a cut on the momentum difference between rgpT and rpT
	  double ptratio = rpT/rgpT;

	  // The mom res is 1.08% constant + 2.71 E-04 linear
	  double ptres = pT_res_constant + pT_res_linear * rgpT; 
	  ptres = ptres * pT_sigmas;

	  if( fabs(ptratio - 1.0) < ptres)
	    {
	      hpt_goodreco->Fill(rpT);
	    }

	  // make histograms of the purity for the tracks that passed the quality and dca cuts
	  if(ipurity >= 0 && ipurity <= npurity)
	    {
	      if(use_reco_pt)
		hpurity[ipurity]->Fill(rpT);
	      else
		hpurity[ipurity]->Fill(rgpT);
	    }
	  /*
	  // print out high momentum tracks to look at details
	  if(rpT > 6.0)
	    cout << " rpurity " << rpurity
		 << " rpT " << rpT
		 << " rgpT " << rgpT
		 << " ptratio " << ptratio
		 << " pT sigmas " << (ptratio - 1.0)/ptres
		 << " rgfx " << rgfx
		 << " rgfy " << rgfy
		 << " rgfz " << rgfz
		 << " rnhits " << rnhits
		 << " quality " << rquality
		 << " dca2d " << rdca2d
		 << " dca2dsigma " << rdca2dsigma
		 << endl; 
	  */

	  double eta = asinh(rpz/sqrt(rpx*rpx+rpy*rpy));
	  heta->Fill(eta);

	}  // end loop over reco'd tracks
      
      nr += ntracks;

    }  // end loop over events
  
  TCanvas *cc=new TCanvas("cc","cc",5,5,700,500);  
  cc->Divide(1,2);
  cc->cd(1);
  //gPad->SetLogy(1);
  hdca2dsigma->Draw();
  cc->cd(2);
  hquality->Draw();
  cc->cd(3);

  TCanvas *cdca = new TCanvas("cdca","cdca",5,20,1200,800);
  cdca->Divide(3,1);

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


  TCanvas *cZdca = new TCanvas("cZdca","cZdac",30,10,600,600);
  hZdca->SetMinimum(0.5);
  gPad->SetLogy(1);
  hZdca->Draw();

  TCanvas *ceta = new TCanvas("ceta","ceta",30,10,600,600);
  heta->SetMinimum(0.0);
  //gPad->SetLogy(1);
  heta->Draw();

  TCanvas *cpt = new TCanvas("cpt","cpt",20,20,600,600);
  cpt->SetLeftMargin(0.13);
  char tstr[500];
  sprintf(tstr,"Reconstructed track purity;p_{T} (GeV/c);Reco'd pT < %.1f #sigma",pT_sigmas);
  TEfficiency* pEff = 0;
  //hpt_goodreco->Rebin(4);
  //hpt_allreco->Rebin(4);
  if(TEfficiency::CheckConsistency(*hpt_goodreco,*hpt_allreco))
    {
      pEff = new TEfficiency(*hpt_goodreco,*hpt_allreco);
      pEff->SetTitle(tstr);
      pEff->SetMarkerStyle(20);
      pEff->SetMarkerColor(kBlack);
      pEff->SetMarkerSize(1);
      pEff->Draw();
      gPad->Update();
      pEff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
      pEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0,1.1);
    }
  

  // Make the purity plot for all tracks 
  //==================================== 

  TCanvas *cpurity = new TCanvas("cpurity","cpurity",50,50,1000,700);
  cpurity->Divide(2,1);
  cpurity->SetLeftMargin(0.12);
  
  int col[8] = {kCyan+1,kGreen+1,kMagenta,kRed,kOrange,kBlack,kViolet,kBlue};
  //int col[7] = {kCyan+1,kGreen+1,kMagenta,kRed,kOrange,kBlack,kBlue};
  
  for(int i=npurity-1;i>-1;i--)
    {
      hpurity[i]->SetMarkerColor(col[i]);
      hpurity[i]->SetMarkerStyle(20);
      hpurity[i]->SetMarkerSize(1);
      hpurity[i]->SetLineColor(col[i]);
      hpurity[i]->SetLineWidth(2);
      hpurity[i]->SetMinimum(0.6);

      hpurity_quality[i]->SetMarkerColor(col[i]);
      hpurity_quality[i]->SetMarkerStyle(20);
      hpurity_quality[i]->SetMarkerSize(1);
      hpurity_quality[i]->SetLineColor(col[i]);
      hpurity_quality[i]->SetLineWidth(2);
      hpurity_quality[i]->SetMinimum(1.0);
      //hpurity_quality[i]->SetMaximum(70000000.0);

      hpurity_dca[i]->SetMarkerColor(col[i]);
      hpurity_dca[i]->SetMarkerStyle(20);
      hpurity_dca[i]->SetMarkerSize(1);
      hpurity_dca[i]->SetLineColor(col[i]);
      hpurity_dca[i]->SetLineWidth(2);
      hpurity_dca[i]->SetMinimum(10.0);
      //hpurity_dca[i]->SetMaximum(4000000.0);

    }

  // plot the purity histograms vs pT
  cpurity->cd(1);
  gPad->SetLogy(1);
  for(int i=npurity-1;i>-1;i--)
    {
      if(i==npurity-1)
        {
          hpurity[i]->GetXaxis()->SetRangeUser(0.3,ptmax);
          hpurity[i]->GetYaxis()->SetTitle("Tracks");
          hpurity[i]->GetYaxis()->SetTitleOffset(1.3);
	  if(use_reco_pt)
	    hpurity[i]->GetXaxis()->SetTitle("p_{T}");
	  else
	    hpurity[i]->GetXaxis()->SetTitle("G4 p_{T}");                                                                                  
          hpurity[i]->DrawCopy("p");
        }
      else
        hpurity[i]->DrawCopy("p same");
    }

  TLegend *lpurity = new TLegend(0.15,0.2,0.45,0.4,"track purity","NDC");
  lpurity->SetBorderSize(1);
  lpurity->SetFillColor(0);
  lpurity->SetFillStyle(0);
  lpurity->AddEntry(hpurity[3],"3 noise tracks","p");
  lpurity->AddEntry(hpurity[2],"2 noise tracks","p");
  lpurity->AddEntry(hpurity[1],"1 noise track","p");
  lpurity->AddEntry(hpurity[0],"0 noise tracks","p");

  lpurity->Draw();

  char cutstr[500];
  if(use_dca_sigmas)
    sprintf(cutstr,"#splitline{quality < %.2f}{dca cut < %.2f #sigma}",quality_cut,dca_sigma_cut);
  else
    sprintf(cutstr,"#splitline{quality < %.2f}{dca cut < %.2f cm}",quality_cut,dca_cut);

  TLatex *labp1 = new TLatex(0.60,0.2,cutstr);
  labp1->SetNDC();
  labp1->SetTextSize(0.035);
  labp1->Draw();

  // make the plots of purity fraction vs pT                                                                                             

  cpurity->cd(2);
  gPad->SetLogy(1);

  TH1D *htotal2 = (TH1D*) hpurity[npurity-1]->Clone("htotal2");

  for(int i=npurity-2;i>-1;i--)
    {
      htotal2->Add(hpurity[i]);
    }

  // get the fractional purity                                                                                                           
  for(int i=npurity-1;i>-1;i--)
    {
      hpurity[i]->Divide(htotal2);
    }

  for(int i=0;i<npurity;i++)
    {
      if(i==0)
        {
          hpurity[i]->GetXaxis()->SetRangeUser(0.3,ptmax);
          hpurity[i]->GetYaxis()->SetTitle("fraction");
          hpurity[i]->GetYaxis()->SetTitleOffset(1.3);

	  hpurity[i]->SetMaximum(1.02);
	  hpurity[i]->SetMinimum(0.000005);
	  if(use_reco_pt)
	    hpurity[i]->GetXaxis()->SetTitle("p_{T}");
	  else
	    hpurity[i]->GetXaxis()->SetTitle("G4 p_{T}");                                                                                  
          hpurity[i]->DrawCopy("p");
        }
      else
        hpurity[i]->DrawCopy("p same");
    }


  lpurity->Draw();

  TLatex *labp2 = new TLatex(0.60,0.25,cutstr);
  labp2->SetNDC();
  labp2->SetTextSize(0.035);
  labp2->Draw();


  TCanvas *cpq = new TCanvas("cpq","cpq",10,10,600,600);
  cpq->SetLeftMargin(0.12);
  gPad->SetLogy(1);
  
  for(int i=npurity-1;i>-1;i--)
    {
      if(i==max_layers-1)
        {
	  hpurity_quality[i]->GetXaxis()->SetTitle("quality");
          hpurity_quality[i]->Draw("p");
        }
      else
        hpurity_quality[i]->Draw("p same");
    }

  TLegend *lpq = new TLegend(0.42,0.55,0.90,0.90,"track purity - 5000 0-4fm Hijing","NDC");
  lpq->SetBorderSize(1);
  lpq->SetFillColor(0);
  lpq->SetFillStyle(0);
  lpq->AddEntry(hpurity[3],"3noise hits","p");
  lpq->AddEntry(hpurity[2],"2 noise hits","p");
  lpq->AddEntry(hpurity[1],"1 noise hit","p");
  lpq->AddEntry(hpurity[0],"0 noise hits","p");
  lpq->Draw();

  TCanvas *cptreco = new TCanvas("cptreco","cptreco",5,5,600,600);
  cptreco->Divide(2,1);

  cptreco->cd(1);
  hpt_compare->Draw();

  cptreco->cd(2);
  hpt_dca2d->Draw();

  /*
  TF1 * fit = new TF1("fit","gaus");
  for(int ipt=0;ipt<NVARBINS;ipt++)
    {
      cptreco->cd(ipt+1);
      hptreco[ipt]->Draw();

      fit->SetParameter(1,hptreco[ipt]->GetMean());
      hptreco[ipt]->Fit(fit);
    }
  */




  /*
  TCanvas *cpd = new TCanvas("cpd","cpd",600,10,600,600);
  cpd->SetLeftMargin(0.12);
  gPad->SetLogy(1);
  
  for(int i=nlayers-1;i>-1;i--)
    {
      if(i==nlayers-1)
        {
	  hpurity_dca[i]->GetXaxis()->SetTitle("|dca2d|/dca2dsigma");
          hpurity_dca[i]->Draw("p");
        }
      else
        hpurity_dca[i]->Draw("p same");
    }

  TLegend *lpd = new TLegend(0.42,0.55,0.90,0.90,"track purity - 5000 0-4fm Hijing","NDC");
  lpd->SetBorderSize(1);
  lpd->SetFillColor(0);
  lpd->SetFillStyle(0);
  if(nlayers == 8)
    lpd->AddEntry(hpurity[7],"purity = 8","p");
  lpd->AddEntry(hpurity[6],"purity = 7","p");
  lpd->AddEntry(hpurity[5],"purity = 6","p");
  lpd->AddEntry(hpurity[4],"purity = 5","p");
  lpd->AddEntry(hpurity[3],"purity = 4","p");
  lpd->AddEntry(hpurity[2],"purity = 3","p");
  lpd->AddEntry(hpurity[1],"purity = 2","p");
  lpd->AddEntry(hpurity[0],"purity = 1","p");
  lpd->Draw();

  TLatex *labpd = new TLatex(0.68,0.68,"8 mm strips");
  //TLatex *labpd = new TLatex(0.67,0.68,"16 mm strips");
  labpd->SetNDC();
  labpd->SetTextSize(0.04);
  labpd->Draw();
  */




  // Output the results

  TFile *fout = new TFile("purity_tpc_out.root","recreate");

  hdca2d[0]->Write();
  hdca2d[1]->Write();
  hdca2d[2]->Write();

  hpt_compare->Write();
  hpt_dca2d->Write();

  for(int ipt=0;ipt<NVARBINS;ipt++)
    {
      hptreco[ipt]->Write();
    }

  hquality->Write();

  for(int i=0;i<npurity;i++)
    {
      hpurity[i]->Write();
      hpurity_dca[i]->Write();
      hpurity_quality[i]->Write();
    }

  hZdca->Write();

  //  for(int ipt=0;ipt<NPTEMBED;ipt++)
  //hdca_embed[ipt]->Write();
  
  fout->Close();


}

