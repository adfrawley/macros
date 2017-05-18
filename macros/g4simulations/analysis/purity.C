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

void purity()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  //=============================================================================
  // Important! Change the resolution so that the 3-sigma cut is correct
  //=============================================================================  

  // revised MIE reference design  
  //double pT_res_constant = 0.0108;
  //double pT_res_linear = 0.000271;

  //===============================================================================

  bool verbose = true;

  // Setup parameters
  //=================

  // should normally be true
  bool use_reco_pt = true;
  double mom_rescale = 1.0;
  
  // This should be inner maps layers (3) + INTT(4)+TPC (60)
  static const int nlayers = 67;  // maximum number of tracking layers for tpc40+intt4+maps3
  static const int nmissed = 4;  // maximum number of missed layers to allow
  double ptmax = 12.2;
  
  //double pT_sigmas = 4.0;  // number of sigmas in pT to use for track evaluation
  
  //based on plots of purity vs DCA  and vs quality
  //double quality_cut = 1.5;
  double quality_cut = 40.0;
  double dca_cut = 0.1; // 1mm
  //===============================

  /*  
  // Open the evaluator output file 
  cout << "Reading ntuple " << endl;

  TChain* ntp_track = new TChain("ntp_track","reco tracks");
  TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
  TChain* ntp_vertex = new TChain("ntp_vertex","events");
  TChain *ntp_cluster = new TChain("ntp_cluster","clusters");
  */
    

  //=======================================
  // define histograms

  TH1D *hnhits = new TH1D("hnhits","hnhits",100,0,99);    

  // Get the purity for all tracks in ntp_track
  
  static const int npurity = nmissed;
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

  TH1D *hdca2dsigma = new TH1D("hdca2dsigma","hdca2dsigma",2000,-5.0,5.0);
  hdca2dsigma->GetXaxis()->SetTitle("dca2d / dca2dsigma");

  TH1D *hZdca = new TH1D("hZdca","hZdca",2000,-0.5,0.5);
  hZdca->GetXaxis()->SetTitle("Z DCA (cm)");

  TH1D *hquality = new TH1D("hquality","hquality",2000,0.0,20.0);
  hquality->GetXaxis()->SetTitle("quality");

  TH1D *hgeta = new TH1D("hgeta","hgeta",100,-1.2,1.2);
  TH1D *hreta = new TH1D("hreta","hreta",100,-1.2,1.2);

  static const int NVARBINS = 36;
  double xbins[NVARBINS+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
			    1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
			    2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
			      4.5,5.0,5.5,6.0,7.0,8.0,10.0};

  TH1D *hpt_truth = new TH1D("hpt_truth","hpt_truth",500, 0.0, 50.0);
  TH2D *hpt_compare = new TH2D("hpt_compare","hpt_compare",500, 0.0, 50.0, 2000, 0.0, 2.0);
  hpt_compare->GetXaxis()->SetTitle("p_{T}");
  hpt_compare->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");

  TH2D *hpt_dca2d = new TH2D("hpt_dca2d","hpt_dca2d",500, 0.0, 50.0,500, -0.04, 0.04);
  hpt_dca2d->GetXaxis()->SetTitle("p_{T}");
  hpt_dca2d->GetYaxis()->SetTitle("dca2d");

  TH1D *hpt_hijing_truth = new TH1D("hpt_hijing_truth","hpt_hijing_truth",500, 0.0, 10.0);
  TH2D *hpt_hijing_compare = new TH2D("hpt_hijing_compare","hpt_hijing_compare",500, 0.0, 10.0, 2000, 0.0, 2.0);
  hpt_hijing_compare->GetXaxis()->SetTitle("p_{T}");
  hpt_hijing_compare->GetYaxis()->SetTitle("#Delta p_{T}/p_{T}");

  TH1D *hptreco[NVARBINS];
  for(int ipt=0;ipt<NVARBINS;ipt++)
    {
      char hn[1000];
      sprintf(hn,"hptreco%i",ipt);
      hptreco[ipt] = new TH1D(hn, hn, NVARBINS, xbins);
    }

  //TH1D *h_evt_dca2d = new TH1D("h_evt_dca2d","h_evt_dca2d",1000,-0.001,0.001);

  // end of histogram definitions
  //===================================================

  //============================================================
  // Loop over events 
  //   Loop over reco'd tracks
  //     Make quality cut
  //       Fill dca2d histos
  //       Make loose dca2d cut (usually < 1 mm)
  //         Fill Z dca histo
  //         Drop tracks outside 1 mm in Z dca from evt vertex
  //           Add this reco'd track 2D histo of (rpT-pT) vs pT
  //============================================================

  // The condor job output files
  for(int i=0;i<1000;i++)
    {
      // Open the evaluator output file 
      cout << "Reading ntuple " << endl;
      
      TChain* ntp_track = new TChain("ntp_track","reco tracks");
      TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
      TChain* ntp_vertex = new TChain("ntp_vertex","events");
      TChain *ntp_cluster = new TChain("ntp_cluster","clusters");
      
      // This include file contains the definitions of the ntuple variables, and the chain definitions
#include "ntuple_variables.C"

      char name[500];
      sprintf(name,"../eval_output/g4svtx_eval_%i.root",i);
      cout << "Adding file number " << i << " with name " << name << endl;
      ntp_vertex->Add(name);
      ntp_track->Add(name);
      ntp_gtrack->Add(name);

      // skip this file if there are no tracks 
      if(! ntp_gtrack->GetEntries())
	continue;

      if(verbose> 0)  
	cout <<  " ntp_vertex entries: " << ntp_vertex->GetEntries()
	     << " ntp_gtrack entries: " << ntp_gtrack->GetEntries()
	     << " ntp_track entries: " << ntp_track->GetEntries()
	     << endl;
      
      int nr = 0;
      int ng = 0;
      //for(int iev=0;iev<ntp_vertex->GetEntries();iev++)
      for(int iev=0;iev<1;iev++)
	{
	  if(verbose) cout << " iev = " << iev << " ng " << ng << " nr " << nr << endl;  
	  /*
	    int recoget = ntp_vertex->GetEntry(iev);
	    if(!recoget)
	    {
	    cout << "Failed to get ntp_vertex entry " << iev << endl;
	    exit(1);
	    }
	    
	    if(iev%20 == 0)
	    cout << "Get event " << iev << " with vertex x " << evx
	    << " vertex y " << evy << " vertex z " << evz 
	    << " ngtracks " << ngtracks << " ntracks " << ntracks 
	    << " ng " << ng << " nr " << nr 
	    << endl;
	  */

	  // expect 100 tracks or less per event
	  ngtracks = 110;

	  int naccept = 0;
	 	  
	  // ngtracks is defined in ntuple_variables.C and is the number of g4 tracks
	  // ntracks is defined in ntuple_variables.C and is the number of reco'd tracks
	  
	  int n_embed_gtrack = 0;            
	  int n_embed_gtrack_1 = 0;            
	  if(verbose > 0) cout << "Process truth tracks:" << endl;
	  for(int ig=ng;ig<ng+ngtracks;ig++)
	    {
	      int recoget = ntp_gtrack->GetEntry(ig);
	      if(recoget == 0)
		{
		  cout << "Failed to get ntp_gtrack entry " << ig << " in ntp_gtrack" << endl;
		  //exit(1);
		  //continue;
		  break;
		}
	      
	      if(tevent != iev)
		{
		  if(verbose) cout << " change of event, new tevent = " << tevent << endl;
		  break;
		}
	      
	      ng++;
	      
	      // only embedded tracks
	      if(!tembed == 1)
		continue;
	      
	      if(verbose > 0) cout << " ig = " << ig << " tevent " << tevent << " tgtrackid " << tgtrackid << endl;
	      
	      if(tembed == 1)
		n_embed_gtrack++;
	      
	      // skip tracks that do not pass through all layers (judged using truth track)
	      if(tnhits < nlayers-nmissed)
		continue;
	      
	      if(tembed == 1)
		n_embed_gtrack_1++;
	      
	      // get the truth pT
	      double tgpT = sqrt(tpx*tpx+tpy*tpy);	  	  
	      
	      // record embedded pions and Hijing tracks separately
	      if(tembed == 0)
		{
		  hpt_hijing_truth->Fill(tgpT);
		}
	      else
		{
		  hpt_truth->Fill(tgpT);
		}
	      
	      if(tembed == 1)
		{
		  double geta = asinh(tpz/sqrt(tpx*tpx+tpy*tpy));
		  hgeta->Fill(geta);
		}
	      
	    }

	  //cout << "n_embed_gtrack = " << n_embed_gtrack << " n_embed_gtrack_1 " << n_embed_gtrack_1 << endl;
	  
	  if(verbose > 0) cout << "Process reco tracks:" << endl;
	  // overestimate number of tracks
	  ntracks = 300;
	  int n_embed_rtrack = 0;                  
	  for(int ir=nr;ir<nr+ntracks;ir++)
	    {
	      
	      int recoget = ntp_track->GetEntry(ir);
	      if(!recoget)
		{
		  cout << "Failed to get ntp_track entry " << ir << endl;
		  //exit(1);
		  //continue;
		  break;
		}
	      
	      if(revent != iev)
		{
		  if(verbose) cout << " change of event, new revent = " << revent << endl;
		  break;
		}
	      
	      nr++;
	      
	      // only embedded tracks
	      if(rgembed != 1)
		continue;
	      
	      if(verbose > 0) cout << " ir = " << ir << " revent " << revent << " rquality " << rquality << "  rgtrackid = " << rgtrackid << " rgnhits " << rgnhits << endl;
	      
	      if(rgembed == 1)
		n_embed_rtrack++;
	      
	      //skip tracks that do not pass through all layers in truth
	      if(rgnhits < nlayers-nmissed)
		{
		  cout << "   --- skip because rgnhits too small " << endl;
		  continue;
		}

	      double rgpT = sqrt(rgpx*rgpx+rgpy*rgpy);	  	  
	      double rpT = sqrt(rpx*rpx+rpy*rpy);
	      
	      if( isnan(rdca2d) )
		rdca2d = 0.0;
	      
	      double corrected_rdca2d = rdca2d;
	      bool calculate_dca2d = true;    // subtract event vertex position from dca values to get dca2d
	      if(calculate_dca2d)
		{
		  corrected_rdca2d = sqrt(pow(rpcax - evx, 2) + pow(rpcay-evy, 2));		  

		  if(rpcax < 0)
		    corrected_rdca2d *= -1.0;
		}
	      
	      /*
		cout << "dca2d " << rdca2d  
		<< " rpcax " << rpcax 
		<< " rpcay " << rpcay
		<< " rpcaz " << rpcaz
		<< endl;
		cout << " corrected rpcax " << rpcax - evx
		<< " rpcay " << rpcay - evy
		<< " rpcaz " << rpcaz - evz
		<< " rdca2d " << corrected_rdca2d 
		<< endl;
	      */
	      
	      // Does the track pass the dca2d cut? 
	      // If so, add to hquality
	      if(rgembed == 1)
		//if(rpurity > nlayers-nmissed && rgpT > 0.5 && fabs(corrected_rdca2d) < dca_cut)
		{
		  hquality->Fill(rquality);
		}
	      
	      // Make a histogram of purity vs quality or dca sigmas with no cuts
	      int ipurity = nlayers - (int) rpurity;
	      //cout << " nlayers " << nlayers << " rpurity " << rpurity << " ipurity " << ipurity << endl;
	      if(rgembed != 1)
		{
		  if(ipurity >= 0 && ipurity < nmissed)
		    {
		      //cout << "ipurity = " << ipurity << endl;
		      hpurity_quality[ipurity]->Fill(rquality);
		      hpurity_dca[ipurity]->Fill(corrected_rdca2d);
		    }
		}
	      
	      // Now make some cuts on dca and quality
	      
	      if(rquality > quality_cut)
		{
		  if(verbose > 0) cout << "   --- failed quality cut - rejected " << endl;
		  continue;
		}
	      
	      // The track passed the quality cut, add it to the 3 panel  dca2d histos
	      // no embedded particles for this
	      if(rgembed != 1)
		if(rquality < quality_cut)
		  {
		    if(rgpT > 0.5 && rgpT <= 1.0)
		      hdca2d[0]->Fill(corrected_rdca2d);
		    else if(rgpT > 1.0 && rgpT <= 2.0)
		      hdca2d[1]->Fill(corrected_rdca2d);
		    else if(rgpT > 2.0)
		      hdca2d[2]->Fill(corrected_rdca2d);
		    
		    if(rgpT > 0.5)
		      hdca2dsigma->Fill(corrected_rdca2d/rdca2dsigma);
		  }
	      
	      if(fabs(corrected_rdca2d) > dca_cut)
		{
		  cout << "   --- skip because failed dca cut" << endl;
		  continue;
		}

	      // Histogram the track Z DCA for Hijing tracks
	      if(rgembed != 1)
		hZdca->Fill(rvz - evz);
	      
	      // drop tracks outside 1 mm in Z dca from the event vertex    
	      if(fabs(rvz - evz) > 0.1)
		{	      
		  cout << "skip because failed z vertex cut" << endl;
		  continue;
		}
	      	      	      
	      //cout << " rgembed = " << rgembed << " rgpT = " << rgpT << endl;
	      // histogram the reco pT for the embedded tracks that pass the cuts
	      if(verbose > 0) cout << "    rnhits = " << rnhits << " rgembed " << rgembed << endl;
	      //if(rgembed == 1 && rnhits > nlayers-nmissed)
		{
		  if(verbose > 0) cout << "   accepted: rgembed = " << rgembed << " rgpT = " << rgpT << endl;
		  naccept++;
		  hpt_compare->Fill(rgpT,rpT/rgpT);
		  hpt_dca2d->Fill(rgpT, corrected_rdca2d);
		  hnhits->Fill(rnhits);
		}
		//else
		//cout << "   --- rejected because of low rnhits" << endl;
	      
	      if(rgembed != 1)
		{
		  // for the non-embedded tracks from Hijing, we want to be able to do pt_sigma cuts later	  
		  hpt_hijing_compare->Fill(rgpT,rpT/rgpT);
		  
		  for(int ipt=0;ipt<NVARBINS-1;ipt++)
		    {
		      if(rpT > xbins[ipt] && rpT < xbins[ipt+1])
			hptreco[ipt]->Fill(rpT);
		    }
		}
	      
	      // make histograms of the purity for the tracks that passed the quality and dca cuts
	      if(ipurity >= 0 && ipurity <= npurity)
		{
		  if(use_reco_pt)
		    hpurity[ipurity]->Fill(rpT);
		  else
		    hpurity[ipurity]->Fill(rgpT);
		}
	      
	      if(rgembed == 1)
		{	  
		  double reta = asinh(rpz/sqrt(rpx*rpx+rpy*rpy));
		  hreta->Fill(reta);
		}
	      
	    }  // end loop over reco'd tracks
	  
	  if(verbose > 0) cout << " Done with loop: n_embed_rtrack = " << n_embed_rtrack << endl; 
     
	  //nr += ntracks;
	  //ng += ngtracks;
	  cout << "  accepted tracks this event = " << naccept << endl;
      
	}  // end loop over events

      delete ntp_gtrack;
      delete ntp_track;
      delete ntp_cluster;
      delete ntp_vertex;
      
    } // end loop over files

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
  //hgeta->SetMinimum(0.0);
  //gPad->SetLogy(1);
  hgeta->Draw();
  hreta->SetLineColor(kRed);
  hreta->Draw("same");

  cout << " hgeta integral = " << hgeta->Integral()  << " hreta integral = "  << hreta->Integral() << endl;

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
      if(i==nlayers-1)
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

  //==============
  // Output the results

  TFile *fout = new TFile("root_files/purity_out.root","recreate");

  hnhits->Write();

  // Three dca2d histos made from Hijing tracks
  hdca2d[0]->Write();
  hdca2d[1]->Write();
  hdca2d[2]->Write();

  // truth pT distributions for embedded and non-embedded particles
  hpt_truth->Write();
  hpt_hijing_truth->Write();
  
  // 2D histos made with embedded pions
  hpt_compare->Write();
  hpt_dca2d->Write();

  // 2d histo made with hijing tracks only
  hpt_hijing_compare->Write();
  
  // efficiency plot made with Hijing tracks
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

  //h_evt_dca2d->Write();

  hreta->Write();
  hgeta->Write();
  
  fout->Close();
  
  
}

