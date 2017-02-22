// was working but then wanted to change to multimap to make faster


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
#include <map>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <string.h>


using namespace std;

void truth_match_v2()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);//(0)
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1); //(0)
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

  bool verbose = true;
  
  // Setup parameters
  //=================
  
  int n_maps_layers = 3;
  int n_intt_layers = 4;
  
  // This should be inner maps layers (3) + intt(4) + TPC (60)
  static const int nlayers = 67;  // maximum number of tracking layers for IT
  static const int nmissed = 32;  // maximum number of missed layers
  double ptmax = 35.0;
  
  //based on plots of purity vs DCA  and vs quality
  double quality_cut = 1.0;
  double dca_cut = 0.1; // 1mm, used only if use_dca_sigmas is false
  int a = 0;
  int b = 0;
  int c = 0;
  
  //===============================
  
  // Open the evaluator output file 
  cout << "Reading ntuple " << endl;
  
  TCanvas *C1 = new TCanvas("C1","C1",50,50,800,600);
  
  TH1D *h1 = new TH1D("h1","",7,0,7); // maps per layer hits
  TH1D *h1pt = new TH1D("h1pt","",80,0.0,40.0); // maps per layer pt
  TH1D *h2 = new TH1D("h2","",7,0,7); // IT per layer hits
  TH1D *h3 = new TH1D("h3","",40,0,40); // TPC per layer hits
  TH1D *h2pt = new TH1D("h2pt","",80,0,40); // IT per layer pt
  TH2D *t2 = new TH2D("t2","",24, -1, 7, 24, -1, 7); 
  //TH2D *rphi = new TH2D("rphi","MAPS+INTT clusters",1000, -14.0, 14.0, 1000, -14.0,14.0); 
  TH2D *rphi = new TH2D("rphi","MAPS+INTT clusters",2000, -79.0, 79.0, 2000, -79.0,79.0); 
  rphi->GetYaxis()->SetTitle("Y");
  rphi->GetXaxis()->SetTitle("X");
  TH1D *t2x = new TH1D("t2x","",8,-1,7); //for  projection of t2 onto x-axis
  TH1D *t2y = new TH1D("t2y","",8,-1,7); // look at map hits for whichever layer psecified here
  // for cuts
  TH1D *t2xx = new TH1D("t2xx","",7,-2.5,5.5); // look at IT hits for whichever layer specified here
  TH1D *t2yy = new TH1D("t2yy","",7,-2.5,5.5); // look at maps hits for whichever layer specified here

  TH2D *delta_rphi = new TH2D("delta_rphi","MAPS+INTT clusters",70.0, 0.0, 70.0, 2000, -0.05, 0.05); 
  
  // The condor job output files
  for(int i=0;i <50; i++)
    {
      TChain* ntp_track = new TChain("ntp_track","reco tracks");
      TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
      TChain* ntp_vertex = new TChain("ntp_vertex","events");
      TChain *ntp_cluster = new TChain("ntp_cluster","clusters");
      
      char name[500];
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/eval_output/g4svx_eval_%i.root",i);

      // cylinder cell + new TPC
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/Sourav_macro_A_eval_output/g4svx_eval_%i.root",i);
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/Sourav_macro_B_eval_output/g4svx_eval_%i.root",i);
 
      // ladders + new TPC
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/ups1s_lad_improved_tpc_eval_output/g4svx_eval_%i.root",i);

      // cylinder cell silicon + old TPC
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/cylinder_silicon_massfix_intt4_cbfit_eval_output/g4svx_eval_%i.root",i);
      sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/cylinders+old_TPC_refit_eval_output/g4svx_eval_%i.root",i);

      ntp_vertex->Add(name);
      ntp_track->Add(name);
      ntp_gtrack->Add(name);
      ntp_cluster->Add(name);
      
      // This include file contains the definitions of the ntuple variables
#include "ntuple_variables.C"
      
      int ntr =  ntp_track->GetEntries();
      int ngtr = ntp_gtrack->GetEntries();
      std:: cout << "File " << i  << " ntracks entries = " << ntr << " ngtracks entries = " << ngtr << endl;
      
      int gtr_start = 0;
      int gtr_end = 0;
      int tr_start = 0;
      int tr_end = 0;
      
      // Loop over events
      int nevt = ntp_vertex->GetEntries();
      cout << "Number of events in this file = " << nevt << endl;
      
      for(int iev=0;iev<nevt;iev++)
	{
	  ntracks = 0.0;
	  ngtracks = 0.0;
	  
	  int evtget = ntp_vertex->GetEntry(iev);
	  if(!evtget)
	    cout << "Failed to get event " << iev << endl;
	  
	  int ntr_ev = (int) ntracks;
	  int ngtr_ev = (int) ngtracks;
	  
	  tr_end = tr_start + ntr_ev;
	  gtr_end = gtr_start + ngtr_ev;
	  
	  std::multimap<int,int> gtrackID_map;
	  std::multimap<int,int> cluster_map;
	  //std::multimap<int,int> cluster_map_tpc;
	  
	  //Use loop to fill multimap of (gtrackid, ntp_track entry number)
	  for(int j=tr_start; j < tr_end; j++)
	    {
	      int tget = ntp_track->GetEntry(j);
	      if(tget == 0)
		cout << "Did not find entry in ntp_track " << j << endl;
	      /*
		if((rgembed == 0) 
		continue;
	      */
	      //cout << "     adding to gtrackID_map:   rgtrackid: " << rgtrackid << ", entry: " << j << endl;
	      gtrackID_map.insert(std::pair<int,int>(rgtrackid, j));  
	    }
	  
	  // use loop to fill multimap of (gtrackID,layer number) for all layers up to 6
	  //cout << " cluster entries " << ntp_cluster->GetEntries() << endl;
	  
	  for(int p=0;p < ntp_cluster->GetEntries(); p++)
	    {
	      int tget = ntp_cluster->GetEntry(p);
	      if(tget == 0)
		cout << "Did not find entry in ntp_cluster " << p << endl;
	      
	      /*
		if(cgembed == 0) 
		continue;
	      */
	      if(cevent == iev)
		{
		  if(layer < n_maps_layers + n_intt_layers )
		    cluster_map.insert(std::pair<int,int>(gtrackID, layer));  
		  
		  if(layer > n_maps_layers + n_intt_layers - 1)
		    cluster_map.insert(std::pair<int,int>(gtrackID, layer));  
		  
		  rphi->Fill(x, y);

		  // Check that this cluster is matched to one of the electron tracks
		  //      This is written to be run on dielectron decays of single Upsilon events - i.e. two tracks - so we can hard-code trackID values
		  //      The following track selections will need to be modified for running on files obtained from other kinds of events
		  if( (trackID == 0 && ( gtrackID == 1 || gtrackID == 2 ) )  || ( trackID == 1 && (gtrackID == 1) || gtrackID == 2) )
		    {
		      // extract the cluseter r-phi resolution 
		      double dx = x - gx;
		      double dy = y - gy;
		      double drphi = sqrt(dx*dx + dy*dy);
		      // Get the sign of the difference in r-phi space
		      double dphi = atan(dy/dx);
		      double phi = atan(gy/gx);
		      
		      // is the measured hit CW or CCW from the truth? Give it a sign
		      double sign = 0.0;
		      if(atan(y/x) > atan(gy/gx)) 
			sign = +1.0;
		      else
			sign = -1.0;
		      // are these matched Upsilon decay tracks?
		      
		      delta_rphi->Fill( (double) layer, sign * drphi); 
		    }
		}
	    }
	  
 	  std::vector<int> vtrack;
	  std::vector<pair<int,double> > ptrack;
	  
	  // loop over truth tracks in ntp_gtrack
	  // 
	  for(int k = gtr_start; k < gtr_end;k++) 
	    {
	      
	      int tget = ntp_gtrack->GetEntry(k);  
	      if(tget == 0)
		cout << "Did not find entry in ntp_gtrack " << k << endl;
	      
	      /*
		if(!tembed)
		continue;
	      */	
	      
	      double gpt = sqrt(pow(tpx,2) + pow(tpy,2));
	      
	      // declare pair of multimap iterators named "range"
	      std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> range;
	      range = cluster_map.equal_range(tgtrackid); // satisfies if(rgtrackid == tgtrackid)
	      
	      int maps_counter = 0; // reset the counter each time a new trackID comes through the loop to count number of times it repeats in multimap
	      int IT_counter = 0;
	      int tpc_counter = 0;
	      
	      // count number of times a key is repeated in cluster multimap
	      for(std::multimap<int,int>::iterator it = range.first; it != range.second; it++) 
		{
		  if((*it).second < 3) // layers 0,1,2
		    {
		      maps_counter++;
		      a += 1;
		    }
		  if((*it).second > 2)
		    {
		      IT_counter++; // layers 3,4,5,6
		      b += 1;
		    }
		  if((*it).second > 6)
		    {
		      tpc_counter++; // layers > 6
		      c += 1;
		    }
		}
	      t2->Fill((double)maps_counter,(double)IT_counter); // 2D histogram with maps 
	      //t2->Fill(maps_counter,IT_counter); // same as above statement 
	      
	      h1->Fill(maps_counter); // 1D histogram with maps hits 
	      h2->Fill(IT_counter); // 1D histogram with IT layers hits
	      h3->Fill(tpc_counter); 

	      if(maps_counter == 3)
		h1pt->Fill(gpt); // maps momentum
	      if(IT_counter == 4)
		h2pt->Fill(gpt); // IT layers momentum
	      
	    }// for k

	  tr_start += ntr_ev;
	  gtr_start += ngtr_ev;

	} // for iev

      delete ntp_gtrack;
      delete ntp_track;
      delete ntp_cluster;
      delete ntp_vertex;
      
    }// for i
  
  cout << "Total hits maps: " << a << ", Total hits IT: " << b << endl;
  C1->Divide(1,4);
  C1->cd(1);
  h1->SetTitle("Tracks that hit maps layers");
  h1->Draw();
  C1->cd(2);
  h2->SetTitle("Tracks that hit IT layers");  
  h2->Draw();
  C1->cd(3);
  h1pt->SetTitle("Momentum of maps layers"); 
  h1pt->Draw();
  C1->cd(4);
  h2pt->SetTitle("Momentum of IT layers"); 
  h2pt->Draw();

  TCanvas *C2 = new TCanvas("C2","C2",50,50,800,600);
  t2->SetTitle("Layer hits per track INTT vs. MAPS");
  t2->SetMarkerStyle(20);
  t2->SetMarkerSize(0.5);
  t2->GetXaxis()->SetTitle("Maps layers (3)");
  t2->GetYaxis()->SetTitle("IT layers (4)");
  t2->Draw();

  TCanvas *C3 = new TCanvas("C3","C3",50,50,800,800); 
  C3->Divide(1,2);
  C3->cd(1);
  t2->ProjectionX("t2x");
  t2x->SetTitle("MAPS layers hit per track");
  t2x->Draw();
  C3->cd(2);
  t2->ProjectionY("t2y");
  t2y->SetTitle("INTT layers hit per track");
  t2y->Draw();

  /*
  //cuts on 2D
  TCanvas *C4 = new TCanvas("C4","C4",50,50,600,600); 
  C4->Divide(1,2);
  C4->cd(1);
  t2->ProjectionX("t2x");
  t2xx->SetTitle("3 maps hits with 4 IT hits");
  t2xx->SetAxisRange(2.5, 3.5); // look at maps hits for whichever layer is specified here
  t2xx->Draw();
  C4->cd(2);
  t2->ProjectionY("t2y");
  t2yy->SetTitle("3 maps hits with 0 IT hits");
  t2yy->SetAxisRange(-0.5, 0.5); // look at IT hits for whichever layer is specified here
  t2yy->Draw();
  */

  TCanvas *C5 = new TCanvas("C5","C5",50,50,600,600); 
  rphi->Draw("p");

  TCanvas *C6 = new TCanvas("C6","C6",50,50,800,800); 
  delta_rphi->Draw("p");

  TFile *t1 = new TFile("t1","RECREATE");
  h1->Write();
  h2->Write();
  h1pt->Write();
  h2pt->Write();
  t2->Write();
  t2x->Write();
  t2y->Write();
  rphi->Write();

  TCanvas *c7 = new TCanvas("c7","c7",50,50,1200,800); 
  c7->Divide(3,1);

  c7->cd(1);
  TH1D *hpy1 = new TH1D("hpy1","MAPS",500, -0.05, 0.05);
  delta_rphi->ProjectionY("hpy1",1,3);
  hpy1->GetXaxis()->SetRangeUser(-0.0016, 0.0016);
  hpy1->GetXaxis()->SetTitle("cm");
  hpy1->GetXaxis()->SetNdivisions(506);
  hpy1->Draw();
  double rms1 = 10000.0 * hpy1->GetRMS();
  char label[500];
  sprintf(label,"RMS %.1f #mu m",rms1);
  TLatex *l1 = new TLatex(0.55,0.92,label);
  l1->SetNDC(1);
  l1->Draw();

  c7->cd(2);
  TH1D *hpy2 = new TH1D("hpy2","INTT",500, -0.05, 0.05);
  delta_rphi->ProjectionY("hpy2",4,7); // for 2 layers
  hpy2->GetXaxis()->SetRangeUser(-0.011, 0.011);
  hpy2->GetXaxis()->SetTitle("cm");
  hpy2->GetXaxis()->SetNdivisions(506);
  hpy2->Draw();
  double rms2 = 10000 * hpy2->GetRMS();
  sprintf(label,"RMS %.1f #mu m",rms2);
  TLatex *l2 = new TLatex(0.55,0.92,label);
  l2->SetNDC(1);
  l2->Draw();

  c7->cd(3);
  TH1D *hpy3 = new TH1D("hpy3","TPC",500,-0.05, 0.05);
  delta_rphi->ProjectionY("hpy3",8,40);
  hpy3->GetXaxis()->SetRangeUser(-0.045, 0.045);
  hpy3->GetXaxis()->SetTitle("cm");
  hpy3->GetXaxis()->SetNdivisions(506);
  hpy3->Draw();
  double rms3 = 10000 * hpy3->GetRMS();
  sprintf(label,"RMS %.1f #mu m",rms3);
  TLatex *l3 = new TLatex(0.55,0.92,label);
  l3->SetNDC(1);
  l3->Draw();
  
  cout << "MAPS cluster number " << hpy1->Integral() << " RMS = " << 10000 * hpy1->GetRMS() << " microns" << endl;
  cout << "INTT cluster number " << hpy2->Integral() << " RMS = " << 10000 * hpy2->GetRMS() << " microns" << endl;
  cout << "TPC cluster number " << hpy3->Integral() << " RMS = " << 10000 * hpy3->GetRMS() << " microns" << endl;

}
  	     
