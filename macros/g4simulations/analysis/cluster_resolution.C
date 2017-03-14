
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

void cluster_resolution()
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
  static const int nlayers = 67;  
  static const int nmissed = 32;  // maximum number of missed layers
  double ptmax = 35.0;
  
  //based on plots of purity vs DCA  and vs quality
  //double quality_cut = 1.0;
  //double dca_cut = 0.1; // 1mm, used only if use_dca_sigmas is false
  int a = 0;
  int b = 0;
  int c = 0;
  
  //===============================
  
  // Open the evaluator output file 
  cout << "Reading ntuple " << endl;
      
  TH1D *h1 = new TH1D("h1","",7,0,7); // maps per layer hits
  TH1D *h2 = new TH1D("h2","",7,0,7); // IT per layer hits
  TH1D *h3 = new TH1D("h3","",65,0,65); // TPC per layer hits
  TH1D *clusters_per_layer_per_g4_track = new TH1D("","",70,0,70);
  TH1D *clusters_per_layer_per_reco_track = new TH1D("","",70,0,70);

  TH2D *t2 = new TH2D("t2","",24, -1, 7, 24, -1, 7); 

  TH2D *rphi = new TH2D("rphi","Cluster map",2000, -79.0, 79.0, 2000, -79.0,79.0); 
  rphi->GetYaxis()->SetTitle("Y (cm)");
  rphi->GetXaxis()->SetTitle("X (cm)");

  TH2D *delta_rphi = new TH2D("delta_rphi","cluster errors by layer",70.0, 0.0, 70.0, 2000, -0.05, 0.05); 
  delta_rphi->GetYaxis()->SetTitle("Cluster Error (cm)");
  delta_rphi->GetXaxis()->SetTitle("Tracking Layer");

  TH2D *cluster_size = new TH2D("cluster_size","cluster size by layer",70.0, 0.0, 70.0, 200, 0.0, 15.0); 
  cluster_size->GetYaxis()->SetTitle("Cluster Size (hits)");
  cluster_size->GetXaxis()->SetTitle("Tracking Layer");

  int ng4_tracks = 0;
  int nreco_tracks = 0;
  
  // The condor job output files
  for(int i=0;i <500; i++)
    {
      TChain* ntp_track = new TChain("ntp_track","reco tracks");
      TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
      TChain* ntp_vertex = new TChain("ntp_vertex","events");
      TChain *ntp_cluster = new TChain("ntp_cluster","clusters");
      
      char name[500];
      // latest files 
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/eval_output/g4svx_eval_%i.root",i);

      // cylinders and old TPC
      //sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/nmissing_fixed_ups1s_cylinder_2pcIT_eval_output/g4svx_eval_%i.root",i);

      // ladders and old TPC
      sprintf(name,"/sphenix/user/frawley/QTG_simulations/macros/macros/g4simulations/no_charge_sharing_eval_output/g4svx_eval_%i.root",i);

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

	  // count g4 tracks with track ID 1 or 2
	  for(int j=gtr_start; j < gtr_end; j++)
	    {
	      int tget = ntp_gtrack->GetEntry(j);
	      if(tget == 0)
		cout << "Did not find entry in ntp_gtrack " << j << endl;

	      if( tgtrackid == 1 || tgtrackid == 2 )
		ng4_tracks++;
	    }
	  
	  std::multimap<int,int> gtrackID_map;  // map of gtrackID to trackID
	  std::multimap<int,int> cluster_map;     // map of 
	  
	  //Use loop to fill multimap of (gtrackid, ntp_track entry number)
	  for(int j=tr_start; j < tr_end; j++)
	    {
	      int tget = ntp_track->GetEntry(j);
	      if(tget == 0)
		cout << "Did not find entry in ntp_track " << j << endl;
	      gtrackID_map.insert(std::pair<int,int>(rgtrackid, j));  

	      // count reco tracks with rtrackid 0 or 1 and rgtrackid 1 or 2 (i.e. matched reconstructed tracks
	      if( (rtrackid == 0 && ( rgtrackid == 1 || rgtrackid == 2 ) )  || ( rtrackid == 1 && (rgtrackid == 1) || rgtrackid == 2) )
		nreco_tracks++;
	    }

	  // use loop to fill multimap of (gtrackID,layer number) for all clusters in this event
	  for(int p=0;p < ntp_cluster->GetEntries(); p++)
	    {
	      int tget = ntp_cluster->GetEntry(p);
	      if(tget == 0)
		cout << "Did not find entry in ntp_cluster " << p << endl;
	      
	      if(cevent == iev)
		{
		  // Used later to count layers / track
		  cluster_map.insert(std::pair<int,int>(gtrackID, layer));  
		  
		  rphi->Fill(x, y);

		  // Histogram number of clusters per layer per g4 track 
		  if(gtrackID == 1 || gtrackID == 2)
		    clusters_per_layer_per_g4_track->Fill(layer);

		  // Check that this cluster is matched to one of the electron tracks
		  //      This is written to be run on ntuples containing dielectron decays of single Upsilon events - i.e. two tracks - so we can hard-code trackID values
		  //      The following track selections will need to be modified for running on files obtained from events with more than 2 tracks
		  //========================================================================================
		  if( (trackID == 0 && ( gtrackID == 1 || gtrackID == 2 ) )  || ( trackID == 1 && (gtrackID == 1) || gtrackID == 2) )
		    {

		      // extract the cluster r-phi resolution 
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

		      if(layer < 3)
			{
			  //if(size > 1)  // optional cut on hits/cluster for MAPS
			  delta_rphi->Fill( (double) layer, sign * drphi); 
			}
		      else
			delta_rphi->Fill( (double) layer, sign * drphi); 
		      
		      // Extract the number of hits per cluster

		      cluster_size->Fill( (double) layer, size);

		      // histogram the number of clusters per reconstructed track vs layer number for these reco tracks
		      clusters_per_layer_per_reco_track->Fill(layer);
		    }
		}
	    }
	  
	  // Capture the number of hit layers per track for each sysbsystem
	  //=============================================

 	  std::vector<int> vtrack;
	  std::vector<pair<int,double> > ptrack;
	  
	  // loop over truth tracks in ntp_gtrack
	  for(int k = gtr_start; k < gtr_end;k++) 
	    {
	      if( !(gtrackID == 1 || gtrackID == 2) )
		 continue;
	      
	      int tget = ntp_gtrack->GetEntry(k);  
	      if(tget == 0)
		cout << "Did not find entry in ntp_gtrack " << k << endl;
	      
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
		  if((*it).second > 2 && (*it).second < 7)
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
	      
	      h1->Fill(maps_counter); // 1D histogram with maps hits 
	      h2->Fill(IT_counter); // 1D histogram with IT layers hits
	      h3->Fill(tpc_counter); 

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


  TCanvas *C1 = new TCanvas("C1","C1",50,50,1000,700);
  C1->Divide(3,1);
  C1->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.01);
  h1->SetTitle("Maps layers / track");
  //h1->SetTitleSize(0.08);
  h1->GetXaxis()->SetTitle("Layers/track");
  h1->GetXaxis()->SetTitleSize(0.055);
  h1->GetXaxis()->SetLabelSize(0.055);
  h1->GetYaxis()->SetLabelSize(0.055);
  h1->Draw();

  C1->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.01);
  h2->SetTitle("INTT layers / track");  
  h2->GetXaxis()->SetTitle("Layers/track");
  h2->GetXaxis()->SetTitleSize(0.055);
  h2->GetXaxis()->SetLabelSize(0.055);
  h2->GetYaxis()->SetLabelSize(0.055);
  h2->Draw();

  C1->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.01);
  h3->SetTitle("TPC layers / track");  
  h3->GetXaxis()->SetTitle("Layers/track");
  h3->GetXaxis()->SetTitleSize(0.055);
  h3->GetXaxis()->SetLabelSize(0.055);
  h3->GetYaxis()->SetLabelSize(0.055);
  h3->Draw();

  C1->Print("Layers_per_track.pdf","pdf");

  TCanvas *C5 = new TCanvas("C5","C5",50,50,600,600); 
  rphi->Draw("p");

  TCanvas *C6 = new TCanvas("C6","C6",50,50,800,800); 
  C6->SetLeftMargin(0.15);
  delta_rphi->GetYaxis()->SetTitleOffset(1.7);
  delta_rphi->Draw("p");

  TCanvas *c7 = new TCanvas("c7","c7",50,50,1200,800); 
  c7->Divide(3,1);

  c7->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.01);
  TH1D *hpy1 = new TH1D("hpy1","MAPS clusters",500, -0.05, 0.05);
  delta_rphi->ProjectionY("hpy1",1,3);
  hpy1->GetXaxis()->SetRangeUser(-0.0016, 0.0016);
  hpy1->GetXaxis()->SetTitle("cluster error (cm)");
  hpy1->SetTitleOffset(0.1,"X");
  hpy1->GetXaxis()->SetTitleSize(0.05);
  hpy1->GetXaxis()->SetLabelSize(0.06);
  hpy1->GetYaxis()->SetLabelSize(0.06);
  hpy1->GetXaxis()->SetNdivisions(506);
  hpy1->GetXaxis()->SetTitleOffset(1.1);
  hpy1->Draw();
  double rms1 = 10000.0 * hpy1->GetRMS();
  char label[500];
  sprintf(label,"RMS %.1f #mu m",rms1);
  TLatex *l1 = new TLatex(0.55,0.92,label);
  l1->SetNDC(1);
  l1->Draw();

  c7->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.01);
  TH1D *hpy2 = new TH1D("hpy2","INTT clusters",500, -0.05, 0.05);
  delta_rphi->ProjectionY("hpy2",4,7); // for 2 layers
  hpy2->GetXaxis()->SetRangeUser(-0.011, 0.011);
  hpy2->GetXaxis()->SetTitle("cluster error (cm)");
  hpy2->GetXaxis()->SetTitleOffset(0.6);
  hpy2->GetXaxis()->SetTitleSize(0.05);
  hpy2->GetXaxis()->SetLabelSize(0.06);
  hpy2->GetYaxis()->SetLabelSize(0.06);
  hpy2->GetXaxis()->SetNdivisions(506);
  hpy2->GetXaxis()->SetTitleOffset(1.1);
  hpy2->Draw();
  double rms2 = 10000 * hpy2->GetRMS();
  sprintf(label,"RMS %.1f #mu m",rms2);
  TLatex *l2 = new TLatex(0.55,0.92,label);
  l2->SetNDC(1);
  l2->Draw();

  c7->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.01);
  TH1D *hpy3 = new TH1D("hpy3","TPC clusters",500,-0.05, 0.05);
  delta_rphi->ProjectionY("hpy3",8,40);
  hpy3->GetXaxis()->SetRangeUser(-0.045, 0.045);
  hpy3->GetXaxis()->SetNdivisions(506);
  hpy3->GetXaxis()->SetTitle("cluster error (cm)");
  hpy3->GetXaxis()->SetTitleOffset(1.1);
  hpy3->GetXaxis()->SetTitleSize(0.05);
  hpy3->GetXaxis()->SetLabelSize(0.06);
  hpy3->GetYaxis()->SetLabelSize(0.06);
  hpy3->Draw();
  double rms3 = 10000 * hpy3->GetRMS();
  sprintf(label,"RMS %.1f #mu m",rms3);
  TLatex *l3 = new TLatex(0.55,0.92,label);
  l3->SetNDC(1);
  l3->Draw();

  c7->Print("Cluster_errors.pdf","pdf");
  
  cout << "MAPS cluster number " << hpy1->Integral() << " RMS = " << 10000 * hpy1->GetRMS() << " microns" << endl;
  cout << "INTT cluster number " << hpy2->Integral() << " RMS = " << 10000 * hpy2->GetRMS() << " microns" << endl;
  cout << "TPC cluster number " << hpy3->Integral() << " RMS = " << 10000 * hpy3->GetRMS() << " microns" << endl;

  TCanvas *C8 = new TCanvas("C8","C8",50,50,800,800); 
  C6->SetLeftMargin(0.15);
  delta_rphi->GetYaxis()->SetTitleOffset(1.7);
  cluster_size->Draw("p");

  TCanvas *c9 = new TCanvas("c9","c9",50,50,1200,800); 
  c9->Divide(3,1);

  c9->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetLogy(1);
  TH1D *hpc1 = new TH1D("hpc1","MAPS hits/cluster",200, 0, 15.0);
  cluster_size->ProjectionY("hpc1",1,3);
  //hpc1->GetXaxis()->SetRangeUser(0.0, 15.0);
  hpc1->GetXaxis()->SetTitle("cluster size (pixels)");
  hpc1->SetTitleOffset(0.1,"X");
  hpc1->GetXaxis()->SetTitleSize(0.05);
  hpc1->GetXaxis()->SetLabelSize(0.06);
  hpc1->GetYaxis()->SetLabelSize(0.06);
  hpc1->GetXaxis()->SetNdivisions(506);
  hpc1->GetXaxis()->SetTitleOffset(1.1);
  hpc1->Draw();
  double rmsc1 = hpc1->GetMean();
  sprintf(label,"Mean %.1f  (hits)",rmsc1);
  TLatex *lc1 = new TLatex(0.55,0.92,label);
  lc1->SetNDC(1);
  lc1->Draw();

  c9->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetLogy(1);
  TH1D *hpc2 = new TH1D("hpc2","INTT hits/cluster",200, 0, 10.0);
  cluster_size->ProjectionY("hpc2",4,7); // for 2 layers
  hpc2->GetXaxis()->SetRangeUser(0.0, 10.0);
  hpc2->GetXaxis()->SetTitle("cluster size (hits)");
  hpc2->GetXaxis()->SetTitleOffset(0.6);
  hpc2->GetXaxis()->SetTitleSize(0.05);
  hpc2->GetXaxis()->SetLabelSize(0.06);
  hpc2->GetYaxis()->SetLabelSize(0.06);
  hpc2->GetXaxis()->SetNdivisions(506);
  hpc2->GetXaxis()->SetTitleOffset(1.1);
  hpc2->Draw();
  double rmsc2 = hpc2->GetMean();
  sprintf(label,"Mean %.1f (hits)",rmsc2);
  TLatex *lc2 = new TLatex(0.55,0.92,label);
  lc2->SetNDC(1);
  lc2->Draw();

  c9->cd(3);
  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.01);
  gPad->SetLogy(1);
  TH1D *hpc3 = new TH1D("hpc3","TPC hits/cluster",200, 0, 10.0);
  cluster_size->ProjectionY("hpc3",8,40);
  hpc3->GetXaxis()->SetRangeUser(0.0, 5.0);
  hpc3->GetXaxis()->SetNdivisions(506);
  hpc3->GetXaxis()->SetTitle("cluster size (hits)");
  hpc3->GetXaxis()->SetTitleOffset(1.1);
  hpc3->GetXaxis()->SetTitleSize(0.05);
  hpc3->GetXaxis()->SetLabelSize(0.06);
  hpc3->GetYaxis()->SetLabelSize(0.06);
  hpc3->Draw();
  double rmsc3 = hpc3->GetMean();
  sprintf(label,"Mean %.1f (hits)",rmsc3);
  TLatex *lc3 = new TLatex(0.55,0.92,label);
  lc3->SetNDC(1);
  lc3->Draw();

  c9->Print("Hits_per_cluster.pdf","pdf");

  /*
  TCanvas *c10 = new TCanvas("c10","c10",50,50,1200,800); 
  c10->Divide(2,1);

  c10->cd(1);
  double gnorm = 1.0 / (double) ng4_tracks;
  clusters_per_layer_per_g4_track->Scale(gnorm);
  clusters_per_layer_per_g4_track->GetYaxis()->SetTitle("Clusters per g4 track");
  clusters_per_layer_per_g4_track->GetYaxis()->SetTitleOffset(1.2);
  clusters_per_layer_per_g4_track->GetXaxis()->SetTitle("Layer");
  clusters_per_layer_per_g4_track->SetMaximum(1.2);
  clusters_per_layer_per_g4_track->Draw();

  c10->cd(2);
  //char temp[500];
  double norm = 1.0 / (double) nreco_tracks;
  //sprintf(temp,"%.2f",norm);
  cout << "nreco_tracks " << nreco_tracks << " ng4_tracks = " << ng4_tracks  << " norm " << norm << endl;
  //TF1 *m = new TF1(temp);
  //clusters_per_layer_per_reco_track->Multiply(m);
  clusters_per_layer_per_reco_track->Scale(norm);
  clusters_per_layer_per_reco_track->GetYaxis()->SetTitle("Clusters per reco track");
  clusters_per_layer_per_reco_track->GetYaxis()->SetTitleOffset(1.2);
  clusters_per_layer_per_reco_track->GetXaxis()->SetTitle("Layer");
  clusters_per_layer_per_reco_track->SetMaximum(1.2);
  clusters_per_layer_per_reco_track->Draw();
  */

  TCanvas *c11 = new TCanvas("c11","c11",5,5,600,800);
  
  hpy1->Draw();
  c11->SetLeftMargin(0.15);  
  c11->SetBottomMargin(0.12);  
  //rms3 = 10000 * hpy1->GetRMS();
  //sprintf(label,"RMS %.1f #mu m",rms3);
  //TLatex *l3 = new TLatex(0.55,0.92,label);
  //l3->SetNDC(1);
  l1->Draw();


  TCanvas *c12 = new TCanvas("c12","c12",500,5,600,800);
  c12->SetLeftMargin(0.15);  
  c12->SetBottomMargin(0.12);  
  hpc1->Draw();
  //rms3 = 10000 * hpy1->GetRMS();
  //sprintf(label,"RMS %.1f #mu m",rms3);
  //TLatex *l3 = new TLatex(0.55,0.92,label);
  //l3->SetNDC(1);
  lc1->Draw();


}
  	     
