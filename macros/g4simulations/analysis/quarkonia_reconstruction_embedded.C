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

void quarkonia_reconstruction_embedded()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);

  bool verbose = false;
  
  double quality_cut = 3.0;
  double dca_cut = 0.1;

  char lepton[100];
  sprintf(lepton,"electron");
  
  double decaymass=0.000511;
  cout << "Assuming decay particle mass is " << decaymass << endl;
  
  // Open the g4 evaluator output file
  
#ifdef TEST
  bool ups1s = true;
  bool ups2s = false;
  bool ups3s = false;

  cout << "Reading test ntuple " << endl;
  TChain* ntp_track = new TChain("ntp_track","reco tracks");
  ntp_track->Add("../g4svtx_eval.root");

  ntp_track->Print();
  TChain* ntp_gtrack = new TChain("ntp_gtrack","g4 tracks");
  ntp_gtrack->Add("../g4svtx_eval.root");
  ntp_gtrack->Print();

  TChain* ntp_vertex = new TChain("ntp_vertex","events");
  ntp_vertex->Add("../g4svtx_eval.root");
  ntp_gtrack->Print();

  TChain *ntp_cluster = new TChain("ntp_cluster","clusters");
#endif

#ifndef TEST
  cout << "Reading electron ntuples " << endl; 

  bool ups1s = true;
  bool ups2s = false;
  bool ups3s = false;
  
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

  // # include "ntuple_files.C"
  //#include "ntuple_files_maps_5layer.C"
    
#endif
  
  // Ntuple access variables
  // This include file contains the definitions of the ntuple variables                                                                        
#include "ntuple_variables.C"
  
  // Define some histograms
  
  int nbpt = 20;
  double ptmax = 10.0;

  TH1F* hrquality = new TH1F("hrquality", "Reconstructed track Quality", 1000, 0, 10);
  TH1F* hrdca2d = new TH1F("hrdca2d", "Reconstructed track dca2d", 1000, 0, 0.05);
  TH1F* hrpt = new TH1F("hrpt"," pT", nbpt, 0.0, ptmax);
  TH1F* hgpt = new TH1F("hgpt","g4 pT", nbpt, 0.0, ptmax);

  /*  
  int nby = 50;
  int nbphi = 50;
  int nbeta = 50;
  int nbth = 50;

  TH1F* hrcharge = new TH1F("hrcharge","reconstructed charge", 300, 0, 30);
  TH1F* hrpx = new TH1F("hrpx","reconstructed px; p_{x}", 300, 0, 30);
  TH1F* hrpy = new TH1F("hrpy","reconstructed px; p_{y}", 300, 0, 30);
  TH1F* hrpz = new TH1F("hrpz","reconstructed px; p_{z}", 300, 0, 30);
  
  TH1F* hrphi = new TH1F("hrphi"," phi", nbphi, -3.5,3.5);
  TH1F* hry = new TH1F("hry"," y", nby, -2.5, 2.5);
  TH1F* hreta = new TH1F("hreta"," reta", nbeta, -3.0, 3.0);
  TH1F* hrdeta = new TH1F("hrdeta"," rdeta", nbeta, -3.0, 3.0);
  TH1F* hrdphi = new TH1F("hrdphi"," rdphi", nbphi, -3.5, 3.5);
  TH1F* hrdth = new TH1F("hrdth"," rdphi", nbth, 0.0, 3.5);
  TH1F* hrdpt = new TH1F("hrdpt"," rdpt", nbpt, 0.0, 10.0);
  
  TH1F* hrvz = new TH1F("hrvz","rvz", 300, -12.0, 12.0);

  // Define some histograms
  TH1F* htpx = new TH1F("htpx","g4 px; p_{x}", 300, 0, 30);
  TH1F* htpy = new TH1F("htpy","g4 px; p_{y}", 300, 0, 30);
  TH1F* htpz = new TH1F("htpz","g4 px; p_{z}", 300, 0, 30);

  TH1F* hgphi = new TH1F("hgphi","g4 phi", nbphi, -3.5,3.5);
  TH1F* hgy = new TH1F("hgy","g4 y", nby, -2.5, 2.5);
  TH1F* hgeta = new TH1F("hgeta"," geta", nbeta, -3.0, 3.0);
  TH1F* hgdeta = new TH1F("hgdeta"," gdeta", nbeta, -3.0, 3.0);
  TH1F* hgdphi = new TH1F("hgdphi"," gdphi", nbphi, -3.5, 3.5);
  TH1F* hgdth = new TH1F("hgdth"," gdphi", nbth, 0.0, 3.5);
  TH1F* hgdpt = new TH1F("hgdpt"," gdpt", nbpt, 0.0, 10.0);

  TH1F* hgvz = new TH1F("hgvz","gvz", 300, -12.0, 12.0);

  TH2D *hreta12 = new TH2D("hreta12","hreta12",100,-2.2,2.2,100,-4,4);
  TH2D *hrtheta12 = new TH2D("hrtheta12","hrtheta12",100,0.5,2.7,100,0.5,2.7);
  TH2D *hrphi12 = new TH2D("hrphi12","hrphi12",100,-4,4,100,-4,4);
  TH1D *g4radmass = new TH1D("g4radmass","G4 invariant mass after tracking",100,7.0,11.0);
  g4radmass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
  */

  TH1D *g4mass = new TH1D("g4mass","G4 input invariant mass",100,7.0,11.0);
  g4mass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
  TH1D *g4mass_primary = new TH1D("g4mass_primary","G4 input invariant mass",100,7.0,11.0);
  g4mass_primary->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");

  TH1D *recomass = new TH1D("recomass","Reconstructed invariant mass",200,7.0,11.0);
  recomass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
  TH1D *recomass_primary = new TH1D("recomass_primary","Reconstructed invariant mass",200,7.0,11.0);
  recomass_primary->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");

  int nups_requested = 20000;
  cout << "Upsilons requested = " << nups_requested << endl;

  //=======================
  // Loop over events
  //=======================

  int nr = 0;
  int ng = 0;
  int nrecog4mass = 0;
  int nrecormass = 0;

  for(int iev=0;iev<ntp_vertex->GetEntries();iev++)
    {
      // drop out when the requested number of reco'd Upsilons has been reached
      if(nrecormass > nups_requested)
	break;

      int recoget = ntp_vertex->GetEntry(iev);


      if(verbose)
	cout << "iev " << iev
	     << " event " << event
	     << " ntracks  (reco) " << ntracks
	     << " ngtracks (g4) " << ngtracks
	     << " gvz " << egvz
	     << " vz " << evz
	     << endl;
      
      //============================
      // process G4 tracks
      // for this event
      //============================

      int ng4trevt_elec = -1;
      int ng4trevt_pos = -1;
      int g4trnum_elec[1000];
      int g4trnum_pos[1000];

      for(int ig=ng;ig<ng+ngtracks;ig++)
        {
          int recoget1 = ntp_gtrack->GetEntry(ig);
      

	  // we want only electrons or positrons
	  if(tflavor != 11 && tflavor != -11)
	    continue;
	  
	  if(tflavor == 11)
	    {
	      // electron
	      ng4trevt_elec++;
	      if(ng4trevt_elec > 999)
		continue;

	      if(verbose)
		cout << " Found electron:" << endl
		     << "  ig " << ig
		     << " ng4trevt_elec " << ng4trevt_elec
		     << " gtrackID " << tgtrackid
		     << " gflavor " << tflavor
		     << " tpx " << tpx
		     << " tpy " << tpy
		     << " tpz " << tpz
		     << endl;
	      
	      g4trnum_elec[ng4trevt_elec] = ig;
	    }
	  else
	    {
	      // positron
	      ng4trevt_pos++;
	      if(ng4trevt_pos > 999)
		continue;

	      if(verbose)
		cout << " Found positron:" << endl
		     << "  ig " << ig
		     << " ng4trevt_pos " << ng4trevt_pos
		     << " gtrackID " << tgtrackid
		     << " gflavor " << tflavor
		     << " tpx " << tpx
		     << " tpy " << tpy
		     << " tpz " << tpz
		     << endl;
	      
	      g4trnum_pos[ng4trevt_pos] = ig;
	    }
	}	  
      ng4trevt_elec++;
      ng4trevt_pos++;

      if(verbose)
	cout << "For this event found " << ng4trevt_elec << " g4 electrons and " << ng4trevt_pos << " g4 positrons"  << endl;
 
      // make all pairs of g4 electrons and positrons
      for(int ielec=0;ielec<ng4trevt_elec;ielec++)
	{
	  int recoelec = ntp_gtrack->GetEntry(g4trnum_elec[ielec]);
	  
	  double elec_pT = sqrt(tpx*tpx+tpy*tpy);
	  double elec_eta = asinh(tpz/sqrt(tpx*tpx+tpy*tpy));

	  int gtrid1 = tgtrackid;
 
	  TLorentzVector t1;
	  double E1 = sqrt( pow(tpx,2) + pow(tpy,2) + pow(tpz,2) + pow(decaymass,2));	  
	  t1.SetPxPyPzE(tpx,tpy,tpz,E1);	  
	  
	  // print out track details
	  if(verbose)
	    cout << "  Pair electron:  iev " << iev << " ielec " << ielec
		 << " g4trnum_elec " << g4trnum_elec[ielec]
		 << " tgtrackid " << tgtrackid
		 << " tflavor " << tflavor
		 << " tpx " << tpx
		 << " tpy " << tpy
		 << " tpz " << tpz
		 << " elec_eta " << elec_eta
		 << " elec_gpT " << elec_pT
		 << endl;
	  
	  for(int ipos =0;ipos<ng4trevt_pos;ipos++)
	    {
	      int recopos = ntp_gtrack->GetEntry(g4trnum_pos[ipos]);

	      int gtrid2 = tgtrackid;

	      double pos_pT = sqrt(tpx*tpx+tpy*tpy);
	      double pos_eta = asinh(tpz/sqrt(tpx*tpx+tpy*tpy));
	      
	      // print out track details
	      if(verbose)
		cout << "  Pair positron: iev " << iev << " ipos " << ipos
		     << " g4trnum_pos " << g4trnum_pos[ipos]
		     << " tgtrackid " << tgtrackid
		     << " tflavor " << tflavor
		     << " tpx " << tpx
		     << " tpy " << tpy
		     << " tpz " << tpz
		     << " pos_eta " << pos_eta
		     << " pos_gpT " << pos_pT
		     << endl;
	      	      
	      // Make G4 invariant mass 
	      
	      TLorentzVector t2;
	      double E2 = sqrt( pow(tpx,2) + pow(tpy,2) + pow(tpz,2) + pow(decaymass,2));
	      t2.SetPxPyPzE(tpx,tpy,tpz,E2);	  
	      
	      TLorentzVector t = t1+t2;
	      
	      if(verbose)
		cout << "                       reco'd g4 mass = " << t.M() << endl << endl;	    
	      
	      if(t.M() > 7.0 && t.M() < 11.0)
		{
		  nrecog4mass++;
		  g4mass->Fill(t.M());	
		  hgpt->Fill(t.Pt());

		  // Capture the mass spectrum where both tracks are the primary Upsilon decay electrons
		  if( (gtrid1 == 1 || gtrid1 == 2) && (gtrid2 == 1 || gtrid2 == 2) ) 
		    g4mass_primary->Fill(t.M());	  

		}
	    }  // end of ipos loop
	} // end of ielec loop
      
      
      if(verbose)
	{
	  cout << " # of g4 electron tracks = " << ng4trevt_elec 
	       << " # of g4 positron tracks = " << ng4trevt_pos << endl;
	}
	  
	  
      //=============================
      // process reconstructed tracks
      // for this event
      //=============================
      
      int nrtrevt = 0;
      int nr_elec = -1;
      int nr_pos = -1;
      int rectrnum_elec[1000];
      int rectrnum_pos[1000];

      for(int ir=nr;ir<nr+ntracks;ir++)
        {
          int recoget = ntp_track->GetEntry(ir);
	  
	  hrquality->Fill(rquality);
	  hrdca2d->Fill(rdca2d);

	  // track quality cuts	
	  if(rquality > 3 || fabs(rdca2d) > 0.1)
	    continue;

	  // need to select electrons and positrons - for now we cheat
	  if(rgflavor != 11 && rgflavor != -11)
	    continue;

	  // make a list of electrons and positrons
	  if(rgflavor == 11)
	    {
	      nr_elec++;
	      rectrnum_elec[nr_elec] = ir;	      
	    }

	  if(rgflavor == -11)
	    {
	      nr_pos++;
	      rectrnum_pos[nr_pos] = ir;	      
	    }


	  double rpT = sqrt(rpx*rpx+rpy*rpy);
	  double reta = asinh(rpz/sqrt(rpx*rpx+rpy*rpy));
	  
	  // print out track details
	  if(verbose)
	    cout << "     revent " << revent << " ir " << ir
		 << " rgtrackid " << rgtrackid
		 << " rgflavor " << rgflavor
		 << " rvz " << rvz
		 << " reta " << reta
		 << " rpT " << rpT
		 << endl;
	}

      nr_elec++;
      nr_pos++;  

      for(int ielec = 0;ielec<nr_elec;ielec++)
	{

	  TLorentzVector t1;
	  
	  int recoget1 = ntp_track->GetEntry(rectrnum_elec[ielec]);

	  int trid1 = rgtrackid;	  

	  double E1 = sqrt( pow(rpx,2) + pow(rpy,2) + pow(rpz,2) + pow(decaymass,2));
	  t1.SetPxPyPzE(rpx,rpy,rpz,E1);	  
	  
	  for(int ipos = 0;ipos<nr_pos;ipos++)
	    {
	      int recoget2 = ntp_track->GetEntry(rectrnum_pos[ipos]);

	      int trid2 = rgtrackid;	  
	  
	      TLorentzVector t2;
	      double E2 = sqrt( pow(rpx,2) + pow(rpy,2) + pow(rpz,2) + pow(decaymass,2));
	  
	      t2.SetPxPyPzE(rpx,rpy,rpz,E2);	  
	  
	      TLorentzVector t = t1+t2;
	      
	      if(verbose)
		cout << " reco'd track mass = " << t.M() << endl;	    
	      
	      if(t.M() > 7.0 && t.M() < 11.0)
		{
		  nrecormass++;
		  recomass->Fill(t.M());	  
		  hrpt->Fill(t.Pt());

		  // Capture the mass spectrum where both tracks are the primary Upsilon decay electrons
		  if( (trid1 == 1 || trid1 == 2) && (trid2 == 1 || trid2 == 2) ) 
		    recomass_primary->Fill(t.M());	  
		}
	    }
	}
      
      // Increment the track entry numbers so we go to the 1st track in the next event
      nr += ntracks;
      ng += ngtracks;
    }

  cout << "nrecog4mass = " << nrecog4mass << endl;

  cout << "nrecormass = " << nrecormass << endl;

  //======================================================
  // End of loop over events and generation of mass histos
  //=====================================================

  TCanvas *cq = new TCanvas("cq","cq",5,5,600,600 );
  cq->Divide(1,2);
  cq->cd(1);
  hrquality->Draw();
  cq->cd(2);
  gPad->SetLogy(1);
  hrdca2d->Draw();

  // Mass histos
  
  TCanvas *cmass = new TCanvas("cmass","cmass",10,10,800,600);

  //cmass->cd(2);
  recomass_primary->SetLineColor(kRed);
  recomass_primary->DrawCopy();  

  recomass->SetLineColor(kBlack);
  recomass->DrawCopy("same");  


  TCanvas *cm_comp = new TCanvas("cm_comp","cm_comp",10,10,800,800);
  cm_comp->Divide(2,1);
  cm_comp->cd(1);
  recomass_primary->Draw();
  recomass->Draw("same");
  // we want from 7 to 11 GeV/c^2 - the whole range
  double yreco = recomass->Integral();
 double yreco_primary = recomass_primary->Integral();
 cout << "Reconstructed mass spectrum has " << yreco_primary << " entries from primary tracks and " << yreco << " entries total " << endl;

  cm_comp->cd(2);
  
  g4mass_primary->SetLineColor(kRed);
  g4mass_primary->Draw();
  g4mass->Draw("same");
  double yg4_primary = g4mass_primary->Integral();
  double yg4 = g4mass->Integral();
  cout << "G4 mass spectrum has " << yg4_primary << " entries from primary tracks and  " << yg4 << " entries total" << endl;

  cout << "Reconstruction efficiency is " << yreco_primary/yg4_primary << endl;


  // Output mass histos for individual Upsilon states


  char fname[500];

  // The Upsilon state is read from the input file
  int state = 3;   // Upsilon state = 1 or 2 or 3
  if(ups1s)
    state = 1;
  if(ups2s)
    state = 2;

  sprintf(fname,"ups%is_qual%.2f_dca2d%.2f.root",state,quality_cut,dca_cut);

  cout << "Create output file " << fname << endl;

  TFile *fout1 = new TFile(fname,"recreate");
  recomass->Write();
  recomass_primary->Write();
  g4mass->Write();
  g4mass_primary->Write();
  hrpt->Write();
  hrquality->Write();
  hrdca2d->Write();
  fout1->Close();

  cout << "Finished write to file " << fname << endl;

}
