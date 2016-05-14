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

#define TEST

void quarkonia_reconstruction()
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
  //TChain *ntp_vertex = new TChain("ntp_vertex","events");
  TChain *ntp_cluster = new TChain("ntp_cluster","clusters");

#include "ntuple_files/ntuple_files_ups3s_FPHX.C"

#endif
  
  // Ntuple access variables
  // This include file contains the definitions of the ntuple variables                                                                        
#include "ntuple_variables.C"
  
  // Define some histograms
  
  int nby = 50;
  int nbphi = 50;
  int nbeta = 50;
  int nbth = 50;
  int nbpt = 20;
  double ptmax = 10.0;

  TH1F* hrquality = new TH1F("hrquality", "Reconstructed track Quality", 1000, 0, 10);
  TH1F* hrdca2d = new TH1F("hrdca2d", "Reconstructed track dca2d", 1000, 0, 0.05);
  TH1F* hrcharge = new TH1F("hrcharge","reconstructed charge", 300, 0, 30);
  TH1F* hrpx = new TH1F("hrpx","reconstructed px; p_{x}", 300, 0, 30);
  TH1F* hrpy = new TH1F("hrpy","reconstructed px; p_{y}", 300, 0, 30);
  TH1F* hrpz = new TH1F("hrpz","reconstructed px; p_{z}", 300, 0, 30);
  
  TH1F* hrpt = new TH1F("hrpt"," pT", nbpt, 0.0, ptmax);
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

  TH1F* hgpt = new TH1F("hgpt","g4 pT", nbpt, 0.0, ptmax);
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
  TH1D *g4mass = new TH1D("g4mass","G4 input invariant mass",100,7.0,11.0);
  g4radmass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
  g4mass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");

  TH1D *recomass = new TH1D("recomass","Reconstructed invariant mass",200,7.0,11.0);
  recomass->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");


  //================================================================================
  // This section determines how many pairs to accumulate before quitting
  // The yields are from yield_revised_MIE.C
  // Note that all of these numbers must be inflated by the correct
  // factor so that the area of each peak is correct after subtraction
  // of the exponential background
  // The numbers to inflate by are the ratio of quarkonia_reconstruction.C yields divided by the CBfitter.C peak areas

  //=================================================================================
 
  double upsilon_inflation[3] = {1.015,1.06,1.14};

  int nups_requested;

  bool do_all = false;

  bool pp = true;
  bool cms_raa = false;
  bool strickland_raa = false;
  bool pAu = false;
  bool AuAu20pc = false;

  if(pp)
    {
      // new luminosity of 175 /pb, from yield_revised_MIE.C
      // with new acceptance of 38% with quality cut 3.0, |dca2d| cut 0.1 cm
      //from yield_revised_MIE.C, with pair eID eff of 0.9 for p+p and a trigger efficiency of 98%

      if(ups1s)
	nups_requested = 8769 * upsilon_inflation[0];
      else if(ups2s)
	nups_requested = 2205 * upsilon_inflation[1];
      else
	nups_requested = 1156 * upsilon_inflation[2];
    }
  else if (strickland_raa)
    {
      // central collisions only
      // These yields are from yield_revised_MIE.C, assuming 100B events MB for Au+Au
      // out of date: assumes 0.511 pair efficieincy and is for 0-20%
     
      if(ups1s)
	nups_requested = 6823;
      else if(ups2s)
	nups_requested = 579;
      else
	nups_requested = 64;
    }
  else if(pAu)
    {
      // from yield_revised_MIE.C
      // pAu 0-20%, 1200 /nb sampled, no suppression 
      // 0.38 pair efficiency
      // assumes eID pair eff of 0.8 for central, 0.9 for peripheral

      if(ups1s)
	nups_requested = 2945 * upsilon_inflation[0];
      else if(ups2s)
	nups_requested = 741 * upsilon_inflation[1];
      else
	nups_requested = 388 * upsilon_inflation[2];
    }
  else if(AuAu20pc)
    {
     // AuAu 0-20% no suppression
      // From yield_revised_MIE.C, assuming 100B events MB for Au+Au
      // 0.38 pair efficiency
      // assumes eID pair eff of 0.49 for central, 0.9 for peripheral

      // to add inefficiency
      double tempfactor = 1.0;

      if(ups1s)
	//nups_requested = 12313 * upsilon_inflation[0] * tempfactor;  // incorrect 51% efficiency number
	nups_requested = 9188 * upsilon_inflation[0] * tempfactor;
      else if(ups2s)
	//nups_requested = 3096 * upsilon_inflation[1] * tempfactor;
	nups_requested = 2311 * upsilon_inflation[1] * tempfactor;
      else
	//nups_requested = 1623 * upsilon_inflation[2] * tempfactor;
	nups_requested = 1212 * upsilon_inflation[2] * tempfactor;
    }
  else
    {
      // AuAu 0-10% no suppression
      // from yield_revised_MIE.C, assuming 100B MB events for Au+Au
      // 0.38 pair efficiency

      // to add inefficiency
      double tempfactor = 1.0;

      if(ups1s)
	nups_requested = 5626 * upsilon_inflation[0] * tempfactor;
      else if(ups2s)
	nups_requested = 1414 * upsilon_inflation[1] * tempfactor;
      else if(ups3s)
	nups_requested = 742 * upsilon_inflation[2] * tempfactor;
    }

  if(do_all)
    nups_requested = 200000;

  cout << "Upsilons requested = " << nups_requested << endl;


  //=======================
  // Loop over events
  //=======================

  int nr = 0;
  int ng = 0;
  int ng4evtgood = 0;
  int ng4trgood = 0;
  int nrecog4mass = 0;
  int nrevtgood = 0;
  int nrtrgood = 0;
  int nrecormass = 0;

  for(int iev=0;iev<ntp_vertex->GetEntries();iev++)
    {
      // drop out when the requested number of reco'd Upsilons has been reached
      if(nrecormass > nups_requested)
	break;

      int recoget = ntp_vertex->GetEntry(iev);

      int Nskip = 0;
      // Skip the first N vertexs
      if(iev < Nskip)
	{
	  nr += ntracks;
	  ng += ngtracks;
	  continue;
	}


      if(verbose)
	cout << "iev " << iev
	     << " event " << event
	     << " ntracks " << ntracks
	     << " ngtracks " << ngtracks
	     << " gvz " << egvz
	     << " vz " << evz
	     << endl;
      
      //============================
      // process G4 tracks
      // for this event
      //============================

      int ng4trevt = 0;
      int nent1 = -1;
      int nent2 = -1;
      for(int ig=ng;ig<ng+ngtracks;ig++)
        {
          int recoget1 = ntp_gtrack->GetEntry(ig);
      
	  double gpT = sqrt(tpx*tpx+tpy*tpy);
	  double eta = asinh(tpz/sqrt(tpx*tpx+tpy*tpy));
	  
	  if(tgtrackid ==1 || tgtrackid == 2)
	    {
	      ng4trgood++;
	      ng4trevt++;

	      if(nent1 < 0)
		nent1 = ig;
	      else
		nent2 = ig;

	      // print out track details
	      if(verbose)
		cout << "     tevent " << tevent << " ig " << ig
		     << " tgtrackid " << tgtrackid
		     << " tflavor " << tflavor
		     << " gvx " << tvx
		     << " gvy " << tvy
		     << " gvz " << tvz
		     << " eta " << eta
		     << " gpT " << gpT
		     << endl;
	    }
	}

      if(verbose)
	cout << " # of g4 tracks = " << ng4trevt << endl;

      if(ng4trevt == 2)
	ng4evtgood++;

      // Make G4 invariant mass if possible
      if(nent1 >= 0 && nent2 >+ 0)
	{
	  TLorentzVector t1, t2;

	  int recoget1 = ntp_gtrack->GetEntry(nent1);
	  
	  double E1 = sqrt( pow(tpx,2) + pow(tpy,2) + pow(tpz,2) + pow(decaymass,2));
	    
	  t1.SetPxPyPzE(tpx,tpy,tpz,E1);	  

	  int recoget2 = ntp_gtrack->GetEntry(nent2);
	  
	  double E2 = sqrt( pow(tpx,2) + pow(tpy,2) + pow(tpz,2) + pow(decaymass,2));
	    
	  t2.SetPxPyPzE(tpx,tpy,tpz,E2);	  

	  TLorentzVector t = t1+t2;

	  if(verbose)
	    cout << " reco'd g4 mass = " << t.M() << endl;	    

	  if(t.M() > 7.0 && t.M() < 11.0)
	    {
	      nrecog4mass++;
	      g4mass->Fill(t.M());	
	      hgpt->Fill(t.Pt());
	    }
	}

      //=============================
      // process reconstructed tracks
      // for this event
      //=============================

      int nrtrevt = 0;
      nent1 = -1;
      nent2 = -1;
      for(int ir=nr;ir<nr+ntracks;ir++)
        {
          int recoget = ntp_track->GetEntry(ir);

	  hrquality->Fill(rquality);
	  hrdca2d->Fill(rdca2d);

	  if(rquality > 3 || fabs(rdca2d) > 0.1)
	    continue;
      
	  double rpT = sqrt(rpx*rpx+rpy*rpy);
	  double reta = asinh(rpz/sqrt(rpx*rpx+rpy*rpy));

	  if(rgtrackid ==1 || rgtrackid == 2)
	    {
	      nrtrgood++;
	      nrtrevt++;

	      if(nent1 < 0)
		nent1 = ir;
	      else
		nent2 = ir;

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
	}

      if(verbose)
	cout << " # of reco'd tracks = " << nrtrevt << endl;

      if(nrtrevt == 2)
	nrevtgood++;

      // Make reco'd invariant mass if possible
      if(nent1 >= 0 && nent2 >= 0)
	{
	  TLorentzVector t1, t2;

	  int recoget1 = ntp_track->GetEntry(nent1);
	  
	  double E1 = sqrt( pow(rpx,2) + pow(rpy,2) + pow(rpz,2) + pow(decaymass,2));
	    
	  t1.SetPxPyPzE(rpx,rpy,rpz,E1);	  

	  int recoget2 = ntp_track->GetEntry(nent2);
	  
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
	    }
	}



      nr += ntracks;
      ng += ngtracks;
      
    }

  cout << "ng4evtgood = " << ng4evtgood << endl;
  cout << "ng4trgood = " << ng4trgood << endl;
  cout << "nrecog4mass = " << nrecog4mass << endl;

  cout << "nrevtgood = " << nrevtgood << endl;
  cout << "nrtrgood = " << nrtrgood << endl;
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
  recomass->SetLineColor(kBlack);
  recomass->DrawCopy();  

  TCanvas *cm_comp = new TCanvas("cm_comp","cm_comp",10,10,800,800);
  cm_comp->Divide(2,1);
  cm_comp->cd(1);
  recomass->Draw();
  // we want from 7 to 11 GeV/c^2 - the whole range
  double yreco = recomass->Integral();
  cout << "Reconstructed mass spectrum has " << yreco << " entries " << endl;

  cm_comp->cd(2);
  g4mass->Draw();
  double yg4 = g4mass->Integral();
  cout << "G4 mass spectrum has " << g4mass->Integral() << " entries " << endl;

  cout << "Reconstruction efficiency is " << yreco/yg4 << endl;


  // Output mass histos for individual Upsilon states


  char fname[500];

  // The Upsilon state is read from the input file
  int state = 3;   // Upsilon state = 1 or 2 or 3
  if(ups1s)
    state = 1;
  if(ups2s)
    state = 2;

  if(pp)
    {
      sprintf(fname,"ups%is_qual%.2f_dca2d%.2f_pp.root",state,quality_cut,dca_cut);
    }
  else if(pAu)
    {
      sprintf(fname,"ups%is_qual%.2f_dca2d%.2f_pAu.root",state,quality_cut,dca_cut);
    }
  else if(AuAu20pc)
    {
      sprintf(fname,"ups%is_qual%.2f_dca2d%.2f_AuAu_20pc.root",state,quality_cut,dca_cut);
    }
  else if(do_all)
    {
      sprintf(fname,"ups%is_qual%.2f_dca2d%.2f_allfiles.root",state,quality_cut,dca_cut);
    }
  else
    {
      sprintf(fname,"ups%is_qual%.2f_dca2d%.2f_AuAu_10pc.root",state,quality_cut,dca_cut);
    }

  cout << "Create output file " << fname << endl;

  TFile *fout1 = new TFile(fname,"recreate");
  recomass->Write();
  hrpt->Write();
  fout1->Close();

  cout << "Finished write to file " << fname << endl;

  // Output pT histo for Upsilons

  /*
  char fname1[500];

  sprintf(fname1,"accepted_upsilon_pT_%i_20pc_central.root",state);
  cout << "Create output file " << fname1 << endl;
  TFile *fout = new TFile(fname1,"recreate");
  cout << "Write output histos for pT distributions to " << fname1 << endl;
  hgpt->Write();
  hrpt->Write();

  fout->Close();
  */

}
