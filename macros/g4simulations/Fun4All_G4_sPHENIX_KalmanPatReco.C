int Fun4All_G4_sPHENIX_KalmanPatReco (
				const int nproc = 5000,
				const int nEvents = 1,
				const int which_tracking = 15,
				//const char * inputFile = NULL,
				const char * outputFile = "SvtxTracks.root",
				const char * embed_input_file = "Hijing_G4Hits.root",
				const bool do_embedding = false
				)
{
  char inputFile[500];
  
  //======================
  // What to run
  //======================

	bool output_tracks = false;

	// pions for momentum, dca performance
	bool pion_momentum = false;
	bool embed_pions = false;  // throw single pions in a Hijing event if true
	// Upsilons
	bool upsilons = true;           // throw single Upsilons if true
	int istate = 1;  // Upsilon state = 1,2,3
	bool embed_upsilons = false;           // if true, throw single Upsilons inside a Hijing event
	// Hijing events only
	bool hijing_events = false;  // if true, throw hijing events only, or with embedded pions or upsilons


	//===============
	// Input options
	//===============

	// Either:
	// read previously generated g4-hits files, in this case it opens a DST and skips
	// the simulations step completely. The G4Setup macro is only loaded to get information
	// about the number of layers used for the cell reco code
	//
	// In case reading production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
	// E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
	const bool readhits = false;
	// Or:

	// read files in HepMC format (typically output from event generators like hijing or pythia)
	bool readhepmc = false; // read HepMC files
	if(hijing_events || embed_upsilons || embed_pions)
	  readhepmc = true;

	// Or:
	// Use pythia
	const bool runpythia8 = false;
	const bool runpythia6 = false;
	// else
	// Use particle generator (default simple generator)
	// or gun/ very simple generator
	const bool usegun = false;
	// And
	// Further choose to embed newly simulated events to a previous simulation. Not compatible with `readhits = true`
	// In case embedding into a production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
	// E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
	//const bool do_embedding = true;

	if(hijing_events || embed_upsilons || embed_pions)
	  {
	    // get the Hijing input file name
	    sprintf(inputFile,"/phenix/hhj/frawley/tracking/stage1_jobs/in/hijing_%.5i.txt.bz2",nproc);
	    
	    cout << "Reading Hijing events from file: " << endl << inputFile << endl; 
	  }
	
	bool do_bbc = true;

	bool do_pipe = true;

	bool do_svtx = true;
	bool do_svtx_cell = do_svtx && true;
	bool do_svtx_track = do_svtx_cell && true;
	bool do_svtx_eval = do_svtx_track && true;
	bool do_svtx_eval = do_svtx_track && true;

	cout << "Tracking option: which_tracking = " << which_tracking  << endl;

	bool do_preshower = false;

	bool do_cemc = false;
	bool do_cemc_cell = do_cemc && true;
	bool do_cemc_twr = do_cemc_cell && true;
	bool do_cemc_cluster = do_cemc_twr && true;
	bool do_cemc_eval = do_cemc_cluster && true;

	bool do_hcalin = false;
	bool do_hcalin_cell = do_hcalin && true;
	bool do_hcalin_twr = do_hcalin_cell && true;
	bool do_hcalin_cluster = do_hcalin_twr && true;
	bool do_hcalin_eval = do_hcalin_cluster && true;

	bool do_magnet = false;

	bool do_hcalout = false;
	bool do_hcalout_cell = do_hcalout && true;
	bool do_hcalout_twr = do_hcalout_cell && true;
	bool do_hcalout_cluster = do_hcalout_twr && true;
	bool do_hcalout_eval = do_hcalout_cluster && true;

	bool do_global = true;
	bool do_global_fastsim = false;

	bool do_jet_reco = true;
	bool do_jet_eval = false;

	bool do_dst_compress = false;

	//Option to convert DST to human command readable TTree for quick poke around the outputs
	bool do_DSTReader = false;
	//---------------
	// Load libraries
	//---------------

	gSystem->Load("libfun4all.so");
	gSystem->Load("libg4detectors.so");
	gSystem->Load("libphhepmc.so");
	gSystem->Load("libg4testbench.so");
	gSystem->Load("libg4hough.so");
	gSystem->Load("libcemc.so");
	gSystem->Load("libg4eval.so");

	// establish the geometry and reconstruction setup
	gROOT->LoadMacro("G4Setup_sPHENIX_KalmanPatReco.C");
	G4Init(do_svtx,do_preshower,do_cemc,do_hcalin,do_magnet,do_hcalout,do_pipe,which_tracking);

	int absorberactive = 1; // set to 1 to make all absorbers active volumes
	//  const string magfield = "1.5"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
	const string magfield = "/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
	const float magfield_rescale = 1.4/1.5; // scale the map to a 1.4 T field

	//---------------
	// Fun4All server
	//---------------

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(0);
	// just if we set some flags somewhere in this macro
	recoConsts *rc = recoConsts::instance();
	// By default every random number generator uses
	// PHRandomSeed() which reads /dev/urandom to get its seed
	// if the RANDOMSEED flag is set its value is taken as seed
	// You ca neither set this to a random value using PHRandomSeed()
	// which will make all seeds identical (not sure what the point of
	// this would be:
	//  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
	// or set it to a fixed value so you can debug your code
	//  rc->set_IntFlag("RANDOMSEED", 12345);

	//-----------------
	// Event generation
	//-----------------
	
	if (readhits)
	  {
	    // Get the hits from a file
	    // The input manager is declared later
	    
	    if (do_embedding)
	      {
		cout <<"Do not support read hits and embed background at the same time."<<endl;
		exit(1);
	      }
	    
	  }
	else if (readhepmc)
	  {
	    // this module is needed to read the HepMC records into our G4 sims
	    // but only if you read HepMC input files
	    HepMCNodeReader *hr = new HepMCNodeReader();
	    se->registerSubsystem(hr);
	  }
	else if (runpythia8)
	  {
	    gSystem->Load("libPHPythia8.so");
	    
	    PHPy8JetTrigger *theTrigger = new PHPy8JetTrigger();
	    //theTrigger->Verbosity(10);
	    theTrigger->SetEtaHighLow(-0.6, 0.6); 
	    theTrigger->SetJetR(.4);
	    theTrigger->SetMinJetPt(20);
	    
	    PHPythia8* pythia8 = new PHPythia8();
	    // see coresoftware/generators/PHPythia8 for example config
	    pythia8->set_config_file("phpythia8.cfg"); 
	    pythia8->register_trigger(theTrigger);
	    se->registerSubsystem(pythia8);
	    
	    HepMCNodeReader *hr = new HepMCNodeReader();
	    hr->Embed(10);
	    se->registerSubsystem(hr);
	  }
	else if (runpythia6)
	  {
	    gSystem->Load("libPHPythia6.so");
	    
	    PHPythia6 *pythia6 = new PHPythia6();
	    pythia6->set_config_file("phpythia6.cfg");
	    se->registerSubsystem(pythia6);
	    
	    HepMCNodeReader *hr = new HepMCNodeReader();
	    se->registerSubsystem(hr);
	  }

	// run pions for momentum, dca performance, alone or embedded in Hijing
	if (pion_momentum || embed_pions)
	  {
	    cout << "Throw 100 pions" << endl;
	    // throw embedded pions to 50 GeV/c in 0.5 GeV/c intervals      
	    for(int i=0; i<100; i++)
	      {
		double pt = (double) i * 0.5 + 0.25;
		
		// toss low multiplicity dummy events
		PHG4SimpleEventGenerator *pgen = new PHG4SimpleEventGenerator();
		pgen->add_particles("pi+",1); // mu-,e-,anti_proton,pi-
		//pgen->add_particles("pi-",1); // mu-,e-,anti_proton,pi-
		
		if (readhepmc) {
		  pgen->set_reuse_existing_vertex(true);
		  pgen->set_existing_vertex_offset_vector(0.0,0.0,0.0);
		} else {
		  pgen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
							 PHG4SimpleEventGenerator::Uniform,
							 PHG4SimpleEventGenerator::Uniform);
		  pgen->set_vertex_distribution_mean(0.0,0.0,0.0);
		  //pgen->set_vertex_distribution_width(0.0,0.0,5.0);
		  pgen->set_vertex_distribution_width(0.0,0.0,0.0);
		}
		pgen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
		pgen->set_vertex_size_parameters(0.0,0.0);
		pgen->set_eta_range(-1.0, 1.0);
		pgen->set_phi_range(-1.0*TMath::Pi(), 1.0*TMath::Pi());
		pgen->set_pt_range(pt, pt);
		
		pgen->Embed(1);
		pgen->Verbosity(0);
		se->registerSubsystem(pgen);	  
	      }    
	    
	  }

	// run upsilons for momentum, dca performance, alone or embedded in Hijing
	if(upsilons || embed_upsilons)
	  {
	    PHG4ParticleGeneratorVectorMeson *vgen = new PHG4ParticleGeneratorVectorMeson();
	    vgen->set_decay_types("e+","e-");    // dielectron decay
	    //vgen->set_vtx_zrange(-10.0, +10.0);
	    vgen->set_vtx_zrange(0.0, 0.0);
	    // Note: this rapidity range completely fills the acceptance of eta = +/- 1 unit
	    vgen->set_rapidity_range(-1.0, +1.0);
	    vgen->set_pt_range(0.0, 10.0);
	    
	    if(istate == 1)
	      {
		// Upsilon(1S)
		vgen->set_mass(9.46);
		vgen->set_width(54.02e-6);
	      }
	    else if (istate == 2)
	      {
		// Upsilon(2S)
		vgen->set_mass(10.0233);
		vgen->set_width(31.98e-6);
	      }
	    else
	      {
		// Upsilon(3S)
		vgen->set_mass(10.3552);
		vgen->set_width(20.32e-6);
	      }
	    
	    vgen->Verbosity(0);
	    se->registerSubsystem(vgen);

	    cout << "Upsilon generator for istate = " << istate << " created and registered "  << endl;	  
	    
	  }

	/*
	  {
	    // toss low multiplicity dummy events
	    PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
	    // mu+,e+,proton,pi+,Upsilon
	    gen->add_particles("pi+",10);
	    //gen->add_particles("pi-",1);
	    if (readhepmc || do_embedding)
	      {
		gen->set_reuse_existing_vertex(true);
		gen->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	      }
	    else
	      {
		gen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
						      PHG4SimpleEventGenerator::Uniform,
						      PHG4SimpleEventGenerator::Uniform);
		gen->set_vertex_distribution_mean(0.0, 0.0, 0.0);
		gen->set_vertex_distribution_width(0.0, 0.0, 0.0);
	      }
	    gen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
	    gen->set_vertex_size_parameters(0.0, 0.0);
	    gen->set_eta_range(-0.5, 0.5);
	    gen->set_phi_range(-1.0 * TMath::Pi(), 1.0 * TMath::Pi());
	    //		gen->set_eta_range(0, 0);
	    //		gen->set_phi_range(0, 0);
	    gen->set_pt_range(0, 40);
	    gen->Embed(10);
	    gen->Verbosity(0);
	    if (! usegun)
	      {
		se->registerSubsystem(gen);
	      }
	    else
	      {
		PHG4ParticleGun *gun = new PHG4ParticleGun();
		//  gun->set_name("anti_proton");
		gun->set_name("geantino");
		gun->set_vtx(0, 0, 0);
		gun->set_mom(10, 0, 0.01);
		// gun->AddParticle("geantino",1.7776,-0.4335,0.);
		// gun->AddParticle("geantino",1.7709,-0.4598,0.);
		// gun->AddParticle("geantino",2.5621,0.60964,0.);
		// gun->AddParticle("geantino",1.8121,0.253,0.);
		//	  se->registerSubsystem(gun);
		PHG4ParticleGenerator *pgen = new PHG4ParticleGenerator();
		pgen->set_name("geantino");
		pgen->set_z_range(0,0);
		pgen->set_eta_range(0.01,0.01);
		pgen->set_mom_range(10,10);
		pgen->set_phi_range(5.3./180.*TMath::Pi(),5.7./180.*TMath::Pi());
		se->registerSubsystem(pgen);
		pgen = new PHG4ParticleGenerator();
		pgen->set_name("geantino");
		pgen->set_z_range(0,0);
		pgen->set_eta_range(0.01,0.01);
		pgen->set_mom_range(10,10);
		pgen->set_phi_range(-0.2./180.*TMath::Pi(),0.2./180.*TMath::Pi());
		se->registerSubsystem(pgen);
	      }
	  }
	*/

	if (!readhits)
	{
		//---------------------
		// Detector description
		//---------------------

		G4Setup(absorberactive, magfield, TPythia6Decayer::kAll,
				do_svtx, do_preshower, do_cemc, do_hcalin, do_magnet, do_hcalout, do_pipe, magfield_rescale);
	}

	//---------
	// BBC Reco
	//---------

	if (do_bbc) 
	{
		gROOT->LoadMacro("G4_Bbc.C");
		BbcInit();
		Bbc_Reco();
	}
	//------------------
	// Detector Division
	//------------------

	if (do_svtx_cell) Svtx_Cells();

	if (do_cemc_cell) CEMC_Cells();

	if (do_hcalin_cell) HCALInner_Cells();

	if (do_hcalout_cell) HCALOuter_Cells();

	//-----------------------------
	// CEMC towering and clustering
	//-----------------------------

	if (do_cemc_twr) CEMC_Towers();
	if (do_cemc_cluster) CEMC_Clusters();

	//-----------------------------
	// HCAL towering and clustering
	//-----------------------------

	if (do_hcalin_twr) HCALInner_Towers();
	if (do_hcalin_cluster) HCALInner_Clusters();

	if (do_hcalout_twr) HCALOuter_Towers();
	if (do_hcalout_cluster) HCALOuter_Clusters();

	if (do_dst_compress) ShowerCompress();

	//--------------
	// SVTX tracking
	//--------------

	if (do_svtx_track) Svtx_Reco();

	//-----------------
	// Global Vertexing
	//-----------------

	if (do_global) 
	{
		gROOT->LoadMacro("G4_Global.C");
		Global_Reco();
	}

	else if (do_global_fastsim) 
	{
		gROOT->LoadMacro("G4_Global.C");
		Global_FastSim();
	}  

	//---------
	// Jet reco
	//---------

	if (do_jet_reco) 
	{
		gROOT->LoadMacro("G4_Jets.C");
		Jet_Reco();
	}

	// HF jet trigger moudle
	//assert (gSystem->Load("libHFJetTruthGeneration") == 0); 
	//{
	//	if (do_jet_reco)
	//	{   
	//		HFJetTruthTrigger * jt = new HFJetTruthTrigger(
	//				"HFJetTruthTrigger.root", 5 , "AntiKt_Truth_r04");
	//		//jt->Verbosity(HFJetTruthTrigger::VERBOSITY_MORE);
	//		jt->set_pt_min(20);
	//		se->registerSubsystem(jt);
	//	}   
	//}

	//----------------------
	// Simulation evaluation
	//----------------------

	char evalFile[500];
	sprintf(evalFile,"eval_output/g4svtx_eval_%i.root",nproc);
	if (do_svtx_eval) Svtx_Eval(evalFile);

	if (do_cemc_eval) CEMC_Eval("g4cemc_eval.root");

	if (do_hcalin_eval) HCALInner_Eval("g4hcalin_eval.root");

	if (do_hcalout_eval) HCALOuter_Eval("g4hcalout_eval.root");

	if (do_jet_eval) Jet_Eval("g4jet_eval.root");

	//-------------- 
	// IO management
	//--------------

	if (readhits)
	{
		// Hits file
		Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
		hitsin->fileopen(inputFile);
		se->registerInputManager(hitsin);
	}
	if (do_embedding)
	{
		if (embed_input_file == NULL)
		{
			cout << "Missing embed_input_file! Exit";
			exit(3);
		}

		Fun4AllDstInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTinEmbed");
		in1->AddFile(embed_input_file); // if one use a single input file
		//in1->AddListFile(embed_input_file); // RecommendedL: if one use a text list of many input files
		se->registerInputManager(in1);
	}
	if (readhepmc)
	{
		Fun4AllInputManager *in = new Fun4AllHepMCInputManager( "DSTIN");
		se->registerInputManager( in );
		se->fileopen( in->Name().c_str(), inputFile );
	}
	else
	{
		// for single particle generators we just need something which drives
		// the event loop, the Dummy Input Mgr does just that
		Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
		se->registerInputManager( in );
	}

	char tracksFile[500];
	sprintf(tracksFile,"eval_output/SvtxTracks_%i.root",nproc);
	if (do_DSTReader)
	{
		//Convert DST to human command readable TTree for quick poke around the outputs
		gROOT->LoadMacro("G4_DSTReader.C");

		G4DSTreader( tracksFile, //
				/*int*/ absorberactive ,
				/*bool*/ do_svtx ,
				/*bool*/ do_preshower ,
				/*bool*/ do_cemc ,
				/*bool*/ do_hcalin ,
				/*bool*/ do_magnet ,
				/*bool*/ do_hcalout ,
				/*bool*/ do_cemc_twr ,
				/*bool*/ do_hcalin_twr ,
				/*bool*/ do_magnet  ,
				/*bool*/ do_hcalout_twr
				);
	}

	if(output_tracks)
	  {
	    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", tracksFile);
	    if (do_dst_compress) DstCompress(out);
	    se->registerOutputManager(out);
	  }

	//-----------------
	// Event processing
	//-----------------
	if (nEvents < 0)
	{
		return;
	}
	// if we run the particle generator and use 0 it'll run forever
	if (nEvents == 0 && !readhits && !readhepmc)
	{
		cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
		cout << "it will run forever, so I just return without running anything" << endl;
		return;
	}

	se->run(nEvents);

	//-----
	// Exit
	//-----

	se->End();

	//getchar();

	std::cout << "All done" << std::endl;
	delete se;
	gSystem->Exit(0);
}
