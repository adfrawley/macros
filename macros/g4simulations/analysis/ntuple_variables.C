

  //===================================
  // access the event ntuple variables
  //===================================

  int nevents = ntp_vertex->GetEntries();
  cout << "Number of events " << nevents << endl;

//int num_g4hits = ntp_g4hit->GetEntries();
//cout << "Number of g4hits " << num_g4hits << endl;

int num_clusters = ntp_cluster->GetEntries();
cout << "Number of clusters " << num_clusters << endl;

  Float_t event;
  Float_t ntracks;
  Float_t ngtracks;
  Float_t ng4hits;
//Float_t nclusters;
  Float_t evx;
  Float_t evy;
  Float_t evz;
  Float_t egvx;
  Float_t egvy;
  Float_t egvz;
  Float_t hit_occupancy_layer0;
  Float_t hit_occupancy_layer1;
  Float_t hit_occupancy_layer2;
  Float_t hit_occupancy_layer3;
  Float_t hit_occupancy_layer4;
  Float_t hit_occupancy_layer5;
//Float_t hit_occupancy_layer6;

  TBranch* b_event;
  TBranch* b_ntracks;
  TBranch* b_ngtracks;
  TBranch* b_ng4hits;
//  TBranch* b_nclusters;
  TBranch* b_evx;
  TBranch* b_evy;
  TBranch* b_evz;
  TBranch* b_egvx;
  TBranch* b_egvy;
  TBranch* b_egvz;
//TBranch* b_hit_occupancy_layer0;
//TBranch* b_hit_occupancy_layer1;
//TBranch* b_hit_occupancy_layer2;
//TBranch* b_hit_occupancy_layer3;
//TBranch* b_hit_occupancy_layer4;
//TBranch* b_hit_occupancy_layer5;
//TBranch* b_hit_occupancy_layer6;

  ntp_vertex->SetBranchAddress("event", &event);
//ntp_vertex->SetBranchAddress("gnhits", &ng4hits);
  ntp_vertex->SetBranchAddress("ntracks", &ntracks);
  ntp_vertex->SetBranchAddress("gntracks", &ngtracks);
//ntp_vertex->SetBranchAddress("nclusters", &nclusters);
  ntp_vertex->SetBranchAddress("vx", &evx);
  ntp_vertex->SetBranchAddress("vy", &evy);
  ntp_vertex->SetBranchAddress("vz", &evz);
  ntp_vertex->SetBranchAddress("gvx", &egvx);
  ntp_vertex->SetBranchAddress("gvy", &egvy);
  ntp_vertex->SetBranchAddress("gvz", &egvz);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer0", &hit_occupancy_layer0);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer1", &hit_occupancy_layer1);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer2", &hit_occupancy_layer2);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer3", &hit_occupancy_layer3);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer4", &hit_occupancy_layer4);
//ntp_vertex->SetBranchAddress("hit_occupancy_layer5", &hit_occupancy_layer5);
//ntp_event->SetBranchAddress("hit_occupancy_layer6", &hit_occupancy_layer6);

  b_event = ntp_vertex->GetBranch("event");
  b_ntracks = ntp_vertex->GetBranch("ntracks");
  b_ngtracks = ntp_vertex->GetBranch("gntracks");
//b_ng4hits = ntp_vertex->GetBranch("gnhits");
//b_nclusters = ntp_vertex->GetBranch("nclusters");
  b_evx = ntp_vertex->GetBranch("vx");
  b_evy = ntp_vertex->GetBranch("vy");
  b_evz = ntp_vertex->GetBranch("vz");
  b_egvx = ntp_vertex->GetBranch("gvx");
  b_egvy = ntp_vertex->GetBranch("gvy");
  b_egvz = ntp_vertex->GetBranch("gvz");
//b_hit_occupancy_layer0 = ntp_vertex->GetBranch("hit_occupancy_layer0");
//b_hit_occupancy_layer1 = ntp_vertex->GetBranch("hit_occupancy_layer1");
//b_hit_occupancy_layer2 = ntp_vertex->GetBranch("hit_occupancy_layer2");
//b_hit_occupancy_layer3 = ntp_vertex->GetBranch("hit_occupancy_layer3");
//b_hit_occupancy_layer4 = ntp_vertex->GetBranch("hit_occupancy_layer4");
//b_hit_occupancy_layer5 = ntp_vertex->GetBranch("hit_occupancy_layer5");
//b_hit_occupancy_layer6 = ntp_vertex->GetBranch("hit_occupancy_layer6");

  //=====================================
  // Access ntp_track (reco'd) variables
  //=====================================

cout << "set up ntp_track access" << endl;

  Float_t rpx;
  Float_t rpy;
  Float_t rpz;
  Float_t rgfx;
  Float_t rgfy;
  Float_t rgfz;
  Float_t rgpx;
  Float_t rgpy;
  Float_t rgpz;
  Float_t rquality;
  Float_t rchisq;
//Float_t rchisqv;
  Float_t rcharge;
  Float_t revent;
  Float_t rtrackid;
  Float_t rgtrackid;
  Float_t rgflavor;
  Float_t rprimary;
  Float_t rpurity;
  Float_t rvz;
  Float_t rpcax;
  Float_t rpcay;
  Float_t rpcaz;
  Float_t rdca2d;
  Float_t rdca2dsigma;
  Float_t rnhits;
  Float_t rgnhits;
  Float_t rgembed;
  
  TBranch* b_px;
  TBranch* b_py;
  TBranch* b_pz;
  TBranch* b_gpx;
  TBranch* b_gpy;
  TBranch* b_gpz;
  TBranch* b_gfx;
  TBranch* b_gfy;
  TBranch* b_gfz;
  TBranch* b_quality;
  TBranch* b_chisq;
//TBranch* b_chisqv;
  TBranch* b_charge;
  TBranch* b_revent;
  TBranch* b_trackid;
  TBranch* b_gtrackid;
  TBranch* b_gflavor;
  TBranch* b_primary;
  TBranch* b_purity;
  TBranch* b_gvz;
  TBranch* b_pcax;
  TBranch* b_pcay;
  TBranch* b_pcaz;
  TBranch* b_dca2d;
  TBranch* b_dca2dsigma;
  TBranch* b_rnhits;
  TBranch* b_rgnhits;
  TBranch* b_rgembed;
  
  //set branches
  ntp_track->SetBranchAddress("px", &rpx);
  ntp_track->SetBranchAddress("py", &rpy);
  ntp_track->SetBranchAddress("pz", &rpz);
  ntp_track->SetBranchAddress("gpx", &rgpx);
  ntp_track->SetBranchAddress("gpy", &rgpy);
  ntp_track->SetBranchAddress("gpz", &rgpz);
  ntp_track->SetBranchAddress("gfx", &rgfx);
  ntp_track->SetBranchAddress("gfy", &rgfy);
  ntp_track->SetBranchAddress("gfz", &rgfz);
  ntp_track->SetBranchAddress("charge", &rcharge);
  ntp_track->SetBranchAddress("quality", &rquality);
  ntp_track->SetBranchAddress("chisq", &rchisq);
//ntp_track->SetBranchAddress("chisqv", &rchisqv);
  ntp_track->SetBranchAddress("event", &revent);
  ntp_track->SetBranchAddress("trackID", &rtrackid);
  ntp_track->SetBranchAddress("gtrackID", &rgtrackid);
  ntp_track->SetBranchAddress("gflavor", &rgflavor);
  ntp_track->SetBranchAddress("gprimary", &rprimary);
  ntp_track->SetBranchAddress("nfromtruth", &rpurity);
  ntp_track->SetBranchAddress("gvz", &rvz);
  ntp_track->SetBranchAddress("pcax", &rpcax);
  ntp_track->SetBranchAddress("pcay", &rpcay);
  ntp_track->SetBranchAddress("pcaz", &rpcaz);
  ntp_track->SetBranchAddress("dca2d", &rdca2d);
  ntp_track->SetBranchAddress("dca2dsigma", &rdca2dsigma);
  ntp_track->SetBranchAddress("nhits", &rnhits);
  ntp_track->SetBranchAddress("gnhits", &rgnhits);
  ntp_track->SetBranchAddress("gembed", &rgembed);
  
  //get Branches
  b_px = ntp_track->GetBranch("px");
  b_py = ntp_track->GetBranch("py");
  b_pz = ntp_track->GetBranch("pz");
  b_gpx = ntp_track->GetBranch("gpx");
  b_gpy = ntp_track->GetBranch("gpy");
  b_gpz = ntp_track->GetBranch("gpz");
  b_gfx = ntp_track->GetBranch("gfx");
  b_gfy = ntp_track->GetBranch("gfy");
  b_gfz = ntp_track->GetBranch("gfz");
  b_charge = ntp_track->GetBranch("charge");
  b_quality = ntp_track->GetBranch("quality");
//b_chisqv = ntp_track->GetBranch("chisqv");
  b_revent = ntp_track->GetBranch("event");
  b_trackid = ntp_track->GetBranch("trackID");
  b_gtrackid = ntp_track->GetBranch("gtrackID");
  b_gflavor = ntp_track->GetBranch("gflavor");
  b_primary = ntp_track->GetBranch("gprimary");
  b_purity = ntp_track->GetBranch("nfromtruth");
  b_gvz = ntp_track->GetBranch("gvz");
  b_pcax = ntp_track->GetBranch("pcax");
  b_pcay = ntp_track->GetBranch("pcay");
  b_pcaz = ntp_track->GetBranch("pcaz");
  b_dca2d = ntp_track->GetBranch("dca2d");
  b_dca2dsigma = ntp_track->GetBranch("dca2dsigma");
  b_rnhits = ntp_track->GetBranch("nhits");
  b_rnhits = ntp_track->GetBranch("gnhits");
  b_rgembed = ntp_track->GetBranch("gembed");
  
  //================================
  // Access the g4track ntuple variables
  //================================

cout << "set up ntp_gtrack access" << endl;

  Float_t tpx;
  Float_t tpy;
  Float_t tpz;
  Float_t tfpx;
  Float_t tfpy;
  Float_t tfpz;
  Float_t tevent;
  Float_t tgtrackid;
  Float_t tflavor;
  Float_t tvx;
  Float_t tvy;
  Float_t tvz;
  Float_t tnhits;
  Float_t tchisq;
  Float_t tprimary;
  Float_t tbestpurity;
  Float_t tembed;

  TBranch* b_tpx;
  TBranch* b_tpy;
  TBranch* b_tpz;
  TBranch* b_tfpx;
  TBranch* b_tfpy;
  TBranch* b_tfpz;
  TBranch* b_tevent;
  TBranch* b_tgtrackid;
  TBranch* b_tflavor;
  TBranch* b_tvx;
  TBranch* b_tvy;
  TBranch* b_tvz;
  TBranch* b_tnhits;
  TBranch* b_tchisq;
  TBranch* b_tprimary;
  TBranch* b_tbestpurity;
  TBranch* b_tembed;

  //set branches
  ntp_gtrack->SetBranchAddress("gpx", &tpx);
  ntp_gtrack->SetBranchAddress("gpy", &tpy);
  ntp_gtrack->SetBranchAddress("gpz", &tpz);
  ntp_gtrack->SetBranchAddress("gfpx", &tfpx);
  ntp_gtrack->SetBranchAddress("gfpy", &tfpy);
  ntp_gtrack->SetBranchAddress("gfpz", &tfpz);
  ntp_gtrack->SetBranchAddress("event", &tevent);
  ntp_gtrack->SetBranchAddress("gtrackID", &tgtrackid);
  ntp_gtrack->SetBranchAddress("gflavor", &tflavor);
  ntp_gtrack->SetBranchAddress("gvx", &tvx);
  ntp_gtrack->SetBranchAddress("gvy", &tvy);
  ntp_gtrack->SetBranchAddress("gvz", &tvz);
  ntp_gtrack->SetBranchAddress("gnhits", &tnhits);
  ntp_gtrack->SetBranchAddress("chisq", &tchisq);
  ntp_gtrack->SetBranchAddress("gprimary", &tprimary);
//ntp_gtrack->SetBranchAddress("bestpurity", &tbestpurity);
  ntp_gtrack->SetBranchAddress("gembed", &tembed);

  //get Branches
  b_tpx = ntp_gtrack->GetBranch("gpx");
  b_tpy = ntp_gtrack->GetBranch("gpy");
  b_tpz = ntp_gtrack->GetBranch("gpz");
  b_tfpx = ntp_gtrack->GetBranch("gfpx");
  b_tfpy = ntp_gtrack->GetBranch("gfpy");
  b_tfpz = ntp_gtrack->GetBranch("gfpz");
  b_tevent = ntp_gtrack->GetBranch("event");
  b_tgtrackid = ntp_gtrack->GetBranch("gtrackID");
  b_tflavor = ntp_gtrack->GetBranch("gflavor");
  b_tvx = ntp_gtrack->GetBranch("gvx");
  b_tvy = ntp_gtrack->GetBranch("gvy");
  b_tvz = ntp_gtrack->GetBranch("gvz");
  b_tnhits = ntp_gtrack->GetBranch("gnhits");
  b_tprimary = ntp_gtrack->GetBranch("gprimary");
//b_tbestpurity = ntp_gtrack->GetBranch("bestpurity");
  b_tembed = ntp_gtrack->GetBranch("gembed");

  // ntp_cluster access variables

cout << "set up ntp_cluster access" << endl;
  
  float_t cevent;
  Float_t hitID;
  Float_t x;
  Float_t y;
  Float_t z;
  Float_t gx;
  Float_t gy;
  Float_t gz;
  float_t layer;
  float_t g4hitID;
  Float_t gtrackID;
  Float_t trackID;
  Float_t gflavor;
  Float_t gpx;
  Float_t gpy;
  Float_t gpz;
  Float_t glast;
  Float_t size;
  Float_t cgembed;

  TBranch* b_cevent;
  TBranch* b_hitID;
  TBranch* b_x;
  TBranch* b_y;
  TBranch* b_z; 
  TBranch* b_gx;
  TBranch* b_gy;
  TBranch* b_gz; 
  TBranch* b_layer;
  TBranch* b_g4hitID;
  TBranch* b_gtrackID;
  TBranch* b_trackID;
  TBranch* b_cgflavor;
  TBranch* b_cgpx;
  TBranch* b_cgpy;
  TBranch* b_cgpz;
  TBranch* b_size;
  TBranch* b_cgembed;
  
  //set branches

  ntp_cluster->SetBranchAddress("event", &cevent);
  ntp_cluster->SetBranchAddress("hitID", &hitID);
  ntp_cluster->SetBranchAddress("x", &x);
  ntp_cluster->SetBranchAddress("y", &y);
  ntp_cluster->SetBranchAddress("z", &z);
  ntp_cluster->SetBranchAddress("gx", &gx);
  ntp_cluster->SetBranchAddress("gy", &gy);
  ntp_cluster->SetBranchAddress("gz", &gz);
  ntp_cluster->SetBranchAddress("layer", &layer);
  ntp_cluster->SetBranchAddress("g4hitID", &g4hitID);
  ntp_cluster->SetBranchAddress("gtrackID", &gtrackID);
  ntp_cluster->SetBranchAddress("trackID", &trackID);
  ntp_cluster->SetBranchAddress("gflavor", &gflavor);
  ntp_cluster->SetBranchAddress("gpx", &gpx);
  ntp_cluster->SetBranchAddress("gpy", &gpy);
  ntp_cluster->SetBranchAddress("gpz", &gpz);
  ntp_cluster->SetBranchAddress("size", &size);
  ntp_cluster->SetBranchAddress("gembed", &cgembed);
  
  //get Branches
  b_cevent = ntp_cluster->GetBranch("event");
  b_hitID = ntp_cluster->GetBranch("hitID");
  b_x = ntp_cluster->GetBranch("x");
  b_y = ntp_cluster->GetBranch("y");
  b_z = ntp_cluster->GetBranch("z");
  b_x = ntp_cluster->GetBranch("gx");
  b_y = ntp_cluster->GetBranch("gy");
  b_z = ntp_cluster->GetBranch("gz");
  b_layer = ntp_cluster->GetBranch("layer");
  b_g4hitID = ntp_cluster->GetBranch("g4hitID");
  b_gtrackID = ntp_cluster->GetBranch("gtrackID");
  b_cgflavor = ntp_cluster->GetBranch("gflavor");
  b_cgpx = ntp_cluster->GetBranch("gpx");
  b_cgpy = ntp_cluster->GetBranch("gpy");
  b_cgpz = ntp_cluster->GetBranch("gpz");
  b_size = ntp_cluster->GetBranch("size");
  b_cgembed = ntp_cluster->GetBranch("gembed");


