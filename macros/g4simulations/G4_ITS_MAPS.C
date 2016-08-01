//  macro for realistic ITS MAPS tracker geometry

int Min_maps_layer = 0;
int Max_maps_layer = 6;
 
void MapsInit(int verbosity = 0)
{
  Min_maps_layer = 0;
  Max_maps_layer = 6;
}

double Maps(PHG4Reco* g4Reco, double radius, 
	    const int absorberactive = 0,
	    int verbosity = 0)
{

  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4testbench.so");

  bool overlapcheck = false; // set to true if you want to check for overlaps

  //---------------------------------
  // Inner Cylinder layers for pixels
  //--------------------------------- 
 
  //double maps_layer_radius[7] = {23.0, 31.0, 39.0, 194.0, 247.0, 353.0, 405.0};   // mm
  double maps_layer_radius[7] = {23.635, 31.5, 39.385, 217.6775, 272.07, 363.179, 415.0775};   // mm  - adjusted for closest fit
 // type 1 = inner barrel stave, 2 = middle barrel stave, 3 = outer barrel stave
  int stave_type[7] = {0, 0, 0, 1, 1, 2, 2};
  double max_radius = 0.0;

  for (int ilayer = Min_maps_layer; ilayer <= Max_maps_layer; ilayer++)
    {
      cout << " ilayer = " << ilayer << endl;

      cout << "Create Maps layer " << ilayer  << " with radius " << maps_layer_radius[ilayer] << " stave type " << stave_type[ilayer] << endl;
      PHG4MapsSubsystem  *lyr = new PHG4MapsSubsystem("MAPS", ilayer, stave_type[ilayer]);
      lyr->Verbosity(2);
      lyr->set_nominal_layer_radius(maps_layer_radius[ilayer]);
      // The cell size is used only during pixilization of sensor hits, but it is convemient to set it now 
      lyr->set_pixel_x(0.0020);  // 20 microns in cm
      lyr->set_pixel_z(0.0020);  // 20 microns in cm
      lyr->set_pixel_thickness(0.0018);  // 18 microns in cm
      lyr->SetActive();
      lyr->OverlapCheck(overlapcheck);
      
      cout << "Created Maps layer for layer " << ilayer << " with radius  " << maps_layer_radius[ilayer]  << endl;
      
      g4Reco->registerSubsystem( lyr );      

      cout << "Registered  Maps layer for layer " << ilayer <<  endl;

      max_radius = maps_layer_radius[ilayer];
    }
     return max_radius * 0.1;
}

void Maps_Cells(int verbosity = 0)
{
  // runs the cellularization of the energy deposits (g4hits) 
  // into detector hits (g4cells)

  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  //-----------
  // Maps cells
  //-----------

  // the hits we get now are the entry and exit points in the Maps sensor
  // Later they have to be assigned to Maps pixels.

  PHG4MapsCellReco *maps_cells = new PHG4MapsCellReco("MAPS");
  maps_cells->Verbosity(0);
  se->registerSubsystem(maps_cells);

  return;
}

void Maps_Reco(int verbosity = 0)
{

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  //----------------------------------
  // Digitize the cell energy into ADC
  //----------------------------------
  // defaults to 8-bit ADC with MIP at 0.25% dynamic range
  PHG4SvtxDigitizer* digi = new PHG4SvtxDigitizer();
  digi->Verbosity(0);
  digi->set_adc_scale(0, 255, 1.0e-6); // 1.0 keV / bit
  digi->set_adc_scale(1, 255, 1.0e-6); // 1.0 keV / bit
  digi->set_adc_scale(2, 255, 1.6e-6); // 1.6 keV / bit
  digi->set_adc_scale(3, 255, 1.6e-6); // 1.6 keV / bit
  digi->set_adc_scale(4, 255, 1.6e-6); // 1.6 keV / bit
  digi->set_adc_scale(5, 255, 1.6e-6); // 1.6 keV / bit
  digi->set_adc_scale(6, 255, 1.6e-6); // 1.6 keV / bit
  se->registerSubsystem( digi );

 
  //------------------------------------------
  // Apply Live Area Inefficiency to Hit Cells
  //------------------------------------------
  // defaults to 1.0 (fully active)
  PHG4SvtxDeadArea* deadarea = new PHG4SvtxDeadArea();
  deadarea->Verbosity(0);
  deadarea->set_hit_efficiency(0,0.99); // Leo says use 1% inefficiency
  deadarea->set_hit_efficiency(1,0.99);
  deadarea->set_hit_efficiency(2,0.99);
  deadarea->set_hit_efficiency(3,0.99);
  deadarea->set_hit_efficiency(4,0.99);
  deadarea->set_hit_efficiency(5,0.99);
  deadarea->set_hit_efficiency(6,0.99);
  se->registerSubsystem( deadarea );


  //----------------------------------
  // Apply MIP thresholds to Hit Cells
  //----------------------------------
  PHG4SvtxThresholds* thresholds = new PHG4SvtxThresholds();
  thresholds->Verbosity(0);
  thresholds->set_threshold(0,0.25);
  thresholds->set_threshold(1,0.25);
  thresholds->set_threshold(2,0.25);
  thresholds->set_threshold(3,0.25);
  thresholds->set_threshold(4,0.25);
  thresholds->set_threshold(5,0.25);
  thresholds->set_threshold(6,0.25);
  //thresholds->set_use_thickness_mip(0, true);
  se->registerSubsystem( thresholds );


  
 //---------------------
  // Make SVTX clusters
  //---------------------
  PHG4SvtxClusterizer* clusterizer = new PHG4SvtxClusterizer();
  clusterizer->Verbosity(0);
  clusterizer->set_threshold(0.33);
  se->registerSubsystem( clusterizer );

  //---------------------
  // Track reconstruction
  //---------------------
  PHG4HoughTransform* hough = new PHG4HoughTransform(7,7);
  hough->set_mag_field(1.4);
  hough->Verbosity(0);
  // ALICE ITS upgrade values for total thickness in X_0
  hough->set_material(0, 0.003);
  hough->set_material(1, 0.003);
  hough->set_material(2, 0.003);
  hough->set_material(3, 0.008);
  hough->set_material(4, 0.008);
  hough->set_material(5, 0.008);
  hough->set_material(6, 0.008);
  hough->setPtRescaleFactor(0.9972);
  hough->set_chi2_cut_init(5.0);
  //hough->set_chi2_cut_fast(60.0,0.0,100.0); // 10.0, 50.0, 75.0
  hough->set_chi2_cut_fast(10.0,50.0,75.0); // 10.0, 50.0, 75.0
  hough->set_chi2_cut_full(5.0);
  hough->set_ca_chi2_cut(5.0);
  hough->setMaxClusterError(3.0);
  hough->setRejectGhosts(false);
  hough->setRemoveHits(false);
  hough->setCutOnDCA(true);
  se->registerSubsystem( hough );

  //---------------------
  // Ghost rejection
  //---------------------
  PHG4TrackGhostRejection* rejection = new PHG4TrackGhostRejection(7);
  rejection->Verbosity(0);
  rejection->set_max_shared_hits(3);
  se->registerSubsystem( rejection );

  //------------------
  // Track Projections
  //------------------
  PHG4SvtxTrackProjection* projection = new PHG4SvtxTrackProjection();
  projection->Verbosity(verbosity);
  se->registerSubsystem( projection );

  //----------------------
  // Beam Spot Calculation
  //----------------------
  PHG4SvtxBeamSpotReco* beamspot = new PHG4SvtxBeamSpotReco();
  beamspot->Verbosity(verbosity);
  se->registerSubsystem( beamspot );

  return;
}

void Maps_Eval(std::string outputfile, int verbosity = 0)
{
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libg4eval.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  SubsysReco* eval = new SvtxEvaluator("SvtxEVALUATOR", outputfile.c_str());
  eval->Verbosity(0);
  //eval->scan_for_embedded(false); 
  se->registerSubsystem( eval );
  
  return;
}
