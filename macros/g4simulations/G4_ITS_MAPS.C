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
 
  double maps_layer_radius[7] = {23.0, 31.0, 39.0, 194.0, 247.0, 353.0, 405.0};   // mm
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
  maps_cells->set_pixel_x(0.0028);  // 28 microns in cm
  maps_cells->set_pixel_y(0.0028);  
  maps_cells->Verbosity(verbosity);
  se->registerSubsystem(maps_cells);

  return;
}

void Maps_Reco(int verbosity = 0)
{

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();

  // still need to digitize the pixel energy

  // still need to apply live area efficiency to hits
  
  // still need to apply MIP threshoilds to hits

  // still need to make clusters
 
  // still need to reconstruct tracks
 
  // still need to run ghost rejection

  // still need to make track projections

  // still need to run beam spot reco

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

  // still need to implement Maps evaluator

  return;
}
