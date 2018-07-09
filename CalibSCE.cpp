#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TView.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TVector3.h>
#include <TPrincipal.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

#include "normVec.hpp"

#include "/uboone/app/users/joelam/SpaceChargeStudy/CalibSCE/Eigen/Dense.h"

using namespace std;

//const Char_t *inputFileLaser = "data/laserDataSCE.root";
const Char_t *inputFileLaser = "data/laserDataSCE_NEW.root";
//const Char_t *inputFileLaser = "data/laserScan_data.root";
//const Char_t *inputFileLaser = "data/laserScan_MC.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small_NEW.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small_NEW_withMCS.root";
//const Char_t *inputFileCosmic = "data/mcsample_small_combined.root";
//const Char_t *inputFileCosmic = "data/oldMCsample_50000events.root";
//const Char_t *inputFileCosmic = "data/newMCsample_2Mevents.root";
//const Char_t *inputFileCosmic = "data/MC_Cosmics.root";
const Char_t *inputFileCosmic = "/uboone/data/users/mrmooney/MC_Cosmics.root";
//const Char_t *inputFileCosmic = "data/Data_Run1_EXTBNB.root";
//const Char_t *inputFileCosmic = "data/first_set_of_EXTBNB_production.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_ProtoDUNESP.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_ProtoDUNESP_withMCS.root";

TFile* outputFile = new TFile("output.root","RECREATE");

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 3.6;
//const Double_t Ly = 6.0;
//const Double_t Lz = 7.2;

const Bool_t isMC = true;
const Bool_t doBulk = true;
const Bool_t isSCEon = true;

const Double_t relAngleCut = 20.0;
const Double_t maxXdist = 0.05;
const Double_t maxYdist = 0.20;
const Double_t maxZdist = 0.20;

Int_t minInputTrackNum = 0;
Int_t maxInputTrackNum = 1000000000;
//Int_t maxInputTrackNum = 1150000;

//const Int_t maxCosmicTracks = -1;
const Int_t maxCosmicTracks = 10000000;
//const Int_t maxCosmicTracks = 1000000;
//const Int_t maxCosmicTracks = 30000;
//const Int_t maxCosmicTracks = 10000;

Double_t minTrackMCS_anode;
Double_t minTrackMCS_cathode;
Double_t minTrackMCS_crossing;

Int_t nCalibDivisions = 25; // MICROBOONE
//Int_t nCalibDivisions = 18; // PROTODUNE-SP

const Double_t piVal = 3.14159265;

Int_t nCalibDivisions_x;
Int_t nCalibDivisions_y;
Int_t nCalibDivisions_z;

Double_t calibWeight[101][101][401];
Double_t calibDeltaX[101][101][401];
Double_t calibDeltaY[101][101][401];
Double_t calibDeltaZ[101][101][401];

Double_t trueDeltaX[101][101][401];
Double_t trueDeltaY[101][101][401];
Double_t trueDeltaZ[101][101][401];

Double_t trueFwdDeltaX[101][101][401];
Double_t trueFwdDeltaY[101][101][401];
Double_t trueFwdDeltaZ[101][101][401];

Double_t calibTopDeltaX[101][401];
Double_t calibTopDeltaY[101][401];
Double_t calibTopDeltaZ[101][401];
Double_t calibTopWeight[101][401];

Double_t calibBottomDeltaX[101][401];
Double_t calibBottomDeltaY[101][401];
Double_t calibBottomDeltaZ[101][401];
Double_t calibBottomWeight[101][401];

Double_t calibDownstreamDeltaX[101][101];
Double_t calibDownstreamDeltaY[101][101];
Double_t calibDownstreamDeltaZ[101][101];
Double_t calibDownstreamWeight[101][101];

Double_t calibUpstreamDeltaX[101][101];
Double_t calibUpstreamDeltaY[101][101];
Double_t calibUpstreamDeltaZ[101][101];
Double_t calibUpstreamWeight[101][101];

Double_t calibCathodeDeltaX[101][401];
Double_t calibCathodeDeltaY[101][401];
Double_t calibCathodeDeltaZ[101][401];
Double_t calibCathodeWeight[101][401];

struct elecInfo
{
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t t;
  Double_t x_mod;
  Double_t y_mod;
  Double_t z_mod;
  Double_t t_mod;
  Int_t fate;
};

struct trackInfo
{
  Int_t pdgID;
  Double_t energy;
  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Double_t theta;
  Double_t phi;
  vector<elecInfo> electrons;
};

struct calibTrackInfo
{
  Double_t x0_calib;
  Double_t y0_calib;
  Double_t z0_calib;
  Double_t x1_calib;
  Double_t y1_calib;
  Double_t z1_calib;
  Double_t theta_calib;
  Double_t phi_calib;
  vector<Double_t> DxVec;
  vector<Double_t> DyVec;
  vector<Double_t> DzVec;
  Bool_t calibFlag;
  trackInfo track;
};

struct Point {
  float x;
  float y;
  float z;
};

struct PCAResults {
  TVector3 centroid;
  pair<TVector3,TVector3> endPoints;
  float length;
  TVector3 eVals;
  vector<TVector3> eVecs;
};

typedef vector<Point> PointCloud;

Double_t doCoordTransformX(const Double_t inputX);
Double_t doCoordTransformY(const Double_t inputY);
Double_t doCoordTransformZ(const Double_t inputZ);
vector<Double_t> getParabolaParameters(const vector<elecInfo> &parabola_points_track);
vector<Double_t> findClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
vector<Double_t> findDistortedClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
vector<trackInfo> getLArSoftTrackSet(Int_t inputType);
vector<trackInfo> getTrackSet(Int_t mapType);
vector<calibTrackInfo> makeCalibTracks(const vector<trackInfo> &tracks);
void doLaserLaserCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor);
void doLaserCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode);
void doCosmicCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode);
void updateCalibTrueTrack(calibTrackInfo &calibTrack);
void updateAllCalibTrueTracks(vector<calibTrackInfo> &calibTracks, Int_t iterNum);
void doCalibration(const vector<trackInfo> &laserTracks, const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode, Int_t numIterations);
void saveTrackInfo(const vector<trackInfo> &tracks);
void loadTruthMap();
Double_t getTruthOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp);
void loadTruthFwdMap();
Double_t getTruthFwdOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp);
void doCalibFaces(const vector<trackInfo> &cosmicTracks, Int_t minTrackPoints, Int_t numTrackSegPoints);
PCAResults DoPCA(const PointCloud &points);

Int_t main(Int_t argc, Char_t** argv)
{
  TStopwatch timer;
  timer.Start();

  if(argc > 1) {
    minInputTrackNum = atoi(argv[1]);
    maxInputTrackNum = atoi(argv[2]);
  }
  
  nCalibDivisions_x = nCalibDivisions;
  nCalibDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)nCalibDivisions));
  nCalibDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)nCalibDivisions));

  if (doBulk == false) {
    minTrackMCS_anode = 0.0;
    minTrackMCS_cathode = 0.0;
    minTrackMCS_crossing = 0.0;
  }
  else if (isMC == true) {
    minTrackMCS_anode = 3.3;
    minTrackMCS_cathode = 1.7;
    minTrackMCS_crossing = 1.4;
  }
  else {
    minTrackMCS_anode = 1.5;
    minTrackMCS_cathode = 1.1;
    minTrackMCS_crossing = 0.0;
  }

  loadTruthMap();
  loadTruthFwdMap();
  
  //vector<trackInfo> laserTracks = getTrackSet(1);
  vector<trackInfo> laserTracks = getLArSoftTrackSet(1);
  //vector<trackInfo> cosmicTracks = getTrackSet(2);
  vector<trackInfo> cosmicTracks = getLArSoftTrackSet(2);

  saveTrackInfo(cosmicTracks);

  doCalibFaces(cosmicTracks,50,15); // was 50,15
  ////doCalibration(laserTracks,cosmicTracks,0.05,3,1,1);
  ////doCalibration(laserTracks,cosmicTracks,0.02,3,1,1);
  //////doCalibration(laserTracks,cosmicTracks,0.01,3,0,2);
  //doCalibration(laserTracks,cosmicTracks,0.01,3,1,1); // Nominal Configuration
  ////doCalibration(laserTracks,cosmicTracks,0.01,1,1,1); // Nominal Configuration #2
  
  timer.Stop();
  cout << "Calibration Time:  " << timer.CpuTime() << " sec." << endl;
  
  outputFile->Write();
  outputFile->Close();
  
  return 0;
}

Double_t doCoordTransformX(const Double_t inputX)
{
  Double_t outputX;

  //outputX = Lx - (Lx/2.56)*inputX/100.0;
  if(isMC) {
    outputX = Lx - (Lx/2.5524)*inputX/100.0;
  }
  else {
    outputX = Lx - (Lx/2.58)*inputX/100.0;
  }

  return outputX;
}

Double_t doCoordTransformY(const Double_t inputY)
{
  Double_t outputY;

  if(isMC) {
    outputY = (Ly/(1.173+1.154))*(inputY+115.4)/100.0;
  }
  else {
    outputY = (Ly/(1.170+1.151))*(inputY+115.1)/100.0;
  }

  return outputY;
}

Double_t doCoordTransformZ(const Double_t inputZ)
{
  Double_t outputZ;

  if(isMC) {
    outputZ = (Lz/(10.368-0.003))*(inputZ-0.3)/100.0;
  }
  else {
    outputZ = (Lz/(10.365+0.007))*(inputZ+0.7)/100.0;
  }
  
  return outputZ;
}

vector<Double_t> getParabolaParameters(const vector<elecInfo> &parabola_points_track)
{
  if (parabola_points_track.size() < 3) cout << "Less than three points provided for the parameters of parabola." << endl;
 
  Double_t x_middle = parabola_points_track.at(1).x_mod;
  Double_t y_middle = parabola_points_track.at(1).y_mod;
  Double_t z_middle = parabola_points_track.at(1).z_mod;
 
  //first/linear transformation
  Double_t x_0 = parabola_points_track.at(0).x_mod-x_middle;
  Double_t x_1 = parabola_points_track.at(1).x_mod-x_middle;
  Double_t x_2 = parabola_points_track.at(2).x_mod-x_middle;
 
  Double_t y_0 = parabola_points_track.at(0).y_mod-y_middle;
  Double_t y_1 = parabola_points_track.at(1).y_mod-y_middle;
  Double_t y_2 = parabola_points_track.at(2).y_mod-y_middle;
 
  Double_t z_0 = parabola_points_track.at(0).z_mod;
  Double_t z_1 = parabola_points_track.at(1).z_mod;
  Double_t z_2 = parabola_points_track.at(2).z_mod;
 
  //angle of rotation where y_0_2 == y_2_2
  Double_t phi = atan((y_2-y_0)/(x_0-x_2));
 
  //preform the second transformation
  Double_t x_0_2 = x_0*cos(phi)-y_0*sin(phi);
  Double_t x_1_2 = x_1*cos(phi)-y_1*sin(phi);
  Double_t x_2_2 = x_2*cos(phi)-y_2*sin(phi);
 
  Double_t y_0_2 = x_0*sin(phi)+y_0*cos(phi);
  Double_t y_1_2 = x_1*sin(phi)+y_1*cos(phi);
  Double_t y_2_2 = x_2*sin(phi)+y_2*cos(phi);
 
  //since x_1_2 = 0 and y_1_2 = 0, because that's the middle point, c = 0
  Double_t a = ((y_0_2*x_2_2)-(y_2_2*x_0_2))/(x_0_2*x_2_2*(x_0_2-x_2_2));
  Double_t b = (y_0_2-(a*pow(x_0_2,2)))/x_0_2;
 
  //find the plane defined by the three points without any transformations: z = d*x + e*y + f
  //any one of the three points dotted with the normal vector will give the plane: P_1*n = 0
 
  //comps. for a vector connecting point 2 to point 1
  Double_t p1p2_x = parabola_points_track.at(1).x_mod-parabola_points_track.at(0).x_mod;
  Double_t p1p2_y = parabola_points_track.at(1).y_mod-parabola_points_track.at(0).y_mod;
  Double_t p1p2_z = parabola_points_track.at(1).z_mod-parabola_points_track.at(0).z_mod;
 
  //comps. for a vector connecting point 3 to point 1
  Double_t p1p3_x = parabola_points_track.at(2).x_mod-parabola_points_track.at(0).x_mod;
  Double_t p1p3_y = parabola_points_track.at(2).y_mod-parabola_points_track.at(0).y_mod;
  Double_t p1p3_z = parabola_points_track.at(2).z_mod-parabola_points_track.at(0).z_mod;
 
  // normal vector, n, = P1P2xP1P3
  Double_t norm_x = (p1p2_y*p1p3_z)-(p1p3_y*p1p2_z);
  Double_t norm_y = (p1p3_x*p1p2_z)-(p1p2_x*p1p3_z);
  Double_t norm_z = (p1p2_x*p1p3_y)-(p1p3_x*p1p2_y);
 
  // z = d*x + e*y + f
  Double_t d = -norm_x/norm_z;
  Double_t e = -norm_y/norm_z;
  Double_t f = ((parabola_points_track.at(0).x_mod*norm_x)+(parabola_points_track.at(0).y_mod*norm_y)+(parabola_points_track.at(0).z_mod*norm_z))/norm_z;
 
  vector<Double_t> return_vector;
  return_vector.push_back(a);
  return_vector.push_back(b);
  return_vector.push_back(d);
  return_vector.push_back(e);
  return_vector.push_back(f);
  return_vector.push_back(phi);
  return_vector.push_back(x_middle);
  return_vector.push_back(y_middle);
  return return_vector;
}

vector<Double_t> findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB)
{ 
  //   (xA,yA,zA)+t(xA_step,yA_step,zA_step)
  Double_t xA = calibTrackA.x0_calib, yA = calibTrackA.y0_calib, zA = calibTrackA.z0_calib;
  Double_t xB = calibTrackB.x0_calib, yB = calibTrackB.y0_calib, zB = calibTrackB.z0_calib;

  Double_t xA_step = -1.*sin(calibTrackA.theta_calib)*sin(calibTrackA.phi_calib);
  Double_t yA_step = cos(calibTrackA.theta_calib);
  Double_t zA_step = sin(calibTrackA.theta_calib)*cos(calibTrackA.phi_calib);

  Double_t xB_step = -1.*sin(calibTrackB.theta_calib)*sin(calibTrackB.phi_calib);
  Double_t yB_step = cos(calibTrackB.theta_calib);
  Double_t zB_step = sin(calibTrackB.theta_calib)*cos(calibTrackB.phi_calib);

  //perpendicular line between the two tracks
  Double_t x_prep = (yA_step*zB_step)-(yB_step*zA_step);
  Double_t y_prep = (xB_step*zA_step)-(xA_step*zB_step);
  Double_t z_prep = (xA_step*yB_step)-(xB_step*yA_step);
 
  // if cross product is zero then the lines are parallel so return distance = -2
  if (x_prep == 0 && y_prep == 0 && z_prep == 0) {
    vector<Double_t> return_vector;
    return_vector.push_back(-2.);
    return_vector.push_back((xA+xB)/2.);
    return_vector.push_back((yA+yB)/2.);
    return_vector.push_back((zA+zB)/2.);   
    return return_vector;
  }
  //normalize the perpendicular line
  Double_t mag_prep = sqrt(pow(x_prep,2)+pow(y_prep,2)+pow(z_prep,2));
  Double_t x_prep_norm = x_prep / mag_prep;
  Double_t y_prep_norm = y_prep / mag_prep;
  Double_t z_prep_norm = z_prep / mag_prep;
 
  //defined to make the math simplier
  Double_t a = y_prep_norm*(xA-xB);
  Double_t b = x_prep_norm*(yA-yB);
  Double_t c = z_prep_norm*(xA-xB);
  Double_t d = x_prep_norm*(zA-zB);
 
  Double_t g = y_prep_norm*xA_step;
  Double_t h = y_prep_norm*xB_step;
  Double_t i = x_prep_norm*yA_step;
  Double_t j = x_prep_norm*yB_step;
  Double_t k = z_prep_norm*xA_step;
  Double_t l = z_prep_norm*xB_step;
  Double_t m = x_prep_norm*zA_step;
  Double_t n = x_prep_norm*zB_step;
 
  Double_t chi = (l-n)/(h-j);
 
  //alpha: "t" for the first line //beta:  "t" for the second line
  Double_t alpha = (chi*(a-b)-c+d)/(k-m-(chi*(g-i)));
  Double_t beta = (c-d+alpha*(k-m))/(l-n);

  Double_t cpa_xA = xA+alpha*xA_step;
  Double_t cpa_yA = yA+alpha*yA_step;
  Double_t cpa_zA = zA+alpha*zA_step;

  Double_t cpa_xB = xB+beta*xB_step;
  Double_t cpa_yB = yB+beta*yB_step;
  Double_t cpa_zB = zB+beta*zB_step;
 
  //distance between the closest points on the lines\tracks
  //Double_t distance = sqrt(pow((xA+alpha*xA_step)-(xB+beta*xB_step),2)+pow((yA+alpha*yA_step)-(yB+beta*yB_step),2)+pow((zA+alpha*zA_step)-(zB+beta*zB_step),2));
  Double_t distance = sqrt(pow(cpa_xA-cpa_xB,2)+pow(cpa_yA-cpa_yB,2)+pow(cpa_zA-cpa_zB,2));

  //midpoint between points of closest approach on both lines
  Double_t x_mid = (cpa_xA+cpa_xB)/2.;
  Double_t y_mid = (cpa_yA+cpa_yB)/2.;
  Double_t z_mid = (cpa_zA+cpa_zB)/2.;
 
  //check to see if this point is outside the detector
  //  TODO: < or <=, Is 0 or Lx "outside"?
  //if (x_mid > Lx || x_mid < 0 || y_mid > Ly || y_mid < 0 || z_mid > Lz || z_mid < 0) distance = -1;
  if (cpa_xA > Lx || cpa_xA < 0 || cpa_yA > Ly || cpa_yA < 0 || cpa_zA > Lz || cpa_zA < 0 ||
      cpa_xB > Lx || cpa_xB < 0 || cpa_yB > Ly || cpa_yB < 0 || cpa_zB > Lz || cpa_zB < 0) distance = -1;
 
  vector<Double_t> return_vector;

  return_vector.push_back(distance);
  return_vector.push_back(x_mid);
  return_vector.push_back(y_mid);
  return_vector.push_back(z_mid);

  ////// NEW 12/5/2017 //////
  return_vector.push_back(cpa_xB-cpa_xA);
  return_vector.push_back(cpa_yB-cpa_yA);
  return_vector.push_back(cpa_zB-cpa_zA);
  ///////////////////////////
  
  return return_vector;
}

vector<Double_t> findDistortedClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB)
{
  vector<Double_t> return_vector;
 
  if (calibTrackA.track.electrons.size() < 3) {
    cout << "Less than three track points are inside the detector for calibTrackA." << endl;
    return return_vector;
  }
  if (calibTrackB.track.electrons.size() < 3) {
    cout << "Less than three track points are inside the detector for calibTrackB." << endl;
    return return_vector;
  }

  Double_t min_distance = -1.;//-1 means no value has yet been set; for the first iteration
  Double_t min_a, min_b;
 
  for(int iter_a = 0; iter_a < calibTrackA.track.electrons.size(); iter_a++)
    for(int iter_b = 0; iter_b < calibTrackB.track.electrons.size(); iter_b++) {
      Double_t distance = sqrt(pow(calibTrackA.track.electrons.at(iter_a).x_mod-calibTrackB.track.electrons.at(iter_b).x_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).y_mod-calibTrackB.track.electrons.at(iter_b).y_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).z_mod-calibTrackB.track.electrons.at(iter_b).z_mod,2));
      if (min_distance == -1 || distance < min_distance) {
	min_distance = distance;
	min_a = iter_a; min_b = iter_b;
      }
    }
    
  //TODO if min_distance is zero then no need to fit a parabola

  vector<elecInfo> parabola_points_calibTrackA, parabola_points_calibTrackB;
 
  if (min_a != 0) parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a-1));
  else parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a+2));
  parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a));
  if (min_a != calibTrackA.track.electrons.size() - 1) parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a+1));
  else parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a-2));
 
  if (min_b != 0) parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b-1));
  else parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b+2));
  parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b));
  if (min_b != calibTrackB.track.electrons.size()-1) parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b+1));
  else parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b-2));
 
  vector<Double_t> parabolaParameters_calibTrackA = getParabolaParameters(parabola_points_calibTrackA);
  vector<Double_t> parabolaParameters_calibTrackB = getParabolaParameters(parabola_points_calibTrackB);
 
  Double_t aA = parabolaParameters_calibTrackA.at(0);
  Double_t bA = parabolaParameters_calibTrackA.at(1);
  Double_t dA = parabolaParameters_calibTrackA.at(2);
  Double_t eA = parabolaParameters_calibTrackA.at(3);
  Double_t fA = parabolaParameters_calibTrackA.at(4);
  Double_t phiA = parabolaParameters_calibTrackA.at(5);
  Double_t x_midA = parabolaParameters_calibTrackA.at(6);
  Double_t y_midA = parabolaParameters_calibTrackA.at(7);
 
  Double_t aB = parabolaParameters_calibTrackB.at(0);
  Double_t bB = parabolaParameters_calibTrackB.at(1);
  Double_t dB = parabolaParameters_calibTrackB.at(2);
  Double_t eB = parabolaParameters_calibTrackB.at(3);
  Double_t fB = parabolaParameters_calibTrackB.at(4);
  Double_t phiB = parabolaParameters_calibTrackB.at(5);
  Double_t x_midB = parabolaParameters_calibTrackB.at(6);
  Double_t y_midB = parabolaParameters_calibTrackB.at(7);
 
  Double_t stepSize_x2 = 2.5*pow(10,-4); //0.25mm
 
  Double_t min_distance_parabola = -1.;
 
  Double_t min_Ax = 0., min_Ay = 0., min_Az = 0., min_Bx = 0., min_By = 0., min_Bz = 0.;
  if (min_a != 0 && min_a != calibTrackA.track.electrons.size()-1 && min_b != 0 && min_b != calibTrackB.track.electrons.size()-1) {
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    Double_t A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    Double_t A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    Double_t B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    Double_t B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
   
    //starting and ending points for the loop //values determined below
    Double_t A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
    if (A_x2_0 < A_x2_2) {
      A_x2_start = A_x2_0;
      A_x2_end = A_x2_2;
    } else {
      A_x2_start = A_x2_2;
      A_x2_end = A_x2_0;
    }
    if (B_x2_0 < B_x2_2) {
      B_x2_start = B_x2_0;
      B_x2_end = B_x2_2;
    } else {
      B_x2_start = B_x2_2;
      B_x2_end = B_x2_0;
    }
   
    for (Double_t A_x2 = A_x2_start; A_x2 < A_x2_end;  A_x2 += stepSize_x2) {
      Double_t A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      //get x,y,z
      Double_t A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      Double_t A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      Double_t A_z = dA*A_x+eA*A_y+fA;
      for (Double_t B_x2 = B_x2_start; B_x2 < B_x2_end; B_x2 += stepSize_x2) {
	Double_t B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	Double_t B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	Double_t B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	Double_t B_z = dB*B_x+eB*B_y+fB;
   
	Double_t distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;
     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }
  else {
    Double_t stepSize_x2_A = stepSize_x2;
    Double_t stepSize_x2_B = stepSize_x2;
   
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    Double_t A_x2_0 = min_a == 0 ? -1. : (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    Double_t A_x2_2 = min_a == calibTrackA.track.electrons.size()-1 ? -1. : (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    Double_t B_x2_0 = min_b == 0 ? -1. : (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    Double_t B_x2_2 = min_b == calibTrackB.track.electrons.size()-1 ? -1. : (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);

    Int_t nearBoundary_A = 1, nearBoundary_B = 1;
   
    Double_t A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
   
    //If the first point is set to be -1 then start the loop form the final point and move towards the first point until the detector boundary is reached.
    if (A_x2_0 == -1.) A_x2_start = A_x2_2;
    else if (A_x2_2 == -1.) A_x2_start = A_x2_0;
    else {
      nearBoundary_A = 0;
     
      A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
      A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
      if (A_x2_0 < A_x2_2) {
	A_x2_start = A_x2_0;
	A_x2_end = A_x2_2;
      } else {
	A_x2_start = A_x2_2;
	A_x2_end = A_x2_0;
      }
    }
    if (B_x2_0 == -1.) B_x2_start = B_x2_2;
    else if (B_x2_2 == -1.) B_x2_start = B_x2_0;
    else {
      nearBoundary_B = 0;
     
      B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
      B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
      if (B_x2_0 < B_x2_2) {
	B_x2_start = B_x2_0;
	B_x2_end = B_x2_2;
      } else {
	B_x2_start = B_x2_2;
	B_x2_end = B_x2_0;
      }
    }

    if (nearBoundary_A == 1) if (A_x2_start > 0) stepSize_x2_A *= -1.;
    if (nearBoundary_B == 1) if (B_x2_start > 0) stepSize_x2_B *= -1.;

    Double_t A_x = 0., A_y = 0., A_z = 0.;
    Int_t firstIter_A = 1; //1 means it's the first iteration
    for (Double_t A_x2 = A_x2_start; nearBoundary_A*((A_x > 0 && A_x < Lx && A_y > 0 && A_y < Ly && A_z > 0 && A_z < Lz)+firstIter_A)+!nearBoundary_A*(A_x2 < A_x2_end); A_x2 += stepSize_x2_A) {
      firstIter_A = 0;
      
      Double_t A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      A_z = dA*A_x+eA*A_y+fA;
     
      Double_t B_y = 0., B_x = 0., B_z = 0.;
      Int_t firstIter_B = 1;
      for (Double_t B_x2 = B_x2_start; nearBoundary_B*((B_x > 0 && B_x < Lx && B_y > 0 && B_y < Ly && B_z > 0 && B_z < Lz)+firstIter_B)+!nearBoundary_B*(B_x2 < B_x2_end); B_x2 += stepSize_x2_B) {
	firstIter_B = 0;
   
	Double_t B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	B_z = dB*B_x+eB*B_y+fB;

	Double_t distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }

  //midpoint between points of closest approach on both lines
  Double_t x_mid = (min_Ax+min_Bx)/2.;
  Double_t y_mid = (min_Ay+min_By)/2.;
  Double_t z_mid = (min_Az+min_Bz)/2.;

  if(isnan(min_distance_parabola)) {
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);    
  }
  else {
    return_vector.push_back(min_distance_parabola);
    return_vector.push_back(x_mid);
    return_vector.push_back(y_mid);
    return_vector.push_back(z_mid);
  }

  ////// NEW 12/5/2017 //////
  return_vector.push_back(min_Bx-min_Ax);
  return_vector.push_back(min_By-min_Ay);
  return_vector.push_back(min_Bz-min_Az);
  ///////////////////////////

  return return_vector;
}

vector<trackInfo> getLArSoftTrackSet(Int_t inputType)
{
  vector<trackInfo> tracks;

  TFile* mapFile;

  if(inputType == 1)
  {
    mapFile = new TFile(inputFileLaser,"READ");

    TTreeReader readerLasers("lasers", mapFile);
    TTreeReader readerTracks("tracks", mapFile);

    TTreeReaderValue<TVector3> track_start(readerLasers, "entry");
    TTreeReaderValue<TVector3> track_end(readerLasers, "exit");
    TTreeReaderArray<TVector3> track_points(readerTracks, "track");

    Int_t nTracks = 0;
    
    while (readerLasers.Next())
    {
      readerTracks.Next();
    
      nTracks++;
    
      Double_t x0, y0, z0;
      Double_t x1, y1, z1;
    
      if ((*track_start).y() > (*track_end).y()) {
        x0 = doCoordTransformX((*track_start).x());
        y0 = doCoordTransformY((*track_start).y());
        z0 = doCoordTransformZ((*track_start).z());
    
        x1 = doCoordTransformX((*track_end).x());
        y1 = doCoordTransformY((*track_end).y());
        z1 = doCoordTransformZ((*track_end).z());
      }
      else {
        x0 = doCoordTransformX((*track_end).x());
        y0 = doCoordTransformY((*track_end).y());
        z0 = doCoordTransformZ((*track_end).z());
    
        x1 = doCoordTransformX((*track_start).x());
        y1 = doCoordTransformY((*track_start).y());
        z1 = doCoordTransformZ((*track_start).z());
      }
      
      Double_t trackLength = sqrt(pow(x0-x1,2.0)+pow(y0-y1,2.0)+pow(z0-z1,2.0));
      
      Double_t theta = acos((y1-y0)/trackLength);
      Double_t phi = acos((z1-z0)/(trackLength*sin(theta)));
      if (x1 > x0) {
        phi = -1.0*fabs(phi);
      }
      else {
        phi = fabs(phi);
      }
    
      trackInfo track;
      elecInfo electron;
    
      track.energy = 1000000.0;
      track.pdgID = 22;
      track.x0 = x0;
      track.y0 = y0;
      track.z0 = z0;
      track.x1 = x1;
      track.y1 = y1;
      track.z1 = z1;
      track.theta = theta;
      track.phi = phi;
      track.electrons.clear();
    
      for(Int_t j = 0; j < track_points.GetSize(); j++)
      {
        if (((j % 10) != 0) && (j != track_points.GetSize()-1)) continue;
        
        electron.x = -1.0; // Dummy (currently don't save this info in file)
        electron.y = -1.0; // Dummy (currently don't save this info in file)
        electron.z = -1.0; // Dummy (currently don't save this info in file)
        electron.t = -1.0; // Dummy (currently don't save this info in file)
        electron.x_mod = min(max(0.0,doCoordTransformX(track_points[j].x())),Lx);
        electron.y_mod = min(max(0.0,doCoordTransformY(track_points[j].y())),Ly);
      	electron.z_mod = min(max(0.0,doCoordTransformZ(track_points[j].z())),Lz);
        electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
        electron.fate = 0; // Dummy (currently don't save this info in file)

        track.electrons.push_back(electron);
      }
      
      tracks.push_back(track);
    }

    return tracks;
  }
  else if(inputType == 2)
  {
    mapFile = new TFile(inputFileCosmic,"READ");
    
    char* filestring = (char*) "";
    char* filestring2 = (char*) "";
    if (isSCEon == false) {
      filestring = (char*) "_copy";
      filestring2 = (char*) "t0ana/";
    }
    
    TTreeReader reader(Form("%sSCEtree%s",filestring2,filestring), mapFile);
    TTreeReaderValue<Double_t> track_startX(reader, Form("track_startX%s",filestring));
    TTreeReaderValue<Double_t> track_startY(reader, Form("track_startY%s",filestring));
    TTreeReaderValue<Double_t> track_startZ(reader, Form("track_startZ%s",filestring));
    TTreeReaderValue<Double_t> track_endX(reader, Form("track_endX%s",filestring));
    TTreeReaderValue<Double_t> track_endY(reader, Form("track_endY%s",filestring));
    TTreeReaderValue<Double_t> track_endZ(reader, Form("track_endZ%s",filestring));
    TTreeReaderValue<Int_t> nPoints(reader, Form("track_nPoints%s",filestring));
    TTreeReaderArray<Double_t> pointX(reader, Form("track_pointX%s",filestring));
    TTreeReaderArray<Double_t> pointY(reader, Form("track_pointY%s",filestring));
    TTreeReaderArray<Double_t> pointZ(reader, Form("track_pointZ%s",filestring));
    TTreeReaderValue<Double_t> track_MCS(reader, Form("track_MCS_measurement%s",filestring));
    TTreeReaderValue<Double_t> track_t0(reader, Form("track_t0%s",filestring));
    
    trackInfo track;
    elecInfo electron;

    TH1F *Zhist_1s = new TH1F("Zhist_1s","",50,0.0,10.0);
    TH1F *Zhist_2s = new TH1F("Zhist_2s","",50,0.0,10.0);
    TH1F *Zhist_3s = new TH1F("Zhist_3s","",50,0.0,10.0);
    TH1F *Zhist_4s = new TH1F("Zhist_4s","",50,0.0,10.0);
    TH1F *Zhist_5s = new TH1F("Zhist_5s","",50,0.0,10.0);
    TH1F *Zhist_6s = new TH1F("Zhist_6s","",50,0.0,10.0);
    TH1F *Zhist_7s = new TH1F("Zhist_7s","",50,0.0,10.0);

    TH1F *Zhist_1e = new TH1F("Zhist_1e","",50,0.0,10.0);
    TH1F *Zhist_2e = new TH1F("Zhist_2e","",50,0.0,10.0);
    TH1F *Zhist_3e = new TH1F("Zhist_3e","",50,0.0,10.0);
    TH1F *Zhist_4e = new TH1F("Zhist_4e","",50,0.0,10.0);
    TH1F *Zhist_5e = new TH1F("Zhist_5e","",50,0.0,10.0);
    TH1F *Zhist_6e = new TH1F("Zhist_6e","",50,0.0,10.0);
    TH1F *Zhist_7e = new TH1F("Zhist_7e","",50,0.0,10.0);

    TF1 *weightFunc;
    if (isMC == true) {
      weightFunc = new TF1("weightFunc","(1.0/75.0)*(x-5.0)^2 + (2.0/3.0)",0.0,10.0);
    }
    else {
      weightFunc = new TF1("weightFunc","(2.0/75.0)*(x-5.0)^2 + (1.0/3.0)",0.0,10.0);
    }
    
    TRandom3 *rand = new TRandom3(1);

    Int_t inputTrackNum = -1;
    Int_t nTracks = 0;
    int containedTracks = 0;
    while (reader.Next())
    {
      inputTrackNum++;
      if ((inputTrackNum < minInputTrackNum) || (inputTrackNum > maxInputTrackNum)) continue;
      
      if ((maxCosmicTracks != -1) && (nTracks >= maxCosmicTracks)) continue;
      
      if (*nPoints < 3) continue;
      
      Double_t xS, yS, zS;
      Double_t xE, yE, zE;
      Double_t x0, y0, z0;
      Double_t x1, y1, z1;
      
      //Double_t x_offset = 0.002564*(*track_t0); // TEMP OFFSET, CHRIS WILL FIX BUG SOON
      Double_t x_offset = 0.0;
      
      if (*track_startY > *track_endY) {
        xS = doCoordTransformX(*track_startX + x_offset);
        yS = doCoordTransformY(*track_startY);
        zS = doCoordTransformZ(*track_startZ);
    
        xE = doCoordTransformX(*track_endX + x_offset);
        yE = doCoordTransformY(*track_endY);
        zE = doCoordTransformZ(*track_endZ);
      }
      else {
        xS = doCoordTransformX(*track_endX + x_offset);
        yS = doCoordTransformY(*track_endY);
        zS = doCoordTransformZ(*track_endZ);
    
        xE = doCoordTransformX(*track_startX + x_offset);
        yE = doCoordTransformY(*track_startY);
        zE = doCoordTransformZ(*track_startZ);
      }

      Zhist_1s->Fill(zS);
      Zhist_1e->Fill(zE);
      
      
      if (((xS < maxXdist) && (xE < maxXdist)) || ((xS > (Lx - maxXdist)) && (xE > (Lx - maxXdist))) || ((xS < maxXdist) && (yS < maxYdist)) || ((xS < maxXdist) && (yS > (Ly - maxYdist))) || ((xS < maxXdist) && (zS < maxZdist)) || ((xS < maxXdist) && (zS > (Lz -maxZdist))) || ((xE < maxXdist) && (yE < maxYdist)) || ((xE < maxXdist) && (yE > (Ly - maxYdist))) || ((xE < maxXdist) && (zE < maxZdist)) || ((xE < maxXdist) && (zE > (Lz -maxZdist))) || ((yS < maxYdist) && (zS < maxZdist)) || ((yS > (Ly - maxYdist)) && (zS < maxZdist)) || ((yS < maxYdist) && (zS > (Lz - maxZdist))) || ((yS > (Ly - maxYdist)) && (zS > (Lz - maxZdist))) || ((yE < maxYdist) && (zE < maxZdist)) || ((yE > (Ly - maxYdist)) && (zE < maxZdist)) || ((yE < maxYdist) && (zE > (Lz - maxZdist))) || ((yE > (Ly - maxYdist)) && (zE > (Lz - maxZdist)))) continue;
      
      

      Zhist_2s->Fill(zS);
      Zhist_2e->Fill(zE);
      
      if (((xS > maxXdist) && (xS < (Lx - maxXdist)) && (yS > maxYdist) && (yS < (Ly - maxYdist)) && (zS > maxZdist) && (zS < (Lz - maxZdist))) || ((xE > maxXdist) && (xE < (Lx - maxXdist)) && (yE > maxYdist) && (yE < (Ly - maxYdist)) && (zE > maxZdist) && (zE < (Lz - maxZdist)))) continue;
      
      

      Zhist_3s->Fill(zS);
      Zhist_3e->Fill(zE);
    
      if ((((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_anode)) continue;
      
      

      Zhist_4s->Fill(zS);
      Zhist_4e->Fill(zE);
      
      if ((((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_cathode)) continue;

      Zhist_5s->Fill(zS);
      Zhist_5e->Fill(zE);

      if ((((xS > (Lx - maxXdist)) && (xE < maxXdist)) || ((xE > (Lx - maxXdist)) && (xS < maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_crossing)) continue;

      Zhist_6s->Fill(zS);
      Zhist_6e->Fill(zE);

      double randNum = rand->Uniform(1.0);
      if (randNum > max(weightFunc->Eval(zS),weightFunc->Eval(zE))) continue;

      Zhist_7s->Fill(zS);
      Zhist_7e->Fill(zE);

      nTracks++;
      ++containedTracks;
      
      

      double SCEfactor = 1.0;
      if (isSCEon == false) {
        SCEfactor = 0.0;
      }
      
      // Correct track end point for cathode-piercing track (end furthest from cathode)
      if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
        if (xS < xE) {
          //std::cout << getTruthFwdOffset(0.0,yS+getTruthOffset(xS,yS,zS,2),zS+getTruthOffset(xS,yS,zS,3),1) << std::endl;
	  xE += SCEfactor*getTruthFwdOffset(0.0,yS+getTruthOffset(xS,yS,zS,2),zS+getTruthOffset(xS,yS,zS,3),1);
	}
	else {
          xS += SCEfactor*getTruthFwdOffset(0.0,yE+getTruthOffset(xE,yE,zE,2),zE+getTruthOffset(xE,yE,zE,3),1);
	  //std::cout << getTruthFwdOffset(0.0,yE+getTruthOffset(xE,yE,zE,2),zE+getTruthOffset(xE,yE,zE,3),1) << std::endl;
	}
      }
      
      if(xS < maxXdist) {
        //std::cout << getTruthOffset(xS,yS,zS,2) << std::endl;
	x0 = 0.0;
        y0 = yS + SCEfactor*getTruthOffset(xS,yS,zS,2); // TEMP FACTOR
        z0 = zS + SCEfactor*getTruthOffset(xS,yS,zS,3); // TEMP FACTOR
      }
      else if (xS > (Lx - maxXdist)) {
        x0 = Lx;
        y0 = yS + SCEfactor*getTruthOffset(xS,yS,zS,2);
        z0 = zS + SCEfactor*getTruthOffset(xS,yS,zS,3);
      }
      else if (min(fabs(yS),fabs(Ly-yS)) < min(fabs(zS),fabs(Lz-zS))) {
        if (fabs(yS) < fabs(Ly-yS)) {
          x0 = xS + SCEfactor*getTruthOffset(xS,yS,zS,1);
    	  y0 = 0.0;
    	  z0 = zS + SCEfactor*getTruthOffset(xS,yS,zS,3);
        }
        else {
          
	  x0 = xS + SCEfactor*getTruthOffset(xS,yS,zS,1);
    	  y0 = Ly;
    	  z0 = zS + SCEfactor*getTruthOffset(xS,yS,zS,3);
        }
      }
      else {
        if (fabs(zS) < fabs(Lz-zS)) {
          x0 = xS + SCEfactor*getTruthOffset(xS,yS,zS,1);
    	  y0 = yS + SCEfactor*getTruthOffset(xS,yS,zS,2);
    	  z0 = 0.0;
        }
        else {
          x0 = xS + SCEfactor*getTruthOffset(xS,yS,zS,1);
    	  y0 = yS + SCEfactor*getTruthOffset(xS,yS,zS,2);
    	  z0 = Lz;
        }
      }
    
      if(xE < maxXdist) {
        x1 = 0.0;
        y1 = yE + SCEfactor*getTruthOffset(xE,yE,zE,2); // TEMP FACTOR
        z1 = zE + SCEfactor*getTruthOffset(xE,yE,zE,3); // TEMP FACTOR
      }
      else if (xE > (Lx - maxXdist)) {
        x1 = Lx;
        y1 = yE + SCEfactor*getTruthOffset(xE,yE,zE,2);
        z1 = zE + SCEfactor*getTruthOffset(xE,yE,zE,3);
      }
      else if (min(fabs(yE),fabs(Ly-yE)) < min(fabs(zE),fabs(Lz-zE))) {
        if (fabs(yE) < fabs(Ly-yE)) {
          x1 = xE + SCEfactor*getTruthOffset(xE,yE,zE,1);
    	  y1 = 0.0;
    	  z1 = zE + SCEfactor*getTruthOffset(xE,yE,zE,3);
        }
        else {
          //std::cout << getTruthOffset(xE,yE,zE,1) << std::endl;
	  x1 = xE + SCEfactor*getTruthOffset(xE,yE,zE,1);
    	  y1 = Ly;
    	  z1 = zE + SCEfactor*getTruthOffset(xE,yE,zE,3);
        }
      }
      else {
        if (fabs(zE) < fabs(Lz-zE)) {
          x1 = xE + SCEfactor*getTruthOffset(xE,yE,zE,1);
    	  y1 = yE + SCEfactor*getTruthOffset(xE,yE,zE,2);
    	  z1 = 0.0;
        }
        else {
          x1 = xE + SCEfactor*getTruthOffset(xE,yE,zE,1);
    	  y1 = yE + SCEfactor*getTruthOffset(xE,yE,zE,2);
    	  z1 = Lz;
        }
      }
    
      Double_t trackLength = sqrt(pow(x0-x1,2.0)+pow(y0-y1,2.0)+pow(z0-z1,2.0));
      //std::cout << x0 << " " << y0 << " " << z0 << " : " << x1 << " " << y1 << " " << z1 << std::endl;
      //std::cout << "Track L: " << trackLength << std::endl;
      
      //Double_t theta = ArcCos((y1-y0)/trackLength);
      //Double_t phi = ArcSin((x0-x1)/(trackLength*Sin(theta)));
      Double_t theta = acos((y1-y0)/trackLength);
      Double_t phi = acos((z1-z0)/(trackLength*sin(theta)));
      if (x1 > x0) {
        phi = -1.0*fabs(phi);
      }
      else {
        phi = fabs(phi);
      }
      
      track.energy = *track_MCS/1000.0;
      track.pdgID = 13;
      track.x0 = x0;
      track.y0 = y0;
      track.z0 = z0;
      track.x1 = x1;
      track.y1 = y1;
      track.z1 = z1;
      track.theta = theta;
      track.phi = phi;
      track.electrons.clear();
    
      for(Int_t j = 0; j < *nPoints; j++)
      {
        //if ((doBulk == true) && ((j % 10) != 0) && (j != *nPoints-1)) continue;
        //if ((doBulk == false) && (j > 52) && (j < *nPoints-53) && (j != *nPoints-1)) continue;
        if (((j % 10) != 0) && (j != *nPoints-1)) continue;
	
        electron.x = -1.0; // Dummy (currently don't save this info in file)
        electron.y = -1.0; // Dummy (currently don't save this info in file)
        electron.z = -1.0; // Dummy (currently don't save this info in file)
        electron.t = -1.0; // Dummy (currently don't save this info in file)
        //electron.x_mod = doCoordTransformX(pointX[j]+x_offset);
        if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
	  if (x0 < x1) {
	    electron.x_mod = doCoordTransformX(pointX[j]+x_offset)+SCEfactor*getTruthFwdOffset(0.0,y0,z0,1);
	  }
	  else {
	    electron.x_mod = doCoordTransformX(pointX[j]+x_offset)+SCEfactor*getTruthFwdOffset(0.0,y1,z1,1);
	  }
        }
        else {
          electron.x_mod = doCoordTransformX(pointX[j]+x_offset);
        }
        electron.y_mod = doCoordTransformY(pointY[j]);
        electron.z_mod = doCoordTransformZ(pointZ[j]);
        electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
        electron.fate = 0; // Dummy (currently don't save this info in file)
    
        track.electrons.push_back(electron);
      }
      
      tracks.push_back(track);
    }

    outputFile->cd();
    
    Zhist_1s->Write();
    Zhist_2s->Write();
    Zhist_3s->Write();
    Zhist_4s->Write();
    Zhist_5s->Write();
    Zhist_6s->Write();
    Zhist_7s->Write();

    Zhist_1e->Write();
    Zhist_2e->Write();
    Zhist_3e->Write();
    Zhist_4e->Write();
    Zhist_5e->Write();
    Zhist_6e->Write();
    Zhist_7e->Write();
    
    std::cout << "N Tracks: " << nTracks << std::endl;
    std::cout << "Contained tracks: " << containedTracks << std::endl;
    std::cout << "Input Track Num: " << inputTrackNum << std::endl;

    return tracks;
  }
  else
  {
    return tracks;
  }
}

vector<trackInfo> getTrackSet(Int_t inputType)
{
  vector<trackInfo> tracks;

  TFile* mapFile;
  Char_t* mapName = (Char_t*)"";

  if(inputType == 1)
  {
    mapFile = new TFile(inputFileLaser,"READ");
    mapName = (Char_t*)"laser";
  }
  else if(inputType == 2)
  {
    mapFile = new TFile(inputFileCosmic,"READ");
    mapName = (Char_t*)"cosmic";
  }
  else
  {
    return tracks;
  }

  TTreeReader reader(Form("SpaCEtree_%s",mapName), mapFile);
  TTreeReaderValue<Double_t> x0(reader, Form("x0_%s",mapName));
  TTreeReaderValue<Double_t> y0(reader, Form("y0_%s",mapName));
  TTreeReaderValue<Double_t> z0(reader, Form("z0_%s",mapName));
  TTreeReaderValue<Double_t> theta(reader, Form("theta_%s",mapName));
  TTreeReaderValue<Double_t> phi(reader, Form("phi_%s",mapName));
  TTreeReaderValue<Int_t> nElec(reader, Form("nElec_%s",mapName));
  TTreeReaderArray<Double_t> elecX(reader, Form("elecX_%s",mapName));
  TTreeReaderArray<Double_t> elecY(reader, Form("elecY_%s",mapName));
  TTreeReaderArray<Double_t> elecZ(reader, Form("elecZ_%s",mapName));

  trackInfo track;
  elecInfo electron;
  
  while (reader.Next())
  {
    track.energy = -1.0;
    track.pdgID = -1;
    track.x0 = *x0;
    track.y0 = *y0;
    track.z0 = *z0;
    track.x1 = -1.0;
    track.y1 = -1.0;
    track.z1 = -1.0;
    track.theta = *theta;
    track.phi = *phi;
    track.electrons.clear();
    
    for(Int_t j = 0; j < *nElec; j++)
    {
      electron.x = -1.0; // Dummy (currently don't save this info in file)
      electron.y = -1.0; // Dummy (currently don't save this info in file)
      electron.z = -1.0; // Dummy (currently don't save this info in file)
      electron.t = -1.0; // Dummy (currently don't save this info in file)
      electron.x_mod = elecX[j];
      electron.y_mod = elecY[j];
      electron.z_mod = elecZ[j];
      electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
      electron.fate = 0; // Dummy (currently don't save this info in file)

      track.electrons.push_back(electron);
    }
    
    tracks.push_back(track);
  }

  return tracks;
}

vector<calibTrackInfo> makeCalibTracks(const vector<trackInfo> &tracks)
{
  vector<calibTrackInfo> calibTracks;
  calibTrackInfo calibTrack;
  trackInfo track;

  for(Int_t i = 0; i < tracks.size(); i++)
  {
    track = tracks.at(i);

    calibTrack.track = track;
    calibTrack.x0_calib = track.x0;
    calibTrack.y0_calib = track.y0;
    calibTrack.z0_calib = track.z0;
    calibTrack.x1_calib = track.x1;
    calibTrack.y1_calib = track.y1;
    calibTrack.z1_calib = track.z1;
    calibTrack.theta_calib = track.theta;
    calibTrack.phi_calib = track.phi;

    calibTrack.DxVec.clear();
    calibTrack.DyVec.clear();
    calibTrack.DzVec.clear();
    for(Int_t j = 0; j < track.electrons.size(); j++)
    {
      calibTrack.DxVec.push_back(0.0);
      calibTrack.DyVec.push_back(0.0);
      calibTrack.DzVec.push_back(0.0);
    }

    calibTrack.calibFlag = false;

    calibTracks.push_back(calibTrack);
  }

  return calibTracks;
}

void doLaserLaserCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numLaserTracks = laserCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 1;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  T_crossings->SetDirectory(outputFile);

  for(Int_t i = 0; i < numLaserTracks; i++)
  {
    cout << "LASER-LASER " << i << endl;

    calibTrackA = laserCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;

    for(Int_t j = i+1; j < numLaserTracks; j++)
    {
      calibTrackB = laserCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;

      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);

      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;
      //if((distVal < 0.0) || (distVal > 0.4)) continue;
      
      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);

      if(distValDistorted > 3.0*distVal+0.005) continue;
      //if(distValDistorted > 0.4) continue;

      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      
      
      trackNum1 = i;
      trackNum2 = j;
      T_crossings->Fill();

      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
  }

  return;
}

void doLaserCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numLaserTracks = laserCalibTracks.size();
  Int_t numCosmicTracks = cosmicCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 2;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  T_crossings->SetDirectory(outputFile);

  for(Int_t i = 0; i < numLaserTracks; i++)
  {
    cout << "LASER-COSMIC " << i << endl;

    calibTrackA = laserCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;

    for(Int_t j = 0; j < numCosmicTracks; j++)
    {
      calibTrackB = cosmicCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;
      if((cosmicTruthMode == 0) && (calibTrackB.calibFlag == false)) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;
      
      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);

      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;

      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);

      if(distValDistorted > 3.0*distVal+0.005) continue;

      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      

      trackNum1 = i;
      trackNum2 = j;
      T_crossings->Fill();

      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);

      if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
  }

  return;
}

void doCosmicCosmicCalib(const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numCosmicTracks = cosmicCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 3;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  T_crossings->SetDirectory(outputFile);

  for(Int_t i = 0; i < numCosmicTracks; i++)
  {
    cout << "COSMIC-COSMIC " << i << endl;

    calibTrackA = cosmicCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;
    if((cosmicTruthMode == 0) && (calibTrackA.calibFlag == false)) continue;

    for(Int_t j = i+1; j < numCosmicTracks; j++)
    {
      calibTrackB = cosmicCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;
      if((cosmicTruthMode == 0) && (calibTrackB.calibFlag == false)) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;
      
      //std::cout << i << " " << calibTrackA.track.x0 << " " << calibTrackA.track.y0 << " " << calibTrackA.track.z0 << " " << calibTrackA.track.theta << " " << calibTrackA.track.phi << std::endl;
      
      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);

      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;
      //if((distVal < 0.0) || (distVal > 0.1)) continue;
      
      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);

      if(distValDistorted > 3.0*distVal+0.005) continue;
      //if(distValDistorted > 0.1) continue;
      std::cout << "Track pair: " << i << ", " << j << std::endl;
      
      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      

      trackNum1 = i;
      trackNum2 = j;
      T_crossings->Fill();

      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);

     
      
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);
      
       if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
  }

  return;
}

void updateCalibTrueTrack(calibTrackInfo &calibTrack)
{
  TPrincipal pointcloud(3,"D");
  Double_t* point = new Double_t[3];
  Int_t numCalibPoints = 0;
  for(Int_t i = 0; i < calibTrack.track.electrons.size(); i++)
  {
    if(((calibTrack.DxVec.at(i) != 0.0) && (calibTrack.DyVec.at(i) != 0.0) && (calibTrack.DzVec.at(i) != 0.0)) || ((TMath::Ceil((calibTrack.track.electrons.at(i).x_mod/Lx)*nCalibDivisions_x) == nCalibDivisions_x))) // if corresponding cell calibrated or if at anode
    {
      point[0] = calibTrack.track.electrons.at(i).x_mod+calibTrack.DxVec.at(i);
      point[1] = calibTrack.track.electrons.at(i).y_mod+calibTrack.DyVec.at(i);
      point[2] = calibTrack.track.electrons.at(i).z_mod+calibTrack.DzVec.at(i);
      pointcloud.AddRow(point);

      numCalibPoints++;
    }
  }
  delete[] point;

  if(numCalibPoints >= 2)
  {
    pointcloud.MakePrincipals();
    
    const TMatrixD* evectors = pointcloud.GetEigenVectors();
    TVectorD unitVecs(3);
    unitVecs = TMatrixDColumn_const(*evectors,0);
    Double_t unitVecX = unitVecs(0);
    Double_t unitVecY = unitVecs(1);
    Double_t unitVecZ = unitVecs(2);
    
    const TVectorD* means = pointcloud.GetMeanValues();
    Double_t x_mean = (*means)(0);
    Double_t y_mean = (*means)(1);
    Double_t z_mean = (*means)(2);

    calibTrack.calibFlag = true;

    if(fabs(calibTrack.track.theta-ArcCos(unitVecY)) < fabs(calibTrack.track.theta-ArcCos(-1.0*unitVecY)))
      calibTrack.theta_calib = ArcCos(unitVecY);
    else
      calibTrack.theta_calib = ArcCos(-1.0*unitVecY);

    calibTrack.phi_calib = ArcTan(-1.0*unitVecX/unitVecZ);
    if(fabs(calibTrack.track.phi-(calibTrack.phi_calib+piVal)) < fabs(calibTrack.track.phi-calibTrack.phi_calib))
      calibTrack.phi_calib += piVal;
    if(fabs(calibTrack.track.phi-(calibTrack.phi_calib+piVal)) < fabs(calibTrack.track.phi-calibTrack.phi_calib))
      calibTrack.phi_calib += piVal;

    Double_t paramValA = -999999999999.0;
    Double_t paramValB = 999999999999.0;

    if((((0.0-x_mean)/unitVecX) > paramValA) && (((0.0-x_mean)/unitVecX) < 0.0))
      paramValA = (0.0-x_mean)/unitVecX;
    if((((0.0-x_mean)/unitVecX) < paramValB) && (((0.0-x_mean)/unitVecX) > 0.0))
      paramValB = (0.0-x_mean)/unitVecX;

    if((((Lx-x_mean)/unitVecX) > paramValA) && (((Lx-x_mean)/unitVecX) < 0.0))
      paramValA = (Lx-x_mean)/unitVecX;
    if((((Lx-x_mean)/unitVecX) < paramValB) && (((Lx-x_mean)/unitVecX) > 0.0))
      paramValB = (Lx-x_mean)/unitVecX;

    if((((0.0-y_mean)/unitVecY) > paramValA) && (((0.0-y_mean)/unitVecY) < 0.0))
      paramValA = (0.0-y_mean)/unitVecY;
    if((((0.0-y_mean)/unitVecY) < paramValB) && (((0.0-y_mean)/unitVecY) > 0.0))
      paramValB = (0.0-y_mean)/unitVecY;

    if((((Ly-y_mean)/unitVecY) > paramValA) && (((Ly-y_mean)/unitVecY) < 0.0))
      paramValA = (Ly-y_mean)/unitVecY;
    if((((Ly-y_mean)/unitVecY) < paramValB) && (((Ly-y_mean)/unitVecY) > 0.0))
      paramValB = (Ly-y_mean)/unitVecY;

    if((((0.0-z_mean)/unitVecZ) > paramValA) && (((0.0-z_mean)/unitVecZ) < 0.0))
      paramValA = (0.0-z_mean)/unitVecZ;
    if((((0.0-z_mean)/unitVecZ) < paramValB) && (((0.0-z_mean)/unitVecZ) > 0.0))
      paramValB = (0.0-z_mean)/unitVecZ;

    if((((Lz-z_mean)/unitVecZ) > paramValA) && (((Lz-z_mean)/unitVecZ) < 0.0))
      paramValA = (Lz-z_mean)/unitVecZ;
    if((((Lz-z_mean)/unitVecZ) < paramValB) && (((Lz-z_mean)/unitVecZ) > 0.0))
      paramValB = (Lz-z_mean)/unitVecZ;

    Double_t xValA = x_mean + paramValA*unitVecX;
    Double_t yValA = y_mean + paramValA*unitVecY;
    Double_t zValA = z_mean + paramValA*unitVecZ;

    Double_t xValB = x_mean + paramValB*unitVecX;
    Double_t yValB = y_mean + paramValB*unitVecY;
    Double_t zValB = z_mean + paramValB*unitVecZ;

    if(sqrt(pow(calibTrack.track.x0-xValA,2)+pow(calibTrack.track.y0-yValA,2)+pow(calibTrack.track.z0-zValA,2)) < sqrt(pow(calibTrack.track.x0-xValB,2)+pow(calibTrack.track.y0-yValB,2)+pow(calibTrack.track.z0-zValB,2)))
    {
      calibTrack.x0_calib = xValA;
      calibTrack.y0_calib = yValA;
      calibTrack.z0_calib = zValA;
      calibTrack.x1_calib = xValB;
      calibTrack.y1_calib = yValB;
      calibTrack.z1_calib = zValB;
    }
    else
    {
      calibTrack.x0_calib = xValB;
      calibTrack.y0_calib = yValB;
      calibTrack.z0_calib = zValB;
      calibTrack.x1_calib = xValA;
      calibTrack.y1_calib = yValA;
      calibTrack.z1_calib = zValA;
    }
  }

  return;
}

void updateAllCalibTrueTracks(vector<calibTrackInfo> &calibTracks, Int_t iterNum)
{
  Int_t numCalibTracks = calibTracks.size();

  TH1F *phiDiff = new TH1F(Form("phiDiff_iter%d",iterNum),Form("Calib. #phi - Actual #phi:  Iteration #%d",iterNum),40,-piVal/4,piVal/4);
  TH1F *thetaDiff = new TH1F(Form("thetaDiff_iter%d",iterNum),Form("Calib. #theta - Actual #theta:  Iteration #%d",iterNum),40,-piVal/8,piVal/8);
  TH1F *x0Diff = new TH1F(Form("x0Diff_iter%d",iterNum),Form("Calib. x_{0} - Actual x_{0}:  Iteration #%d",iterNum),40,-0.25,0.25);
  TH1F *y0Diff = new TH1F(Form("y0Diff_iter%d",iterNum),Form("Calib. y_{0} - Actual y_{0}:  Iteration #%d",iterNum),40,-0.25,0.25);
  TH1F *z0Diff = new TH1F(Form("z0Diff_iter%d",iterNum),Form("Calib. z_{0} - Actual z_{0}:  Iteration #%d",iterNum),40,-0.25,0.25);
  TH1F *x1Diff = new TH1F(Form("x1Diff_iter%d",iterNum),Form("Calib. x_{1} - Actual x_{1}:  Iteration #%d",iterNum),40,-0.25,0.25);
  TH1F *y1Diff = new TH1F(Form("y1Diff_iter%d",iterNum),Form("Calib. y_{1} - Actual y_{1}:  Iteration #%d",iterNum),40,-0.25,0.25);
  TH1F *z1Diff = new TH1F(Form("z1Diff_iter%d",iterNum),Form("Calib. z_{1} - Actual z_{1}:  Iteration #%d",iterNum),40,-0.25,0.25);

  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;

  Double_t corr_xlow_ylow_zlow[3] = {0.0,0.0,0.0};
  Double_t corr_xlow_ylow_zhigh[3] = {0.0,0.0,0.0};
  Double_t corr_xlow_yhigh_zlow[3] = {0.0,0.0,0.0};
  Double_t corr_xlow_yhigh_zhigh[3] = {0.0,0.0,0.0};
  Double_t corr_xhigh_ylow_zlow[3] = {0.0,0.0,0.0};
  Double_t corr_xhigh_ylow_zhigh[3] = {0.0,0.0,0.0};
  Double_t corr_xhigh_yhigh_zlow[3] = {0.0,0.0,0.0};
  Double_t corr_xhigh_yhigh_zhigh[3] = {0.0,0.0,0.0};

  Double_t sumWeight;
  Double_t averageDx;
  Double_t averageDy;
  Double_t averageDz;

  for(Int_t i = 0; i < numCalibTracks; i++)
  {
    cout << "UPDATE COSMIC " << i << endl;

    calibTracks.at(i).DxVec.clear();
    calibTracks.at(i).DyVec.clear();
    calibTracks.at(i).DzVec.clear();

    for(Int_t j = 0; j < calibTracks.at(i).track.electrons.size(); j++)
    {     
      xValDistorted = calibTracks.at(i).track.electrons.at(j).x_mod;
      yValDistorted = calibTracks.at(i).track.electrons.at(j).y_mod;
      zValDistorted = calibTracks.at(i).track.electrons.at(j).z_mod;
      
      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);
      
      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      sumWeight = calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]; 
      if(sumWeight == 0.0)
      {
        averageDx = 0.0;
        averageDy = 0.0;
        averageDz = 0.0;
      }
      else
      {
        averageDx = (calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex])/sumWeight;
        averageDy = (calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex])/sumWeight;
        averageDz = (calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]*calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]*calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]*calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]+calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]*calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex])/sumWeight;
      }

      if(calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex] == 0.0)
      {
        corr_xlow_ylow_zlow[0] = averageDx;
        corr_xlow_ylow_zlow[1] = averageDy;
        corr_xlow_ylow_zlow[2] = averageDz;
      }
      else
      {
        corr_xlow_ylow_zlow[0] = calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex];
        corr_xlow_ylow_zlow[1] = calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex];
        corr_xlow_ylow_zlow[2] = calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex];
      }

      if(calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex] == 0.0)
      {
        corr_xlow_ylow_zhigh[0] = averageDx;
        corr_xlow_ylow_zhigh[1] = averageDy;
        corr_xlow_ylow_zhigh[2] = averageDz;
      }
      else
      {
        corr_xlow_ylow_zhigh[0] = calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex];
        corr_xlow_ylow_zhigh[1] = calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex];
        corr_xlow_ylow_zhigh[2] = calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex];
      }

      if(calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex] == 0.0)
      {
        corr_xlow_yhigh_zlow[0] = averageDx;
        corr_xlow_yhigh_zlow[1] = averageDy;
        corr_xlow_yhigh_zlow[2] = averageDz;
      }
      else
      {
        corr_xlow_yhigh_zlow[0] = calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex];
        corr_xlow_yhigh_zlow[1] = calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex];
        corr_xlow_yhigh_zlow[2] = calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex];
      }

      if(calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex] == 0.0)
      {
        corr_xlow_yhigh_zhigh[0] = averageDx;
        corr_xlow_yhigh_zhigh[1] = averageDy;
        corr_xlow_yhigh_zhigh[2] = averageDz;
      }
      else
      {
        corr_xlow_yhigh_zhigh[0] = calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex];
        corr_xlow_yhigh_zhigh[1] = calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex];
        corr_xlow_yhigh_zhigh[2] = calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex];
      }

      if(calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex] == 0.0)
      {
        corr_xhigh_ylow_zlow[0] = averageDx;
        corr_xhigh_ylow_zlow[1] = averageDy;
        corr_xhigh_ylow_zlow[2] = averageDz;
      }
      else
      {
        corr_xhigh_ylow_zlow[0] = calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex];
        corr_xhigh_ylow_zlow[1] = calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex];
        corr_xhigh_ylow_zlow[2] = calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex];
      }

      if(calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex] == 0.0)
      {
        corr_xhigh_ylow_zhigh[0] = averageDx;
        corr_xhigh_ylow_zhigh[1] = averageDy;
        corr_xhigh_ylow_zhigh[2] = averageDz;
      }
      else
      {
        corr_xhigh_ylow_zhigh[0] = calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex];
        corr_xhigh_ylow_zhigh[1] = calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex];
        corr_xhigh_ylow_zhigh[2] = calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex];
      }

      if(calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex] == 0.0)
      {
        corr_xhigh_yhigh_zlow[0] = averageDx;
        corr_xhigh_yhigh_zlow[1] = averageDy;
        corr_xhigh_yhigh_zlow[2] = averageDz;
      }
      else
      {
        corr_xhigh_yhigh_zlow[0] = calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex];
        corr_xhigh_yhigh_zlow[1] = calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex];
        corr_xhigh_yhigh_zlow[2] = calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex];
      }

      if(calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex] == 0.0)
      {
        corr_xhigh_yhigh_zhigh[0] = averageDx;
        corr_xhigh_yhigh_zhigh[1] = averageDy;
        corr_xhigh_yhigh_zhigh[2] = averageDz;
      }
      else
      {
        corr_xhigh_yhigh_zhigh[0] = calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex];
        corr_xhigh_yhigh_zhigh[1] = calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex];
        corr_xhigh_yhigh_zhigh[2] = calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex]/calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex];
      }

      calibTracks.at(i).DxVec.push_back(((corr_xlow_ylow_zlow[0]*(1.0-xCalibFrac) + corr_xhigh_ylow_zlow[0]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zlow[0]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zlow[0]*xCalibFrac)*yCalibFrac)*(1.0-zCalibFrac) + ((corr_xlow_ylow_zhigh[0]*(1.0-xCalibFrac) + corr_xhigh_ylow_zhigh[0]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zhigh[0]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zhigh[0]*xCalibFrac)*yCalibFrac)*zCalibFrac);
      calibTracks.at(i).DyVec.push_back(((corr_xlow_ylow_zlow[1]*(1.0-xCalibFrac) + corr_xhigh_ylow_zlow[1]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zlow[1]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zlow[1]*xCalibFrac)*yCalibFrac)*(1.0-zCalibFrac) + ((corr_xlow_ylow_zhigh[1]*(1.0-xCalibFrac) + corr_xhigh_ylow_zhigh[1]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zhigh[1]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zhigh[1]*xCalibFrac)*yCalibFrac)*zCalibFrac);
      calibTracks.at(i).DzVec.push_back(((corr_xlow_ylow_zlow[2]*(1.0-xCalibFrac) + corr_xhigh_ylow_zlow[2]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zlow[2]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zlow[2]*xCalibFrac)*yCalibFrac)*(1.0-zCalibFrac) + ((corr_xlow_ylow_zhigh[2]*(1.0-xCalibFrac) + corr_xhigh_ylow_zhigh[2]*xCalibFrac)*(1.0-yCalibFrac) + (corr_xlow_yhigh_zhigh[2]*(1.0-xCalibFrac) + corr_xhigh_yhigh_zhigh[2]*xCalibFrac)*yCalibFrac)*zCalibFrac);
    }

    updateCalibTrueTrack(calibTracks.at(i));

    if(calibTracks.at(i).calibFlag == true)
    {
      phiDiff->Fill(calibTracks.at(i).phi_calib-calibTracks.at(i).track.phi);
      thetaDiff->Fill(calibTracks.at(i).theta_calib-calibTracks.at(i).track.theta);
      x0Diff->Fill(calibTracks.at(i).x0_calib-calibTracks.at(i).track.x0);
      y0Diff->Fill(calibTracks.at(i).y0_calib-calibTracks.at(i).track.y0);
      z0Diff->Fill(calibTracks.at(i).z0_calib-calibTracks.at(i).track.z0);
      x1Diff->Fill(calibTracks.at(i).x1_calib-calibTracks.at(i).track.x1);
      y1Diff->Fill(calibTracks.at(i).y1_calib-calibTracks.at(i).track.y1);
      z1Diff->Fill(calibTracks.at(i).z1_calib-calibTracks.at(i).track.z1);
    }
  }

  outputFile->cd();
  phiDiff->Write();
  thetaDiff->Write();
  x0Diff->Write();
  y0Diff->Write();
  z0Diff->Write();
  x1Diff->Write();
  y1Diff->Write();
  z1Diff->Write();

  delete phiDiff;
  delete thetaDiff;
  delete x0Diff;
  delete y0Diff;
  delete z0Diff;

  return;
}

void doCalibration(const vector<trackInfo> &laserTracks, const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t cosmicTruthMode, Int_t numIterations)
{
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if(x == nCalibDivisions_x)
          calibWeight[x][y][z] = 1.0;
	else
          calibWeight[x][y][z] = 0.0;

        calibDeltaX[x][y][z] = 0.0;
        calibDeltaY[x][y][z] = 0.0;
        calibDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  vector<calibTrackInfo> laserCalibTracks = makeCalibTracks(laserTracks);
  vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);

  Double_t x_true, y_true, z_true;
  Double_t x_reco, y_reco, z_reco;
  Double_t Dx, Dy, Dz;
  Int_t elecFate;

  TTree *T_calib = new TTree("SpaCEtree_calib","SpaCEtree_calib");
  T_calib->Branch("x_true",&x_true,"data_calib/D");
  T_calib->Branch("y_true",&y_true,"data_calib/D");
  T_calib->Branch("z_true",&z_true,"data_calib/D");
  T_calib->Branch("x_reco",&x_reco,"data_calib/D");
  T_calib->Branch("y_reco",&y_reco,"data_calib/D");
  T_calib->Branch("z_reco",&z_reco,"data_calib/D");
  T_calib->Branch("Dx",&Dx,"data_calib/D");
  T_calib->Branch("Dy",&Dy,"data_calib/D");
  T_calib->Branch("Dz",&Dz,"data_calib/D");
  T_calib->Branch("elecFate",&elecFate,"data_calib/I");
  T_calib->SetDirectory(outputFile);

  //doLaserLaserCalib(laserCalibTracks,distScale,maxDistFactor);
  if(cosmicTruthMode == 0)
  {
    for(Int_t i = 0; i < numIterations; i++)
    {
      updateAllCalibTrueTracks(cosmicCalibTracks,2*i+1);
      doLaserCosmicCalib(laserCalibTracks,cosmicCalibTracks,distScale,maxDistFactor,cosmicTruthMode);
      updateAllCalibTrueTracks(cosmicCalibTracks,2*i+2);
      doCosmicCosmicCalib(cosmicCalibTracks,distScale,maxDistFactor,cosmicTruthMode);
    }
  }
  else
  {
    //doLaserCosmicCalib(laserCalibTracks,cosmicCalibTracks,distScale,maxDistFactor,cosmicTruthMode);
    doCosmicCosmicCalib(cosmicCalibTracks,distScale,maxDistFactor,cosmicTruthMode);
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if(calibWeight[x][y][z] > 0.0)
	{
          calibDeltaX[x][y][z] /= calibWeight[x][y][z];
          calibDeltaY[x][y][z] /= calibWeight[x][y][z];
          calibDeltaZ[x][y][z] /= calibWeight[x][y][z];
	}
      }
    }
  }

  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;
      z_reco = -1.0*Lz/nCalibDivisions_z;
  
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        z_reco += Lz/nCalibDivisions_z;
  
        Dx = calibDeltaX[x][y][z];
        Dy = calibDeltaY[x][y][z];
        Dz = calibDeltaZ[x][y][z];

        x_true = x_reco+Dx;
        y_true = y_reco+Dy;
        z_true = z_reco+Dz;

        if(calibWeight[x][y][z] > 0.0)
          elecFate = 1;
	else
          elecFate = 0;
    
        T_calib->Fill();
      }
    }
  }

  return;
}

void saveTrackInfo(const vector<trackInfo> &tracks)
{
  Int_t trackID;

  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Double_t theta;
  Double_t phi;
  Double_t energy;

  Int_t nElec;
  Double_t elecX[10000];
  Double_t elecY[10000];
  Double_t elecZ[10000];

  TTree *T_tracks;
  T_tracks = new TTree("SpaCEtree_tracks","SpaCEtree_tracks");
  T_tracks->Branch("trackID_tracks",&trackID,"trackID_tracks/I");
  T_tracks->Branch("x0_tracks",&x0,"x0_tracks/D");  
  T_tracks->Branch("y0_tracks",&y0,"y0_tracks/D");
  T_tracks->Branch("z0_tracks",&z0,"z0_tracks/D");
  T_tracks->Branch("x1_tracks",&x1,"x1_tracks/D");
  T_tracks->Branch("y1_tracks",&y1,"y1_tracks/D");
  T_tracks->Branch("z1_tracks",&z1,"z1_tracks/D");
  T_tracks->Branch("theta_tracks",&theta,"theta_tracks/D");
  T_tracks->Branch("phi_tracks",&phi,"phi_tracks/D");
  T_tracks->Branch("energy_tracks",&energy,"energy_tracks/D");
  T_tracks->Branch("nElec_tracks",&nElec,"nElec_tracks/I");
  T_tracks->Branch("elecX_tracks",&elecX,"elecX_tracks[nElec_tracks]/D");
  T_tracks->Branch("elecY_tracks",&elecY,"elecY_tracks[nElec_tracks]/D");
  T_tracks->Branch("elecZ_tracks",&elecZ,"elecZ_tracks[nElec_tracks]/D");
  T_tracks->SetDirectory(outputFile);

  Int_t numTracks = tracks.size();
  for(Int_t j = 0; j < numTracks; j++)
  {
    trackID = j;
    
    x0 = tracks.at(j).x0;
    y0 = tracks.at(j).y0;
    z0 = tracks.at(j).z0;
    x1 = tracks.at(j).x1;
    y1 = tracks.at(j).y1;
    z1 = tracks.at(j).z1;
    theta = tracks.at(j).theta;
    phi = tracks.at(j).phi;
    energy = tracks.at(j).energy;
    
    //if() continue; // make check to see if track is of "high momentum" (which is?)

    nElec = 0;
    Int_t numElectrons = tracks.at(j).electrons.size();
    while(nElec < tracks.at(j).electrons.size())
    {        
      elecX[nElec] = tracks.at(j).electrons.at(nElec).x_mod;
      elecY[nElec] = tracks.at(j).electrons.at(nElec).y_mod;
      elecZ[nElec] = tracks.at(j).electrons.at(nElec).z_mod;
      nElec++;
    }

    T_tracks->Fill();
  }

  return;
}

void loadTruthMap()
{
  TFile* fileTruth = new TFile("data/dispOutput_MicroBooNE_E273.root");

  TTreeReader reader("SpaCEtree_bkwdDisp", fileTruth);
  TTreeReaderValue<Double_t> reco_x(reader, "x_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> reco_y(reader, "y_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> reco_z(reader, "z_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dx(reader, "Dx.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dy(reader, "Dy.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dz(reader, "Dz.data_bkwdDisp");
  TTreeReaderValue<Int_t> elecFate(reader, "elecFate.data_bkwdDisp");

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        trueDeltaX[x][y][z] = 0.0;
        trueDeltaY[x][y][z] = 0.0;
        trueDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  while (reader.Next())
  {
    if (*elecFate == 1) {
      trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dx;
      trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dy;
      trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dz;
    }
    else {
      trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
      trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
      trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y-1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-2][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-2][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-2][z];
          }
          else if (y == 1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][2][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][2][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][2][z];
          }
          else if (z == nCalibDivisions_z-1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][nCalibDivisions_z-2];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][nCalibDivisions_z-2];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][nCalibDivisions_z-2];
          }
          else if (z == 1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][2];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][2];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][2];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if ((y == nCalibDivisions_y) && (z == 0)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][1];
          }
          else if ((y == nCalibDivisions_y) && (z == nCalibDivisions_z)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
          }
          else if ((y == 0) && (z == 0)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][1];
          }
          else if ((y == 0) && (z == nCalibDivisions_z)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][nCalibDivisions_z-1];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][z];
          }
          else if (y == 0) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][z];
          }
          else if (z == nCalibDivisions_z) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][nCalibDivisions_z-1];
          }
          else if (z == 0) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][1];
          }
	}
      }
    }
  }
  

  

  return;
}

Double_t getTruthOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.0) {
    xVal = 0.0;
  }
  if (xVal > Lx) {
    xVal = Lx;
  }

  if (yVal < 0.0) {
    yVal = 0.0;
  }
  if (yVal > Ly) {
    yVal = Ly;
  }

  if (zVal < 0.0) {
    zVal = 0.0;
  }
  if (zVal > Lz) {
    zVal = Lz;
  }

  if (comp == 1) {
    offset = trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 2) {
    offset = trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 3) {
    offset = trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  
  if ((comp == 1) && (offset < 0.0)) {
    offset = 0.0;
  }

  return offset;
}

void loadTruthFwdMap()
{
  TFile* fileTruth = new TFile("data/dispOutput_MicroBooNE_E273.root");

  TTreeReader reader("SpaCEtree_fwdDisp", fileTruth);
  TTreeReaderValue<Double_t> true_x(reader, "x_true.data_fwdDisp");
  TTreeReaderValue<Double_t> true_y(reader, "y_true.data_fwdDisp");
  TTreeReaderValue<Double_t> true_z(reader, "z_true.data_fwdDisp");
  TTreeReaderValue<Double_t> Dx(reader, "Dx.data_fwdDisp");
  TTreeReaderValue<Double_t> Dy(reader, "Dy.data_fwdDisp");
  TTreeReaderValue<Double_t> Dz(reader, "Dz.data_fwdDisp");
  TTreeReaderValue<Int_t> elecFate(reader, "elecFate.data_fwdDisp");

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        trueFwdDeltaX[x][y][z] = 0.0;
        trueFwdDeltaY[x][y][z] = 0.0;
        trueFwdDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  while (reader.Next())
  {
    if (*elecFate == 1) {
      trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dx;
      trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dy;
      trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dz;
    }
    else {
      trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
      trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
      trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y-1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-2][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-2][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-2][z];
          }
          else if (y == 1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][2][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][2][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][2][z];
          }
          else if (z == nCalibDivisions_z-1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][nCalibDivisions_z-2];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][nCalibDivisions_z-2];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][nCalibDivisions_z-2];
          }
          else if (z == 1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][2];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][2];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][2];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if ((y == nCalibDivisions_y) && (z == 0)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][1];
          }
          else if ((y == nCalibDivisions_y) && (z == nCalibDivisions_z)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
          }
          else if ((y == 0) && (z == 0)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][1];
          }
          else if ((y == 0) && (z == nCalibDivisions_z)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][nCalibDivisions_z-1];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][z];
          }
          else if (y == 0) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][z];
          }
          else if (z == nCalibDivisions_z) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][nCalibDivisions_z-1];
          }
          else if (z == 0) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][1];
          }
	}
      }
    }
  }
 

  return;
}

Double_t getTruthFwdOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.0) {
    xVal = 0.0;
  }
  if (xVal > Lx) {
    xVal = Lx;
  }

  if (yVal < 0.0) {
    yVal = 0.0;
  }
  if (yVal > Ly) {
    yVal = Ly;
  }

  if (zVal < 0.0) {
    zVal = 0.0;
  }
  if (zVal > Lz) {
    zVal = Lz;
  }

  if (comp == 1) {
    offset = trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 2) {
    offset = trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 3) {
    offset = trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  
  if ((comp == 1) && (offset > 0.0)) {
    offset = 0.0;
  }

  return offset;
}

void doCalibFaces(const vector<trackInfo> &cosmicTracks, Int_t minTrackPoints, Int_t numTrackSegPoints)
{
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      calibTopDeltaX[x][z] = 0.0;
      calibTopDeltaY[x][z] = 0.0;
      calibTopDeltaZ[x][z] = 0.0;
      calibTopWeight[x][z] = 0.0;      

      calibBottomDeltaX[x][z] = 0.0;
      calibBottomDeltaY[x][z] = 0.0;
      calibBottomDeltaZ[x][z] = 0.0;
      calibBottomWeight[x][z] = 0.0;
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      calibUpstreamDeltaX[x][y] = 0.0;
      calibUpstreamDeltaY[x][y] = 0.0;
      calibUpstreamDeltaZ[x][y] = 0.0;
      calibUpstreamWeight[x][y] = 0.0;

      calibDownstreamDeltaX[x][y] = 0.0;
      calibDownstreamDeltaY[x][y] = 0.0;
      calibDownstreamDeltaZ[x][y] = 0.0;
      calibDownstreamWeight[x][y] = 0.0;
    }
  }

  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      calibCathodeDeltaX[y][z] = 0.0;
      calibCathodeDeltaY[y][z] = 0.0;
      calibCathodeDeltaZ[y][z] = 0.0;
      calibCathodeWeight[y][z] = 0.0;
    }
  }

  vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);

  calibTrackInfo calibTrack;
  Int_t numCosmicTracks = cosmicCalibTracks.size();

  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Int_t trackNum;
  Int_t numElec;
  
  for(Int_t i = 0; i < numCosmicTracks; i++)
  {
    cout << "COSMIC-FACE " << i << endl;

    calibTrack = cosmicCalibTracks.at(i);
    if(calibTrack.track.electrons.size() < minTrackPoints) continue;

    trackNum = i;
    numElec = calibTrack.track.electrons.size();

    Int_t numNearAnode = 0;
    if (calibTrack.track.electrons.at(0).x_mod > (Lx - maxXdist)) {
      numNearAnode++;
    }
    if (calibTrack.track.electrons.at(numElec-1).x_mod > (Lx - maxXdist)) {
      numNearAnode++;
    }
    if((numNearAnode == 0) || (numNearAnode == 2)) continue;
	
    PointCloud calibPoints;
    
    Double_t avgX1 = 0.0;
    Double_t avgX2 = 0.0;

    for(Int_t j = 0; j < 5; j++)
    {
      Double_t elecX = calibTrack.track.electrons.at(j).x_mod;
      avgX1 += elecX;
    }
    avgX1 /= 5.0;

    for(Int_t j = numElec-1; j > numElec-6; j--)
    {
      Double_t elecX = calibTrack.track.electrons.at(j).x_mod;
      avgX2 += elecX;
    }
    avgX2 /= 5.0;

    double SCEfactor = 1.0;
    if (isSCEon == false) {
      SCEfactor = 0.0;
    }
    
    Int_t numBadPoints = 0;
    if(avgX1 < avgX2) {
      xValDistorted = calibTrack.track.electrons.at(0).x_mod;
      yValDistorted = calibTrack.track.electrons.at(0).y_mod;
      zValDistorted = calibTrack.track.electrons.at(0).z_mod;

      for (int j = numElec-1; j > numElec-1-numTrackSegPoints; j--) {
        if ((calibTrack.track.electrons.at(j).x_mod <= 0.0) || (calibTrack.track.electrons.at(j).x_mod >= Lx) || (calibTrack.track.electrons.at(j).y_mod <= 0.0) || (calibTrack.track.electrons.at(j).y_mod >= Ly) || (calibTrack.track.electrons.at(j).z_mod <= 0.0) || (calibTrack.track.electrons.at(j).z_mod >= Lz)) {
          numBadPoints++;
	  continue;
	}

	Point tempPoint;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod;
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod;
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod;
        tempPoint.x = calibTrack.track.electrons.at(j).x_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1); // TEMP USE OF TRUTH INFO
        tempPoint.y = calibTrack.track.electrons.at(j).y_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2); // TEMP USE OF TRUTH INFO
        tempPoint.z = calibTrack.track.electrons.at(j).z_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3); // TEMP USE OF TRUTH INFO
 
        calibPoints.push_back(tempPoint);
      }
    }
    else {
      xValDistorted = calibTrack.track.electrons.at(numElec-1).x_mod;
      yValDistorted = calibTrack.track.electrons.at(numElec-1).y_mod;
      zValDistorted = calibTrack.track.electrons.at(numElec-1).z_mod;

      for (int j = 0; j < numTrackSegPoints; j++) {
	if ((calibTrack.track.electrons.at(j).x_mod <= 0.0) || (calibTrack.track.electrons.at(j).x_mod >= Lx) || (calibTrack.track.electrons.at(j).y_mod <= 0.0) || (calibTrack.track.electrons.at(j).y_mod >= Ly) || (calibTrack.track.electrons.at(j).z_mod <= 0.0) || (calibTrack.track.electrons.at(j).z_mod >= Lz)) {
          numBadPoints++;
	  continue;
	}
	
	Point tempPoint;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod;
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod;
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod;
        tempPoint.x = calibTrack.track.electrons.at(j).x_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1); // TEMP USE OF TRUTH INFO
        tempPoint.y = calibTrack.track.electrons.at(j).y_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2); // TEMP USE OF TRUTH INFO
        tempPoint.z = calibTrack.track.electrons.at(j).z_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3); // TEMP USE OF TRUTH INFO

        calibPoints.push_back(tempPoint);
      }
    }
    if (numBadPoints > 5) continue;
    
    PCAResults results = DoPCA(calibPoints);

    Double_t startX;
    Double_t startY;
    Double_t startZ;
    if (results.endPoints.first(0) > results.endPoints.second(0)) {
      startX = results.endPoints.first(0);
      startY = results.endPoints.first(1);
      startZ = results.endPoints.first(2);
    }
    else {
      startX = results.endPoints.second(0);
      startY = results.endPoints.second(1);
      startZ = results.endPoints.second(2);
    }

    Double_t unitX = results.eVecs[0](0);
    Double_t unitY = results.eVecs[0](1);
    Double_t unitZ = results.eVecs[0](2);
    if (unitX > 0.0) {
      unitX *= -1.0;
      unitY *= -1.0;
      unitZ *= -1.0;
    }

    Int_t whichFace = -1;
    
    Double_t t_Top = (Ly-startY)/unitY;
    Double_t x_Top = unitX*t_Top + startX;
    Double_t y_Top = Ly;
    Double_t z_Top = unitZ*t_Top + startZ;
    if((t_Top > 0.0) && (x_Top > 0.0) && (x_Top < Lx) && (z_Top > 0.0) && (z_Top < Lz)) {
      xVal = x_Top;
      yVal = y_Top;
      zVal = z_Top;
      whichFace = 0;
    }

    Double_t t_Bottom = (0.0-startY)/unitY;
    Double_t x_Bottom = unitX*t_Bottom + startX;
    Double_t y_Bottom = 0.0;
    Double_t z_Bottom = unitZ*t_Bottom + startZ;
    if((t_Bottom > 0.0) && (x_Bottom > 0.0) && (x_Bottom < Lx) && (z_Bottom > 0.0) && (z_Bottom < Lz)) {
      xVal = x_Bottom;
      yVal = y_Bottom;
      zVal = z_Bottom;
      whichFace = 1;
    }

    Double_t t_Upstream = (0.0-startZ)/unitZ;
    Double_t x_Upstream = unitX*t_Upstream + startX;
    Double_t y_Upstream = unitY*t_Upstream + startY;
    Double_t z_Upstream = 0.0;
    if((t_Upstream > 0.0) && (x_Upstream > 0.0) && (x_Upstream < Lx) && (y_Upstream > 0.0) && (y_Upstream < Ly)) {
      xVal = x_Upstream;
      yVal = y_Upstream;
      zVal = z_Upstream;
      whichFace = 2;
    }

    Double_t t_Downstream = (Lz-startZ)/unitZ;
    Double_t x_Downstream = unitX*t_Downstream + startX;
    Double_t y_Downstream = unitY*t_Downstream + startY;
    Double_t z_Downstream = Lz;
    if((t_Downstream > 0.0) && (x_Downstream > 0.0) && (x_Downstream < Lx) && (y_Downstream > 0.0) && (y_Downstream < Ly)) {
      xVal = x_Downstream;
      yVal = y_Downstream;
      zVal = z_Downstream;
      whichFace = 3;
    }

    Double_t t_Cathode = (0.0-startX)/unitX;
    Double_t x_Cathode = 0.0;
    Double_t y_Cathode = unitY*t_Cathode + startY;
    Double_t z_Cathode = unitZ*t_Cathode + startZ;
    if((t_Cathode > 0.0) && (y_Cathode > 0.0) && (y_Cathode < Ly) && (z_Cathode > 0.0) && (z_Cathode < Lz)) {
      xVal = x_Cathode;
      yVal = y_Cathode;
      zVal = z_Cathode;
      whichFace = 4;
    }
    
    xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
    xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
    yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
    yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
    zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
    zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);
    
    xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
    yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
    zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);
    
    if(xValDistorted < 0.0) {
      xCalibLowIndex = 0;
      xCalibHighIndex = 1;
      xCalibFrac = 0.0;
    }
    else if (xValDistorted > Lx) {
      xCalibLowIndex = nCalibDivisions_x - 1;
      xCalibHighIndex = nCalibDivisions_x;
      xCalibFrac = 1.0;
    }
    
    if(yValDistorted < 0.0) {
      yCalibLowIndex = 0;
      yCalibHighIndex = 1;
      yCalibFrac = 0.0;
    }
    else if (yValDistorted > Ly) {
      yCalibLowIndex = nCalibDivisions_y - 1;
      yCalibHighIndex = nCalibDivisions_y;
      yCalibFrac = 1.0;
    }
    
    if(zValDistorted < 0.0) {
      zCalibLowIndex = 0;
      zCalibHighIndex = 1;
      zCalibFrac = 0.0;
    }
    else if (zValDistorted > Lz) {
      zCalibLowIndex = nCalibDivisions_z - 1;
      zCalibHighIndex = nCalibDivisions_z;
      zCalibFrac = 1.0;
    }
    

    
    if (whichFace == 0) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      calibTopWeight[xCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibTopDeltaX[xCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibTopDeltaY[xCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibTopDeltaZ[xCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      calibTopWeight[xCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibTopDeltaX[xCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibTopDeltaY[xCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibTopDeltaZ[xCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibTopWeight[xCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibTopDeltaX[xCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibTopDeltaY[xCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibTopDeltaZ[xCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibTopWeight[xCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibTopDeltaX[xCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibTopDeltaY[xCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibTopDeltaZ[xCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
    else if (whichFace == 1) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      calibBottomWeight[xCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibBottomDeltaX[xCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibBottomDeltaY[xCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibBottomDeltaZ[xCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      calibBottomWeight[xCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibBottomDeltaX[xCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibBottomDeltaY[xCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibBottomDeltaZ[xCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibBottomWeight[xCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibBottomDeltaX[xCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibBottomDeltaY[xCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibBottomDeltaZ[xCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibBottomWeight[xCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibBottomDeltaX[xCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibBottomDeltaY[xCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibBottomDeltaZ[xCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
    else if (whichFace == 2) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      calibUpstreamWeight[xCalibLowIndex][yCalibLowIndex] += tempFactor;
      calibUpstreamDeltaX[xCalibLowIndex][yCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibUpstreamDeltaY[xCalibLowIndex][yCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibUpstreamDeltaZ[xCalibLowIndex][yCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      calibUpstreamWeight[xCalibLowIndex][yCalibHighIndex] += tempFactor;
      calibUpstreamDeltaX[xCalibLowIndex][yCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibUpstreamDeltaY[xCalibLowIndex][yCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibUpstreamDeltaZ[xCalibLowIndex][yCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibUpstreamWeight[xCalibHighIndex][yCalibLowIndex] += tempFactor;
      calibUpstreamDeltaX[xCalibHighIndex][yCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibUpstreamDeltaY[xCalibHighIndex][yCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibUpstreamDeltaZ[xCalibHighIndex][yCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibUpstreamWeight[xCalibHighIndex][yCalibHighIndex] += tempFactor;
      calibUpstreamDeltaX[xCalibHighIndex][yCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibUpstreamDeltaY[xCalibHighIndex][yCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibUpstreamDeltaZ[xCalibHighIndex][yCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
    else if (whichFace == 3) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      calibDownstreamWeight[xCalibLowIndex][yCalibLowIndex] += tempFactor;
      calibDownstreamDeltaX[xCalibLowIndex][yCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDownstreamDeltaY[xCalibLowIndex][yCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDownstreamDeltaZ[xCalibLowIndex][yCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      calibDownstreamWeight[xCalibLowIndex][yCalibHighIndex] += tempFactor;
      calibDownstreamDeltaX[xCalibLowIndex][yCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDownstreamDeltaY[xCalibLowIndex][yCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDownstreamDeltaZ[xCalibLowIndex][yCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibDownstreamWeight[xCalibHighIndex][yCalibLowIndex] += tempFactor;
      calibDownstreamDeltaX[xCalibHighIndex][yCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibDownstreamDeltaY[xCalibHighIndex][yCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibDownstreamDeltaZ[xCalibHighIndex][yCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibDownstreamWeight[xCalibHighIndex][yCalibHighIndex] += tempFactor;
      calibDownstreamDeltaX[xCalibHighIndex][yCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibDownstreamDeltaY[xCalibHighIndex][yCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibDownstreamDeltaZ[xCalibHighIndex][yCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }
    else if (whichFace == 4) {
      tempFactor = (1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibCathodeWeight[yCalibLowIndex][zCalibLowIndex] += tempFactor;
      calibCathodeDeltaX[yCalibLowIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibCathodeDeltaY[yCalibLowIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibCathodeDeltaZ[yCalibLowIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = (1.0-yCalibFrac)*zCalibFrac;
      calibCathodeWeight[yCalibLowIndex][zCalibHighIndex] += tempFactor;
      calibCathodeDeltaX[yCalibLowIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibCathodeDeltaY[yCalibLowIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibCathodeDeltaZ[yCalibLowIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = yCalibFrac*(1.0-zCalibFrac);
      calibCathodeWeight[yCalibHighIndex][zCalibLowIndex] += tempFactor;
      calibCathodeDeltaX[yCalibHighIndex][zCalibLowIndex] += tempFactor*(xVal-xValDistorted);
      calibCathodeDeltaY[yCalibHighIndex][zCalibLowIndex] += tempFactor*(yVal-yValDistorted);
      calibCathodeDeltaZ[yCalibHighIndex][zCalibLowIndex] += tempFactor*(zVal-zValDistorted);

      tempFactor = yCalibFrac*zCalibFrac;
      calibCathodeWeight[yCalibHighIndex][zCalibHighIndex] += tempFactor;
      calibCathodeDeltaX[yCalibHighIndex][zCalibHighIndex] += tempFactor*(xVal-xValDistorted);
      calibCathodeDeltaY[yCalibHighIndex][zCalibHighIndex] += tempFactor*(yVal-yValDistorted);
      calibCathodeDeltaZ[yCalibHighIndex][zCalibHighIndex] += tempFactor*(zVal-zValDistorted);
    }

    if (whichFace == 0) {
      cout << "  " << whichFace << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << (xVal - xValDistorted) << " " << (yVal - yValDistorted) << " " << (zVal - zValDistorted) << " " << tempFactor << endl;
    }
  }
  
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      //Looks like a bug if weight < 0.0 ...
      if(calibTopWeight[x][z] > 0.0) {
        calibTopDeltaX[x][z] /= calibTopWeight[x][z];
        calibTopDeltaY[x][z] /= calibTopWeight[x][z];
        calibTopDeltaZ[x][z] /= calibTopWeight[x][z];
      }

      if(calibBottomWeight[x][z] > 0.0) {
        calibBottomDeltaX[x][z] /= calibBottomWeight[x][z];
        calibBottomDeltaY[x][z] /= calibBottomWeight[x][z];
        calibBottomDeltaZ[x][z] /= calibBottomWeight[x][z];
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      if(calibUpstreamWeight[x][y] > 0.0) {
        calibUpstreamDeltaX[x][y] /= calibUpstreamWeight[x][y];
        calibUpstreamDeltaY[x][y] /= calibUpstreamWeight[x][y];
        calibUpstreamDeltaZ[x][y] /= calibUpstreamWeight[x][y];
      }

      if(calibDownstreamWeight[x][y] > 0.0) {
        calibDownstreamDeltaX[x][y] /= calibDownstreamWeight[x][y];
        calibDownstreamDeltaY[x][y] /= calibDownstreamWeight[x][y];
        calibDownstreamDeltaZ[x][y] /= calibDownstreamWeight[x][y];
      }
    }
  }

  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      if(calibCathodeWeight[y][z] > 0.0) {
        calibCathodeDeltaX[y][z] /= calibCathodeWeight[y][z];
        calibCathodeDeltaY[y][z] /= calibCathodeWeight[y][z];
        calibCathodeDeltaZ[y][z] /= calibCathodeWeight[y][z];
      }
    }
  }

  Double_t x_true, y_true, z_true;
  Double_t x_reco, y_reco, z_reco;
  Double_t Dx, Dy, Dz;
  Int_t elecFate;
  Int_t faceNum;

  TTree *T_calibFaces = new TTree("SpaCEtree_calibFaces","SpaCEtree_calibFaces");
  T_calibFaces->Branch("x_true",&x_true,"data_calibFaces/D");
  T_calibFaces->Branch("y_true",&y_true,"data_calibFaces/D");
  T_calibFaces->Branch("z_true",&z_true,"data_calibFaces/D");
  T_calibFaces->Branch("x_reco",&x_reco,"data_calibFaces/D");
  T_calibFaces->Branch("y_reco",&y_reco,"data_calibFaces/D");
  T_calibFaces->Branch("z_reco",&z_reco,"data_calibFaces/D");
  T_calibFaces->Branch("Dx",&Dx,"data_calibFaces/D");
  T_calibFaces->Branch("Dy",&Dy,"data_calibFaces/D");
  T_calibFaces->Branch("Dz",&Dz,"data_calibFaces/D");
  T_calibFaces->Branch("elecFate",&elecFate,"data_calibFaces/I");
  T_calibFaces->Branch("faceNum",&faceNum,"data_calibFaces/I");
  T_calibFaces->SetDirectory(outputFile);

  faceNum = 0;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;

      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibTopDeltaX[x][z];
        Dy = calibTopDeltaY[x][z];
        Dz = calibTopDeltaZ[x][z];
      }

      x_true = x_reco+Dx;
      y_true = Ly;
      y_reco = y_true-Dy;
      z_true = z_reco+Dz;
  
      if((calibTopWeight[x][z] > 0.0) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      T_calibFaces->Fill();
    }
  }

  faceNum = 1;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;

      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibBottomDeltaX[x][z];
        Dy = calibBottomDeltaY[x][z];
        Dz = calibBottomDeltaZ[x][z];
      }

      x_true = x_reco+Dx;
      y_true = 0.0;
      y_reco = y_true-Dy;
      z_true = z_reco+Dz;
  
      if((calibBottomWeight[x][z] > 0.0) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      T_calibFaces->Fill();
    }
  }

  faceNum = 2;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;

      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibUpstreamDeltaX[x][y];
        Dy = calibUpstreamDeltaY[x][y];
        Dz = calibUpstreamDeltaZ[x][y];
      }

      x_true = x_reco+Dx;
      y_true = y_reco+Dy;
      z_true = 0.0;
      z_reco = z_true-Dz;
  
      if((calibUpstreamWeight[x][y] > 0.0) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      T_calibFaces->Fill();
    }
  }

  faceNum = 3;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;

      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibDownstreamDeltaX[x][y];
        Dy = calibDownstreamDeltaY[x][y];
        Dz = calibDownstreamDeltaZ[x][y];
      }

      x_true = x_reco+Dx;
      y_true = y_reco+Dy;
      z_true = Lz;
      z_reco = z_true-Dz;
  
      if((calibDownstreamWeight[x][y] > 0.0) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      T_calibFaces->Fill();
    }
  }

  faceNum = 4;
  y_reco = -1.0*Ly/nCalibDivisions_y;
  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    y_reco += Ly/nCalibDivisions_y;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;

      Dx = calibCathodeDeltaX[y][z];
      Dy = calibCathodeDeltaY[y][z];
      Dz = calibCathodeDeltaZ[y][z];

      x_true = 0.0;
      x_reco = x_true-Dx;
      y_true = y_reco+Dy;
      z_true = z_reco+Dz;
  
      if(calibCathodeWeight[y][z] > 0.0)
        elecFate = 1;
      else
        elecFate = 0;
    
      T_calibFaces->Fill();
    }
  }
  
  return;
}

PCAResults DoPCA(const PointCloud &points) {

  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;

  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;

  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }

  if (nThreeDHits == 0) {
    PCAResults results;
    return results; // FAIL FROM NO INPUT POINTS
  }

  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);

  // Define elements of our covariance matrix
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;

  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);

      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  // Using Eigen package
  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  //if (std::fabs(weightSum) < std::numeric_limits<float>::epsilon())
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - The total weight of three dimensional hits is 0!" << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  //if (eigenMat.info() != Eigen::ComputationInfo::Success)
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - PCA decompose failure, number of three D hits = " << nThreeDHits << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);

    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }

  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));

  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));

  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;

  return results;
}
