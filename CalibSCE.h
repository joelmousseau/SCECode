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



class elecInfo{
  
  public:
    double x;
    double y;
    double z;
    double t;
    double x_mod;
    double y_mod;
    double z_mod;
    double t_mod;
    int    fate;

};

class trackInfo{
  
  public:
    int pdgID;
    double energy;
    double x0;
    double y0;
    double z0;
    double x1;
    double y1;
    double z1;
    double theta;
    double phi;
    std::vector<elecInfo> electrons;

};

class calibTrackInfo{
  
  public:
    double x0_calib;
    double y0_calib;
    double z0_calib;
    double x1_calib;
    double y1_calib;
    double z1_calib;
    double theta_calib;
    double phi_calib;
    std::vector<double> DxVec;
    std::vector<double> DyVec;
    std::vector<double> DzVec;
    bool calibFlag;
    trackInfo track;
};

class  distortionVoxel{
   
   public:
     distortionVoxel();
     distortionVoxel(int x, int y, int z);
     ~distortionVoxel();
     std::vector<int> getVoxelNumber();
     int setVoxelNumber(std::vector<int> index);
     int addDistortion(float dist);
     std::vector<float> getDistortions();
     std::vector<float> getWeights();
     unsigned int getNPairs();
   
   private:
     std::vector<float> distortions;
     std::vector<float> weights;
     int xVoxel;
     int yVoxel;
     int zVoxel;  
               
};


class  distortionMap
{
   public:
     distortionMap();
     ~distortionMap();
     std::vector<  std::vector < std::vector<distortionVoxel> > > getMap(); 
     std::vecor< std::vector < std::vector<float> > > calculateMap()
     /*std::vecor< std::vector < std::vector<float> > > calculateMap();
     std::vecor< std::vector < std::vector<float> > > calculateMapErrors();*/
     float calculateMedian(distortionVoxel vox);
     /*float calculateMode()*/;
     float calculateMean(distortionVoxel vox);
     
     
     private:
       std::vector<  std::vector < std::vector<distortionVoxel> > > theMap;
       int nXVoxels = 26;
       int nYVoxels = 26;
       int nZVoxels = 101;
           
     

};

class SCECalib{
   

   public:
     calibSCE() { Intialize(); };
     ~calibSCE(){outputFile->Write(); outputFile->Close();}
     
     double doCoordTransformX(const double inputX);

     double doCoordTransformY(const double inputY);

     double doCoordTransformZ(const double inputZ);

     std::vector<double> getParabolaParameters(const std::vector<elecInfo> &parabola_points_track);

     std::vector<double> findClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);

     std::vector<double> findDistortedClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);

     std::vector<trackInfo> getLArSoftTrackSet(Int_t inputType);

     std::vector<trackInfo> getTrackSet(Int_t mapType);

     std::vector<calibTrackInfo> makeCalibTracks(const std::vector<trackInfo> &tracks);

     void doLaserLaserCalib(const std::vector<calibTrackInfo> &laserCalibTracks, double distScale, double maxDistFactor);

     void doLaserCosmicCalib(const std::vector<calibTrackInfo> &laserCalibTracks, const std::vector<calibTrackInfo> &cosmicCalibTracks, double distScale, double maxDistFactor, Int_t cosmicTruthMode);

     void doCosmicCosmicCalib(const std::vector<calibTrackInfo> &laserCalibTracks, double distScale, double maxDistFactor, Int_t cosmicTruthMode);

     void updateCalibTrueTrack(calibTrackInfo &calibTrack);

     void updateAllCalibTrueTracks(std::vector<calibTrackInfo> &calibTracks, Int_t iterNum);

     void doCalibration(const std::vector<trackInfo> &laserTracks, const std::vector<trackInfo> &cosmicTracks, double distScale, double maxDistFactor, Int_t cosmicTruthMode, Int_t numIterations);

     void saveTrackInfo(const std::vector<trackInfo> &tracks);

     void loadTruthMap(distortionMap &xMap, distortionMap &yMap, distortionMap &zMap);

     double getTruthOffset(double xVal, double yVal, double zVal, int comp);
  
  private:
     void Intialize();          
     std::string inputFileLaser;
     std::string inputFileCosmic;
     

     TFile *outputFile;

     double Lx;
     double Ly;
     double Lz;


     const bool isMC;
     const bool doBulk;

     double relAngleCut;
     double maxXdist;
     double maxYdist;
     double maxZdist;

     int minInputTrackNum;
     int maxInputTrackNum;

     //Use 30k tracks
     int maxCosmicTracks;


     double minTrackMCS_anode;
     double minTrackMCS_cathode;
     double minTrackMCS_crossing;

     int nCalibDivisions;

     double piVal;

     int nCalibDivisions_x;
     int nCalibDivisions_y;
     int nCalibDivisions_z; 
     
};
