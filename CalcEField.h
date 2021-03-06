//
//  CalcEField.hpp
//  
//
//  Created by Joel Allen Mousseau on 7/9/18.
//

#ifndef CalcEField_h
#define CalcEField_h

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
#include <TError.h>
#include "TROOT.h"
#include <TChain.h>

#include "TBranch.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "math.h"
#include "vector"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "DistortionClass.h"

double calcChi2(const double *vals);

struct Point {
    float x;
    float y;
    float z;
};

struct PCAResults {
    TVector3 centroid;
    std::pair<TVector3,TVector3> endPoints;
    float length;
    TVector3 eVals;
    std::vector<TVector3> eVecs;
};

typedef std::vector<Point> PointCloud;

PCAResults DoPCA(const PointCloud &points);

class eFieldCalculator{

  public:
     
    eFieldCalculator() {gErrorIgnoreLevel = kFatal; useRelErrorCut = false; driftSign = 1.0; laserDriftVScale = 1.0; cosmicDriftVScale = 1.0;}
     
     ~eFieldCalculator(){}
    
     double doCoordTransformX(const Double_t inputX);
    
     double doCoordTransformY(const Double_t inputY);
    
     double doCoordTransformZ(const Double_t inputZ);
    
    double doInvCoordTransformX(const Double_t inputX);
    
     double doInvCoordTransformY(const Double_t inputY);
    
     double doInvCoordTransformZ(const Double_t inputZ);
     
     bool goodLaser(double value, double error){return (value < 1000.0  && error > 0.0 && isRelErrorSmall(value, error));}
     
     bool goodCosmic(double value, double error){return (isRelErrorSmall(value, error) );}
     
     bool isRelErrorSmall(double value, double error){return (useRelErrorCut ? (fabs(error/value) < minRelErr) : true);}
     
     std::vector<double> cosmicToLaser(std::vector<double> inputVec);

     std::vector<double> laserToCosmic(std::vector<double> inputVec);

     std::vector<double> voxToCosmic(std::vector<int> vox);

     std::vector<int> cosmicToVox(std::vector<double> point);
     
     void compareCalib(bool isData = false);
     
     void compareCalibZXPlane(bool isData = false);
     
     void compareTruth(bool sigmaDiff = false);
    
     void makeFwdMapPlots();
     
     double calculateFWHM(TH1F *hist);
     
     std::vector<int> calculateFWHMBins(TH1F *hist);
     
     double getDotProduct(std::vector<double> vecOne, std::vector<double> vecTwo);
     
     double getVectorMagnitude(std::vector<double> vec);
     
     double getAngle(std::vector<double> vecOne, std::vector<double> vecTwo);
     
     void  drawPlanarPlot(TH2F hist, int planeNum, const char *label, const char *filename, axisType axis, double zMax = 0.0);
    
    void  drawPlanarPlot(TH2F *hist, int planeNum, const char *label, const char *filename, axisType axis, double zMax = 0.0);
    
    void  draw1DPlot(TH1F *histOne, TH1F *histTwo, int planeNum, const char *label, const char *filename, axisType axis, double zMax = 0.0);
     
    void combineMaps(bool isData, bool skipLaser, bool skipCosmic);
    
    void combineMaps(int lowX = 0, int highX = 300, int lowY = 0, int highY = 300, int lowZ = 21, int highZ = 87, bool isData = true);
    
    void  compareFaces(bool isData = false);
    
    void combineWeightedMaps();
    
    void MakeDistorionTree();
    
    void MakeDistortionHistograms(bool isFwd = false);
    
    float LinInterp(float x, float x1, float x2, float q00, float q01);
    
    float TrilinInterp(float x, float y, float z, float q000, float q001, float q010, float q011, float q100, float q101, float q110, float q111, float x1, float x2, float y1, float y2, float z1, float z2);
    
    
    void doFits();
    
    void compareDataMC(std::string inputData, std::string inputMC, std::string outputTag);

    std::vector<double> studyResults2(std::string inputMapFileName = "IterativeLaserMap.root", std::string plotName = "NoPlot.png");
    
    std::vector<double> Residual_afterTrackCorr(std::string inputMapFileName = "IterativeLaserMap.root", std::string plotName = "NoPlot.png");
    
    void setDriftVScale(double laser, double cosmic);
    
    bool isInLaserRegion(double x, double y, double z);
    
    void makeCSVMap(std::string inputFileName, std::string outputFileName);
    
    double cosmicDriftVScale;
    double laserDriftVScale;
    
  private:

     const int nCalibDivisions_x = 25;
     const int nCalibDivisions_y = 25;
     const int nCalibDivisions_z = 100;
     
     const double Lx = 2.5;
     const double Ly = 2.5;
     const double Lz = 10.0;
     
    
    /*const double Lx = 256.04;
    const double Ly = 232.50;
    const double Lz = 1036.8;
    */
     const double TPC_X = 256.04;
     const double TPC_Y = 232.50;
     const double TPC_Z = 1036.8;
     
     const double minRelErr = 100.0;
     
     const double maxDiff   = 0.5;
     
     const double xMin = -1.0*(TPC_X/Lx)*Lx/(2.0*((double)nCalibDivisions_x));
     const double xMax = (TPC_X/Lx)*(Lx+Lx/(2.0*((double)nCalibDivisions_x)));
     const double yMin = (TPC_Y/Ly)*(-1.0*(Ly/2.0)-1.0*Ly/(2.0*((double)nCalibDivisions_y)));
     const double yMax = (TPC_Y/Ly)*((Ly/2.0)+Ly/(2.0*((double)nCalibDivisions_y)));
     const double zMin = -1.0*(TPC_Z/Lz)*Lz/(2.0*((double)nCalibDivisions_z));
     const double zMax = (TPC_Z/Lz)*(Lz+Lz/(2.0*((double)nCalibDivisions_z)));
    
    
     
     bool useRelErrorCut;
    
    double driftSign;

     TH3F *laser_dX;
     TH3F *laser_dY;
     TH3F *laser_dZ;
     
     TH3F *cosmic_dX;
     TH3F *cosmic_dY;
     TH3F *cosmic_dZ;
     
     
     TH3F *laser_dX_err;
     TH3F *laser_dY_err;
     TH3F *laser_dZ_err;
     
     TH3F *cosmic_dX_err;
     TH3F *cosmic_dY_err;
     TH3F *cosmic_dZ_err;
     
     
     
     
     
     


};

#endif /* CalcEField_h */
