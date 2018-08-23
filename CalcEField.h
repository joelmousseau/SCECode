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

class eFieldCalculator{

  public:
     
     eFieldCalculator() {}
     
     ~eFieldCalculator(){}
     
     std::vector<double> cosmicToLaser(std::vector<double> inputVec);

     std::vector<double> laserToCosmic(std::vector<double> inputVec);

     std::vector<double> voxToCosmic(std::vector<int> vox);

     std::vector<int> cosmicToVox(std::vector<double> point);
     
     void compareCalib();

  private:

     const int nCalibDivisions_x = 25;
     const int nCalibDivisions_y = 25;
     const int nCalibDivisions_z = 100;
     
     const double Lx = 2.5;
     const double Ly = 2.5;
     const double Lz = 10.0;
     
     const double TPC_X = 256.04;
     const double TPC_Y = 232.50;
     const double TPC_Z = 1036.8;




};

#endif /* CalcEField_h */
