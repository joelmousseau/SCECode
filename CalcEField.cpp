//
//  CalcEField.cpp
//  
//
//  Created by Joel Allen Mousseau on 7/9/18.
//

/*
 KEY: TH3F    Reco_Displacement_X;1    Reco Displacement X
 KEY: TH3F    Reco_Displacement_Y;1    Reco Displacement Y
 KEY: TH3F    Reco_Displacement_Z;1    Reco Displacement Z
 KEY: TH3F    Reco_Displacement_X_Error;1    Reco Deviation of Displacement X
 KEY: TH3F    Reco_Displacement_Y_Error;1    Reco Deviation of Displacement X
 KEY: TH3F    Reco_Displacement_Z_Error;1    Reco Deviation of Displacement X
 */
#include "CalcEField.h"
#inlcude "DistortionClass.h"

int main(int argc, char *argv[]){
    std::string inputCosmic = "BulkOnly.root";
    std::string inputLaser  = "RecoCorr-N3-S50-newselectbroadmuonSim-2side-Anode.root";
    TFile *cosmicFile =  new TFile(inputCosmic.c_str(),"READ");
    TFile *laserFile  =  new TFile(inputLaser.c_str(),"READ");
    TTreeReader reader("SCEtreeBulk_calib", cosmicFile);
    
    TH3F *h_xDist    = (TH3F*)laserFile->GetObject("Reco_Displacement_X");
    TH3F *h_xDistErr = (TH3F*)laserFile->GetObject("Reco_Displacement_X_Error");
    TTreeReaderValue<double> Dx(reader, "Dx");
    TTreeReaderValue<double> DxErr(reader, "DxErr");
    TTreeReaderValue<double> x_reco(reader, "x_reco");
    TTreeReaderValue<double> y_reco(reader, "y_reco");
    TTreeReaderValue<double> z_reco(reader, "z_reco");
    
    //Histograms of Difference
    const int nXBins = 25;
    TH2F *h_xDiff[nXBins];
    for (int i =0; i < nXBins; ++int) {
        h_xDiff = new TH2F(Form("xDiff%d", i), "X Difference", 100, 0, 10.0, 25, 0, 2.5 );
    }
    
    while(reader.Next()){
        std::vector<double> cosmicVox;
        cosmicVox.push_back(x_reco);
        cosmicVox.push_back(y_reco);
        cosmicVox.push_back(z_reco);
        
        std::vector<double> laserVox = cosmicToLaser(cosmicVox);
        
        int xBin = h_xDist->GetXaxis()->FindBin(laserVox[0]);
        int yBin = h_xDist->GetXaxis()->FindBin(laserVox[1]);
        int zBin = h_xDist->GetXaxis()->FindBin(laserVox[2]);
        
        double laserVal = h_xDist->GetBinContent(xBin, yBin, zBin);
        double Difference = (laserVal - 100.0*Dx);
        std::cout << Difference << std::endl;
        h_xDiff[xBin]->Fill(z_reco, y_reco, Difference);
        
    }
    
}

std::vector<double> cosmicToLaser(std::vector<double> inputVec){
    std::vector<double> return_vec;
    const double cathodeX = 250.0;
    const double bottomY  = 125.0;
    const double metersTocm = 100.0;
    //x transform
    return_vec.push_back((cathodeX - metersTocm*inputVec[0]));
    //y transform
    return_vec.push_back((metersTocm*inputVec[1] - bottomY));
    //z transform
    return_vec.push_back(metersTocm*inputVec[2]);
    return return_vec;
}


std::vector<double> laserToCosmic(std::vector<double> inputVec){
    std::vector<double> return_vec;
    const double cathodeX = 2.5;
    const double bottomY  = 1.25;
    const double cmToMeters = 0.01;
    //x transform
    return_vec.push_back(cathodeX - cmToMeters*inputVec[0]);
    //y transform
    return_vec.push_back(cmToMeters*inputVec[1] + bottomY);
    //z transform
    return_vec.push_back(cmToMeters*inputVec[2]);
    return return_vec;
}
