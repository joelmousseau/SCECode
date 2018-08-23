//
//  CalcEField.cpp
//  
//
//  Created by Joel Allen Mousseau on 7/9/18.
//

/*
 Voxel is bin center
 KEY: TH3F    Reco_Displacement_X;1    Reco Displacement X
 KEY: TH3F    Reco_Displacement_Y;1    Reco Displacement Y
 KEY: TH3F    Reco_Displacement_Z;1    Reco Displacement Z
 KEY: TH3F    Reco_Displacement_X_Error;1    Reco Deviation of Displacement X
 KEY: TH3F    Reco_Displacement_Y_Error;1    Reco Deviation of Displacement X
 KEY: TH3F    Reco_Displacement_Z_Error;1    Reco Deviation of Displacement X
 */
#include "CalcEField.h"
#include "DistortionClass.h"

void eFieldCalculator::compareCalib()
{
  const double minRelErr = 100.0;
  
  const double xScale = TPC_X/Lx;
  const double yScale = TPC_Y/Ly;
  const double zScale = TPC_Z/Lz;

  double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TFile* fileLaser = new TFile("RecoCorr-N3-S50-LaserMC-2side-Anode.root");
  TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
  TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
  TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
  TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
  TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
  TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
  
  TFile* fileCosmic = new TFile("output_hists_MC_200k_Aug3.root");
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
  
  const double xMin = -1.0*(TPC_X/Lx)*Lx/(2.0*((double)nCalibDivisions_x));
  const double xMax = (TPC_X/Lx)*(Lx+Lx/(2.0*((double)nCalibDivisions_x)));
  const double yMin = (TPC_Y/Ly)*(-1.0*(Ly/2.0)-1.0*Ly/(2.0*((double)nCalibDivisions_y)));
  const double yMax = (TPC_Y/Ly)*((Ly/2.0)+Ly/(2.0*((double)nCalibDivisions_y)));
  const double zMin = -1.0*(TPC_Z/Lz)*Lz/(2.0*((double)nCalibDivisions_z));
  const double zMax = (TPC_Z/Lz)*(Lz+Lz/(2.0*((double)nCalibDivisions_z)));
  
  std::cout << "XMin: " << xMin << " YMin: " << yMin << std::endl;

  //Convert from cosmic to laser coordiantes
  TH3F* diff_dX = new TH3F("diff_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dY = new TH3F("diff_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dZ = new TH3F("diff_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);

  for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
  {
    for(Int_t j = 1; j <= diff_dX->GetNbinsY(); j++)
    {
      for(Int_t k = 1; k <= diff_dX->GetNbinsZ(); k++)
      {
	if((laser_dX_err->GetBinContent(i,j,k) > 0.0) && (fabs(laser_dX_err->GetBinContent(i,j,k)/laser_dX->GetBinContent(i,j,k)) < minRelErr) && (cosmic_dX_err->GetBinContent(i,j,k) > 0.0) && (fabs(cosmic_dX_err->GetBinContent(i,j,k)/cosmic_dX->GetBinContent(i,j,k)) < minRelErr))
	{
          diff_dX->SetBinContent(i,j,k,laser_dX->GetBinContent(i,j,k)-cosmic_dX->GetBinContent(i,j,k));
	}

	if((laser_dY_err->GetBinContent(i,j,k) > 0.0) && (fabs(laser_dY_err->GetBinContent(i,j,k)/laser_dY->GetBinContent(i,j,k)) < minRelErr) && (cosmic_dY_err->GetBinContent(i,j,k) > 0.0) && (fabs(cosmic_dY_err->GetBinContent(i,j,k)/cosmic_dY->GetBinContent(i,j,k)) < minRelErr))
	{
          diff_dY->SetBinContent(i,j,k,laser_dY->GetBinContent(i,j,k)-cosmic_dY->GetBinContent(i,j,k));
	}

	if((laser_dZ_err->GetBinContent(i,j,k) > 0.0) && (fabs(laser_dZ_err->GetBinContent(i,j,k)/laser_dZ->GetBinContent(i,j,k)) < minRelErr) && (cosmic_dZ_err->GetBinContent(i,j,k) > 0.0) && (fabs(cosmic_dZ_err->GetBinContent(i,j,k)/cosmic_dZ->GetBinContent(i,j,k)) < minRelErr))
	{
          diff_dZ->SetBinContent(i,j,k,laser_dZ->GetBinContent(i,j,k)-cosmic_dZ->GetBinContent(i,j,k));
	}
      }
    }
  }

  gStyle->SetTitleW(0.9);

  for(Int_t k = 1; k <= diff_dX->GetNbinsZ(); k++)
  {
    TH2F laser_2D_dX(Form("laser_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F laser_2D_dY(Form("laser_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F laser_2D_dZ(Form("laser_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);

    TH2F cosmic_2D_dX(Form("cosmic_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F cosmic_2D_dY(Form("cosmic_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F cosmic_2D_dZ(Form("cosmic_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);

    TH2F diff_2D_dX(Form("diff_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F diff_2D_dY(Form("diff_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F diff_2D_dZ(Form("diff_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);

    for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
    {
      for(Int_t j = 1; j <= diff_dX->GetNbinsY(); j++)
      {
        laser_2D_dX.SetBinContent(i,j,laser_dX->GetBinContent(i,j,k));
        laser_2D_dY.SetBinContent(i,j,laser_dY->GetBinContent(i,j,k));
        laser_2D_dZ.SetBinContent(i,j,laser_dZ->GetBinContent(i,j,k));

        cosmic_2D_dX.SetBinContent(i,j,cosmic_dX->GetBinContent(i,j,k));
        cosmic_2D_dY.SetBinContent(i,j,cosmic_dY->GetBinContent(i,j,k));
        cosmic_2D_dZ.SetBinContent(i,j,cosmic_dZ->GetBinContent(i,j,k));
	
        diff_2D_dX.SetBinContent(i,j,diff_dX->GetBinContent(i,j,k));
        diff_2D_dY.SetBinContent(i,j,diff_dY->GetBinContent(i,j,k));
        diff_2D_dZ.SetBinContent(i,j,diff_dZ->GetBinContent(i,j,k));
      }
    }

    TCanvas c_laser_dX(Form("canv_2D_dX_%d",k),"",600,600);
    c_laser_dX.cd();
    laser_2D_dX.SetTitle(Form("Laser #DeltaX [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    laser_2D_dX.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dX.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dX.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dX.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dX.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dX.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dX.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dX.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dX.SetStats(0);
    laser_2D_dX.Draw("COLZ");
    c_laser_dX.SaveAs(Form("plotdump/laser_2D_dX_%d.png",k));

    TCanvas c_laser_dY(Form("canv_2D_dY_%d",k),"",600,600);
    c_laser_dY.cd();
    laser_2D_dY.SetTitle(Form("Laser #DeltaY [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    laser_2D_dY.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dY.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dY.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dY.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dY.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dY.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dY.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dY.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dY.SetStats(0);
    laser_2D_dY.Draw("COLZ");
    c_laser_dY.SaveAs(Form("plotdump/laser_2D_dY_%d.png",k));
    
    TCanvas c_laser_dZ(Form("canv_2D_dZ_%d",k),"",600,600);
    c_laser_dZ.cd();
    laser_2D_dZ.SetTitle(Form("Laser #DeltaZ [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    laser_2D_dZ.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dZ.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dZ.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dZ.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dZ.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dZ.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dZ.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dZ.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dZ.SetStats(0);
    laser_2D_dZ.Draw("COLZ");
    c_laser_dZ.SaveAs(Form("plotdump/laser_2D_dZ_%d.png",k));

    TCanvas c_cosmic_dX(Form("canv_2D_dX_%d",k),"",600,600);
    c_cosmic_dX.cd();
    cosmic_2D_dX.SetTitle(Form("Cosmic #DeltaX [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    cosmic_2D_dX.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dX.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dX.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dX.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dX.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dX.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dX.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dX.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dX.SetStats(0);
    cosmic_2D_dX.Draw("COLZ");
    c_cosmic_dX.SaveAs(Form("plotdump/cosmic_2D_dX_%d.png",k));

    TCanvas c_cosmic_dY(Form("canv_2D_dY_%d",k),"",600,600);
    c_cosmic_dY.cd();
    cosmic_2D_dY.SetTitle(Form("Cosmic #DeltaY [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    cosmic_2D_dY.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dY.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dY.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dY.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dY.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dY.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dY.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dY.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dY.SetStats(0);
    cosmic_2D_dY.Draw("COLZ");
    c_cosmic_dY.SaveAs(Form("plotdump/cosmic_2D_dY_%d.png",k));
    
    TCanvas c_cosmic_dZ(Form("canv_2D_dZ_%d",k),"",600,600);
    c_cosmic_dZ.cd();
    cosmic_2D_dZ.SetTitle(Form("Cosmic #DeltaZ [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    cosmic_2D_dZ.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dZ.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dZ.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dZ.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dZ.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dZ.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dZ.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dZ.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dZ.SetStats(0);
    cosmic_2D_dZ.Draw("COLZ");
    c_cosmic_dZ.SaveAs(Form("plotdump/cosmic_2D_dZ_%d.png",k));

    TCanvas c_diff_dX(Form("canv_2D_dX_%d",k),"",600,600);
    c_diff_dX.cd();
    diff_2D_dX.SetTitle(Form("Laser-Cosmic #DeltaX [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    diff_2D_dX.GetXaxis()->SetTitle("X_{reco} [cm]");
    diff_2D_dX.GetXaxis()->SetTitleOffset(1.0);
    diff_2D_dX.GetXaxis()->SetTitleSize(0.04);
    diff_2D_dX.GetYaxis()->SetTitle("Y_{reco} [cm]");
    diff_2D_dX.GetYaxis()->SetTitleOffset(1.2);
    diff_2D_dX.GetYaxis()->SetTitleSize(0.04);
    diff_2D_dX.GetZaxis()->SetNoExponent(kTRUE);
    diff_2D_dX.GetZaxis()->SetLabelSize(0.025);
    diff_2D_dX.GetZaxis()->SetRangeUser(-3.0,3.0);
    diff_2D_dX.SetStats(0);
    diff_2D_dX.Draw("COLZ");
    c_diff_dX.SaveAs(Form("plotdump/diff_2D_dX_%d.png",k));

    TCanvas c_diff_dY(Form("canv_2D_dY_%d",k),"",600,600);
    c_diff_dY.cd();
    diff_2D_dY.SetTitle(Form("Laser-Cosmic #DeltaY [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    diff_2D_dY.GetXaxis()->SetTitle("X_{reco} [cm]");
    diff_2D_dY.GetXaxis()->SetTitleOffset(1.0);
    diff_2D_dY.GetXaxis()->SetTitleSize(0.04);
    diff_2D_dY.GetYaxis()->SetTitle("Y_{reco} [cm]");
    diff_2D_dY.GetYaxis()->SetTitleOffset(1.2);
    diff_2D_dY.GetYaxis()->SetTitleSize(0.04);
    diff_2D_dY.GetZaxis()->SetNoExponent(kTRUE);
    diff_2D_dY.GetZaxis()->SetLabelSize(0.025);
    diff_2D_dY.GetZaxis()->SetRangeUser(-3.0,3.0);
    diff_2D_dY.SetStats(0);
    diff_2D_dY.Draw("COLZ");
    c_diff_dY.SaveAs(Form("plotdump/diff_2D_dY_%d.png",k));
    
    TCanvas c_diff_dZ(Form("canv_2D_dZ_%d",k),"",600,600);
    c_diff_dZ.cd();
    diff_2D_dZ.SetTitle(Form("Laser-Cosmic #DeltaZ [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    diff_2D_dZ.GetXaxis()->SetTitle("X_{reco} [cm]");
    diff_2D_dZ.GetXaxis()->SetTitleOffset(1.0);
    diff_2D_dZ.GetXaxis()->SetTitleSize(0.04);
    diff_2D_dZ.GetYaxis()->SetTitle("Y_{reco} [cm]");
    diff_2D_dZ.GetYaxis()->SetTitleOffset(1.2);
    diff_2D_dZ.GetYaxis()->SetTitleSize(0.04);
    diff_2D_dZ.GetZaxis()->SetNoExponent(kTRUE);
    diff_2D_dZ.GetZaxis()->SetLabelSize(0.025);
    diff_2D_dZ.GetZaxis()->SetRangeUser(-3.0,3.0);
    diff_2D_dZ.SetStats(0);
    diff_2D_dZ.Draw("COLZ");
    c_diff_dZ.SaveAs(Form("plotdump/diff_2D_dZ_%d.png",k));
  }
  
}

std::vector<double> eFieldCalculator::cosmicToLaser(std::vector<double> inputVec){
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


std::vector<double> eFieldCalculator::laserToCosmic(std::vector<double> inputVec){
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

std::vector<double> eFieldCalculator::voxToCosmic(std::vector<int> vox){
    const double xVoxWidth = 0.1;
    const double yVoxWidth = 0.1;
    const double zVoxWidth = 0.1;
    
    std::vector<double> return_vec;
    
    return_vec.push_back(vox[0]*xVoxWidth);
    return_vec.push_back(vox[1]*xVoxWidth);
    return_vec.push_back(vox[2]*xVoxWidth);
    
    return return_vec;

}

std::vector<int> eFieldCalculator::cosmicToVox(std::vector<double> point){
    
    std::vector<int> return_vec;
    
    return_vec.push_back((int)TMath::Nint( (nCalibDivisions_x-1)*(point[0]/Lx)));
    return_vec.push_back((int)TMath::Nint( (nCalibDivisions_y-1)*(point[1]/Ly)));
    return_vec.push_back((int)TMath::Nint( (nCalibDivisions_z-1)*(point[2]/Lz)));
    
    return return_vec;

}

int main(int argc, char *argv[]){

  eFieldCalculator *calculator = new eFieldCalculator();
  calculator->compareCalib(); 

  return 0;
} //end of main

/*int main(int argc, char *argv[]){
    std::string inputCosmic = "BulkOnlyOutput.root";
    std::string inputLaser  = "RecoCorr-N3-S50-newselectbroadmuonSim-2side-Anode.root";
    std::string inputTruth  = "dispOutput_MicroBooNE_E273.root";
    std::string outputFileName = "DifferenceHistograms.root";
    TFile *cosmicFile =  new TFile(inputCosmic.c_str(),"READ");
    TFile *laserFile  =  new TFile(inputLaser.c_str(),"READ");
    TFile *truthFile  =  new TFile(inputTruth.c_str(), "READ");
    TFile *outputFile =  new TFile(outputFileName.c_str(),"RECREATE");
    
    const double driftV_offset = 2.0;
    const double metersToCM = 100.0;
    //const int nXBins = 26;    
    
    double trueDeltaX[nCalibDivisions_x][nCalibDivisions_y][nCalibDivisions_z];
    double trueDeltaY[nCalibDivisions_x][nCalibDivisions_y][nCalibDivisions_z];
    double trueDeltaZ[nCalibDivisions_x][nCalibDivisions_y][nCalibDivisions_z];
    
    for(Int_t x = 0; x < nCalibDivisions_x; x++){
      for(Int_t y = 0; y < nCalibDivisions_y; y++){
        for(Int_t z = 0; z < nCalibDivisions_z; z++){
          trueDeltaX[x][y][z] = 0.0;
          trueDeltaY[x][y][z] = 0.0;
          trueDeltaZ[x][y][z] = 0.0;
        }
     }
   }
   
   std::cout << "Reading Truth..." << std::endl;
   if(truthFile->IsZombie()){
        std::cout << "Bad input file: " << inputTruth << std::endl;
        return 1;
        
   }
    
   std::cout << "Reading..." << std::endl;
   if(cosmicFile->IsZombie()){
        std::cout << "Bad input file: " << inputCosmic << std::endl;
        return 1;
        
   }
    
    //Histograms of Difference
    
    TH2F *h_xDiffYZ[nCalibDivisions_x];
    TH2F *h_yDiffYZ[nCalibDivisions_x];
    TH2F *h_zDiffYZ[nCalibDivisions_x];
    
    TH2F *h_xDiffXY[nCalibDivisions_z];
    TH2F *h_xTruthDiffXY[nCalibDivisions_z];
    
    TH2F *h_xTruthDiffYZ[nCalibDivisions_x];
    TH2F *h_yTruthDiffYZ[nCalibDivisions_x];
    TH2F *h_zTruthDiffYZ[nCalibDivisions_x];
    
    TH2F *h_xLaser[nCalibDivisions_x];
    for (int i =0; i < nCalibDivisions_x; ++i) {
        h_xDiffYZ[i] = new TH2F(Form("xDiffYZ%d", i), "X Difference", 101, 0, 10.0, 26, 0, 2.5 );
	h_xDiffYZ[i]->SetDirectory(outputFile);
	
	h_xTruthDiffYZ[i] = new TH2F(Form("xTruthDiffYZ%d", i), "X Difference with Truth", 101, 0, 10.0, 26, 0, 2.5 );
	h_xTruthDiffYZ[i]->SetDirectory(outputFile);
	
	h_yDiffYZ[i] = new TH2F(Form("yDiffYZ%d", i), "Y Difference", 101, 0, 10.0, 26, 0, 2.5 );
	h_yDiffYZ[i]->SetDirectory(outputFile);
	
	h_zDiffYZ[i] = new TH2F(Form("zDiffYZ%d", i), "Z Difference", 101, 0, 10.0, 26, 0, 2.5 );
	h_zDiffYZ[i]->SetDirectory(outputFile);
	
	h_xLaser[i] = new TH2F(Form("xLaser%d", i), "X Laser", 101, 0, 10.0, 26, 0, 2.5 );
	h_xLaser[i]->SetDirectory(outputFile);
    }
    
    for (int i=0; i < nCalibDivisions_z; ++i){
      h_xDiffXY[i]      = new TH2F(Form("xDiffXY%d", i), "X Difference", 26, 0, 2.5, 26, 0, 2.5 );
      h_xTruthDiffXY[i] = new TH2F(Form("xTruthDiffXY%d", i), "X Difference", 26, 0, 2.5, 26, 0, 2.5 );
    
    }
    
    TTreeReader readerTruth("SpaCEtree_fwdDisp", truthFile);
    
    TTreeReaderValue<double> true_x(readerTruth, "x_true.data_fwdDisp");
    TTreeReaderValue<double> true_y(readerTruth, "y_true.data_fwdDisp");
    TTreeReaderValue<double> true_z(readerTruth, "z_true.data_fwdDisp");
    TTreeReaderValue<double> Dx_true(readerTruth, "Dx.data_fwdDisp");
    TTreeReaderValue<double> Dy_true(readerTruth, "Dy.data_fwdDisp");
    TTreeReaderValue<double> Dz_true(readerTruth, "Dz.data_fwdDisp");
    TTreeReaderValue<int> elecFate(readerTruth, "elecFate.data_fwdDisp");
    
    while(readerTruth.Next()){
       std::vector<double> truthPoint;
       truthPoint.push_back(*true_x);
       truthPoint.push_back(*true_y);
       truthPoint.push_back(*true_z);
       
       std::vector<int> voxel = cosmicToVox(truthPoint);
       
       if(0) std::cout << truthPoint[0] << " : " << voxel[0] << " " << truthPoint[1] << " : " << voxel[1] << " " << truthPoint[2] << " : " << voxel[2] << std::endl;
       
       trueDeltaX[voxel[0]][voxel[1]][voxel[2]] = *Dx_true;
       trueDeltaY[voxel[0]][voxel[1]][voxel[2]] = *Dy_true;
       trueDeltaZ[voxel[0]][voxel[1]][voxel[2]] = *Dz_true;
    
    }
    
    
    TTreeReader reader("SCEtreeBulk_calib", cosmicFile);
    
    TH3F *h_xDist;
    laserFile->GetObject("Reco_Displacement_X", h_xDist);
    TH3F *h_xDistErr;
    laserFile->GetObject("Reco_Displacement_X_Error", h_xDistErr);
    TH3F *h_yDist;
    laserFile->GetObject("Reco_Displacement_Y", h_yDist);
    TH3F *h_yDistErr;
    laserFile->GetObject("Reco_Displacement_Y_Error", h_yDistErr);
    TH3F *h_zDist;
    laserFile->GetObject("Reco_Displacement_Z", h_zDist);
    TH3F *h_zDistErr;
    laserFile->GetObject("Reco_Displacement_Z_Error", h_zDistErr);
              
    TTreeReaderValue<Double_t> Dx(reader, "Dx");
    TTreeReaderValue<Double_t> DxErr(reader, "DxErr");
    TTreeReaderValue<Double_t> Dy(reader, "Dy");
    TTreeReaderValue<Double_t> DyErr(reader, "DyErr");
    TTreeReaderValue<Double_t> Dz(reader, "Dz");
    TTreeReaderValue<Double_t> DzErr(reader, "DzErr");
    TTreeReaderValue<Double_t> x_reco(reader, "x_reco");
    TTreeReaderValue<Double_t> y_reco(reader, "y_reco");
    TTreeReaderValue<Double_t> z_reco(reader, "z_reco");
    
    
    std::cout << "Entering Loop" << std::endl;
    while (reader.Next()){
	std::vector<double> cosmicPoint;
        cosmicPoint.push_back(*x_reco);
        cosmicPoint.push_back(*y_reco);
        cosmicPoint.push_back(*z_reco);
	
	std::vector<int> cosmicVox = cosmicToVox(cosmicPoint);
        
        
        
        std::vector<double> laserVox = cosmicToLaser(cosmicPoint);
        
        int xBin = h_xDist->GetXaxis()->FindBin(laserVox[0]);
	
        int yBin = h_xDist->GetYaxis()->FindBin(laserVox[1]);
        int zBin = h_xDist->GetZaxis()->FindBin(laserVox[2]);
        
        double laserVal  = h_xDist->GetBinContent(xBin, yBin, zBin);
        double cosmicVal = *Dx;
	double truthVal  = -999.9;
	if(cosmicVox[0] < nCalibDivisions_x && cosmicVox[1] < nCalibDivisions_y && cosmicVox[2] < nCalibDivisions_z) 
	  truthVal  = ((trueDeltaX[cosmicVox[0]][cosmicVox[1]][cosmicVox[2]] + trueDeltaX[cosmicVox[0]+1][cosmicVox[1]+1][cosmicVox[2]+1]) )/ 2;
	else
	  truthVal  = trueDeltaX[cosmicVox[0]][cosmicVox[1]][cosmicVox[2]];  
	
	
        double Difference      = (laserVal - metersToCM*cosmicVal + driftV_offset);
	double LaserTruthDiff  = (metersToCM*truthVal - laserVal); 
        
	if(laserVal < 10e3){
	  h_xLaser[cosmicVox[0]]->Fill(*z_reco, *y_reco, laserVal);
	  h_xDiffYZ[cosmicVox[0]]->Fill(*z_reco, *y_reco, Difference);
	  h_xTruthDiffYZ[cosmicVox[0]]->Fill(*z_reco, *y_reco, LaserTruthDiff);
	  h_xTruthDiffXY[cosmicVox[2]]->Fill(*x_reco, *y_reco, LaserTruthDiff);
	}
	
	laserVal = h_yDist->GetBinContent(xBin, yBin, zBin);
        cosmicVal = *Dy;
	Difference = (laserVal - metersToCM*cosmicVal);
	
	if(laserVal < 10e3){
	  //h_yLaser[xBin]->Fill(*z_reco, *y_reco, laserVal);
	  h_yDiffYZ[xBin]->Fill(*z_reco, *y_reco, Difference);
	}
	
	laserVal = h_zDist->GetBinContent(xBin, yBin, zBin);
        cosmicVal = *Dz;
	Difference = (laserVal - metersToCM*cosmicVal);
	
	if(laserVal < 10e3){
	  //h_zLaser[xBin]->Fill(*z_reco, *y_reco, laserVal);
	  h_zDiffYZ[xBin]->Fill(*z_reco, *y_reco, Difference);
	}
	
	
	
	
	if(0){
	  std::cout << "Cos X: " << cosmicVox[0] << " " << "Cos Y: " << cosmicVox[1] << " " << "Cos Z: " << cosmicVox[2] << std::endl;
	  std::cout << "Las X: " << laserVox[0] << " " << "Las Y: " << laserVox[1] << " " << "Las Z: " << laserVox[2] << std::endl;
	  std::cout << "XBin: " << xBin << " " << "YBin: " << yBin << " " << "ZBin: " << zBin << std::endl;
	  std::cout << "Distortion: " << 100*cosmicVal << " "  << "Laser: " << laserVal << " " << "Diff: " << Difference << std::endl;
	}
        
	
	  
        
    }//end of loop over ttreereader
    
    outputFile->Write();
    cosmicFile->Close();
    laserFile->Close();
    outputFile->Close();
    return 0;
    
}*/
