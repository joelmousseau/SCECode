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




double calcChi2(const double *vals){
   
    const Double_t x = vals[0];
    double weight = 1 - x;
    double cosmicVal = vals[1];
    double laserVal  = vals[2];
    double truthVal  = vals[3];
    double cosmicErr = vals[4];
    double laserErr  = vals[5];
    //std::cout << weight << " " << x << std::endl;
    float numerator = pow( ( ( ( (x/cosmicErr)*cosmicVal+(weight/laserErr)*laserVal)/( (1/cosmicErr) + (1/laserErr) ) )- truthVal),2);
    float denominator = pow(truthVal, 2);
    //std::cout << numerator << " " << denominator << std::endl;
    //float numerator = pow((x+y),2);
    //float denominator = 1.0;
    
    if(denominator > 0.0)
        return numerator/denominator;
    else
        return 0.0;
    
}

void eFieldCalculator::setDriftVScale(double laser, double cosmic){
    laserDriftVScale = laser;
    cosmicDriftVScale = cosmic;
}

double eFieldCalculator::doCoordTransformX(const Double_t inputX)
{
    Double_t outputX;
    outputX = Lx - (Lx/2.56)*inputX/100.0;
    
    return outputX;
}

double eFieldCalculator::doCoordTransformY(const Double_t inputY)
{
    Double_t outputY;
    outputY = (Ly/2.33)*(inputY+116.5)/100.0;
    
    return outputY;
}

double eFieldCalculator::doCoordTransformZ(const Double_t inputZ)
{
    Double_t outputZ;
    outputZ = (Lz/10.37)*(inputZ)/100.0;
    
    return outputZ;
}

double eFieldCalculator::doInvCoordTransformX(const Double_t inputX)
{
    Double_t outputX;
    outputX = 100.0*(2.56/Lx)*(Lx - inputX);
    
    return outputX;
}

double eFieldCalculator::doInvCoordTransformY(const Double_t inputY)
{
    Double_t outputY;
    outputY = 100.0*(2.33/Ly)*inputY - 116.5;
    
    return outputY;
}

double eFieldCalculator::doInvCoordTransformZ(const Double_t inputZ)
{
    Double_t outputZ;
    outputZ = 100.0*(10.37/Lz)*inputZ;
    
    return outputZ;
}

bool eFieldCalculator::isInLaserRegion(double x, double y, double z){
    if(x < 204.83 && x >= 51.21 && y >= -51.15 && y < 51.15 && z >= 207.36 && z < 808.7)
        return true;
    else
        return false;
    
}

void eFieldCalculator::makeSmoothMap(std::string inputMapFileName, std::string outputMapFileName, bool doTriLinSmoothing, bool doEdgeSmoothing){
    
    //const bool doTriLinSmoothing = false;
    //const bool doEdgeSmoothing   = true;
    TFile* fileLaser = new TFile(inputMapFileName.c_str());
    TH3F* combine_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
    TH3F* combine_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
    TH3F* combine_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
    
    TFile *outputHistos = new TFile(outputMapFileName.c_str(), "RECREATE");
    
    TH3F* smooth_dX = new TH3F("Reco_Displacement_X","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* smooth_dY = new TH3F("Reco_Displacement_Y","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* smooth_dZ = new TH3F("Reco_Displacement_Z","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH1F* h_weightsX = new TH1F("weightsX", "Weight Applied to X Distortion", 600, -5.0, 5.0);
    TH1F* h_weightsY = new TH1F("weightsY", "Weight Applied to Y Distortion", 600, -5.0, 5.0);
    TH1F* h_weightsZ = new TH1F("weightsZ", "Weight Applied to Z Distortion", 600, -5.0, 5.0);
    
    TH1F* h_diffX = new TH1F("diffX", "Difference Across X Boundary", 100, -1.0, 1.0);
    
    int lowX = 6;
    int highX  = 21;
    int lowY = 8;
    int highY  = 19;
    int lowZ = 21;
    int highZ  = 79;
    int fitFailures = 0;
    for(int i=1; i <= combine_dX->GetNbinsX(); ++i){
        for(int j=1; j <= combine_dX->GetNbinsY(); ++j){
            for(int k=1; k <= combine_dX->GetNbinsZ(); ++k){
                bool atXLowerBoundary = false;
                bool atYLowerBoundary = false;
                bool atZLowerBoundary = false;
                bool atXUpperBoundary = false;
                bool atYUpperBoundary = false;
                bool atZUpperBoundary = false;
                
                atXLowerBoundary = (i == lowX && j > lowY && j < highY && k > lowZ && k < highZ);
                atYLowerBoundary = (j == lowY && i > lowX && i < highX && k > lowZ && k < highZ);
                atZLowerBoundary = (k == lowZ && i > lowX && i < highX && j > lowY && j < highY);
                
                atXUpperBoundary = (i == highX && j > lowY && j < highY && k > lowZ && k < highZ);
                atYUpperBoundary = (j == highY && i > lowX && i < highX && k > lowZ && k < highZ);
                atZUpperBoundary = (k == highZ && i > lowX && i < highX && j > lowY && j < highY);
                
                if(doEdgeSmoothing && atXLowerBoundary){
                    //std::cout << i << std::endl;
                    
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i+1);
                    pointTwo.push_back(j);
                    pointTwo.push_back(k);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::xAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::xAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::xAxis);
                    
                   // std::cout << "Diff Before: " << (combine_dX->GetBinContent(i,j,k) - combine_dX->GetBinContent(i+1,j,k) ) << std::endl;
                   // std::cout << "Diff After: "  << ( << std::endl;
                    
                                                     h_diffX->Fill(weightX*combine_dX->GetBinContent(i,j,k) - combine_dX->GetBinContent(i+1,j,k) );
                    
                   float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                   float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                   float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                   
                   smooth_dX->SetBinContent(i,j,k, smoothX);
                   smooth_dY->SetBinContent(i,j,k, smoothY);
                   smooth_dZ->SetBinContent(i,j,k, smoothZ);
                   if(smoothX > 100.0)
                       ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                else if(doEdgeSmoothing && atXUpperBoundary){
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i-1);
                    pointTwo.push_back(j);
                    pointTwo.push_back(k);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::xAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::xAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::xAxis);
                    
                    float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                    float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                    float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                    if(smoothX > 100.0)
                        ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                else if(doEdgeSmoothing && atYLowerBoundary){
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i);
                    pointTwo.push_back(j+1);
                    pointTwo.push_back(k);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::yAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::yAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::yAxis);
                    
                    float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                    float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                    float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                    if(smoothX > 100.0)
                        ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                else if(doEdgeSmoothing && atYUpperBoundary){
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i);
                    pointTwo.push_back(j-1);
                    pointTwo.push_back(k);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::yAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::yAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::yAxis);
                    
                    float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                    float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                    float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                    if(smoothX > 100.0)
                        ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                else if(doEdgeSmoothing && atZLowerBoundary){
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i);
                    pointTwo.push_back(j);
                    pointTwo.push_back(k+1);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::zAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::zAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::zAxis);
                    
                    float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                    float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                    float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                    if(smoothX > 100.0)
                        ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                else if(doEdgeSmoothing && atZUpperBoundary){
                    std::vector<int> pointOne;
                    pointOne.push_back(i);
                    pointOne.push_back(j);
                    pointOne.push_back(k);
                    
                    std::vector<int> pointTwo;
                    pointTwo.push_back(i);
                    pointTwo.push_back(j);
                    pointTwo.push_back(k-1);
                    
                    float weightX = SmoothBoundary(combine_dX, combine_dX, pointOne, pointTwo, axisType::zAxis);
                    float weightY = SmoothBoundary(combine_dY, combine_dY, pointOne, pointTwo, axisType::zAxis);
                    float weightZ = SmoothBoundary(combine_dZ, combine_dZ, pointOne, pointTwo, axisType::zAxis);
                    
                    float smoothX = weightX*combine_dX->GetBinContent(i,j,k);
                    float smoothY = weightY*combine_dY->GetBinContent(i,j,k);
                    float smoothZ = weightZ*combine_dZ->GetBinContent(i,j,k);
                    
                    h_weightsX->Fill(weightX);
                    h_weightsY->Fill(weightY);
                    h_weightsZ->Fill(weightZ);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                    if(smoothX > 100.0)
                        ++fitFailures;
                    if(smoothY > 100.0)
                        ++fitFailures;
                    if(smoothZ > 100.0)
                        ++fitFailures;
                }
                
                
                else if(doTriLinSmoothing && (i < highX+1 && i > lowX-1 && k > lowZ-1 && k < highZ+1 && j > lowY-1 && j < highY+1) ){
                    std::vector<int> point;
                    point.push_back(i);
                    point.push_back(j);
                    point.push_back(k);
                    
                    float smoothX = TrilinInterp(combine_dX, point);
                    float smoothY = TrilinInterp(combine_dY, point);
                    float smoothZ = TrilinInterp(combine_dZ, point);
                    
                    smooth_dX->SetBinContent(i,j,k, smoothX);
                    smooth_dY->SetBinContent(i,j,k, smoothY);
                    smooth_dZ->SetBinContent(i,j,k, smoothZ);
                    
                   
                }
                
                else{
                    smooth_dX->SetBinContent(i,j,k, combine_dX->GetBinContent(i,j,k));
                    smooth_dY->SetBinContent(i,j,k, combine_dY->GetBinContent(i,j,k));
                    smooth_dZ->SetBinContent(i,j,k, combine_dZ->GetBinContent(i,j,k));
                }
            }
        }
    }
    
    /*
    for(int i =1; i < startX; i++){
     
        for(int j = 1; j < startY; j++){
            for(int k = 1; k < startZ; k++){
                //std::cout << i << " " << j << " " << k << " " <<  combine_dX->GetBinContent(i,j,k) << std::endl;
                smooth_dX->SetBinContent(i,j,k, combine_dX->GetBinContent(i,j,k));
                smooth_dY->SetBinContent(i,j,k, combine_dY->GetBinContent(i,j,k));
                smooth_dZ->SetBinContent(i,j,k, combine_dZ->GetBinContent(i,j,k));
     
            }
        }
    }
    
    for(int i = startX; i <= stopX; i++){
        
        for(int j = startY; j <= stopY; j++){
            for(int k = startZ; k <= stopZ; k++){
                std::vector<int> point;
                point.push_back(i);
                point.push_back(j);
                point.push_back(k);
                
                float smoothX = TrilinInterp(combine_dX, point);
                float smoothY = TrilinInterp(combine_dY, point);
                float smoothZ = TrilinInterp(combine_dZ, point);
                
                smooth_dX->SetBinContent(i,j,k, smoothX);
                smooth_dY->SetBinContent(i,j,k, smoothY);
                smooth_dZ->SetBinContent(i,j,k, smoothZ);
                
                
            }
        }
    }
    
    for(int i =stopX+1; i <= smooth_dX->GetNbinsX(); i++){
     
        for(int j = stopY+1; j <= smooth_dX->GetNbinsY(); j++){
            for(int k = stopZ+1; k <= smooth_dX->GetNbinsZ(); k++){
                smooth_dX->SetBinContent(i,j,k, combine_dX->GetBinContent(i,j,k));
                smooth_dY->SetBinContent(i,j,k, combine_dY->GetBinContent(i,j,k));
                smooth_dZ->SetBinContent(i,j,k, combine_dZ->GetBinContent(i,j,k));
            }
        }
    }
    */
   
    for(int k = 1; k <= smooth_dX->GetNbinsZ(); k++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
            for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
                combine_2D_dX.SetBinContent(i,j,smooth_dX->GetBinContent(i,j,k));
                combine_2D_dY.SetBinContent(i,j,smooth_dY->GetBinContent(i,j,k));
                combine_2D_dZ.SetBinContent(i,j,smooth_dZ->GetBinContent(i,j,k));
                
                
            }
        }
        drawPlanarPlot(combine_2D_dX, k, "Smoothed #DeltaX", "smooth_2D_dX", axisType::zAxis, 5.0);
        drawPlanarPlot(combine_2D_dY, k, "Smoothed #DeltaY", "smooth_2D_dY", axisType::zAxis);
        drawPlanarPlot(combine_2D_dZ, k, "Smoothed #DeltaZ", "smooth_2D_dZ", axisType::zAxis);
    }
    
    for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                //if(i == 5 && j == 7 && k == 20)
                //    std::cout << i << " " << j << " " << k << " " <<  smooth_dX->GetBinContent(i,j,k) << std::endl;
                combine_2D_dX.SetBinContent(k,i,smooth_dX->GetBinContent(i,j,k));
                combine_2D_dY.SetBinContent(k,i,smooth_dY->GetBinContent(i,j,k));
                combine_2D_dZ.SetBinContent(k,i,smooth_dZ->GetBinContent(i,j,k));
                
                
            }
        }
        drawPlanarPlot(combine_2D_dX, j, "Smoothed #DeltaX", "smooth_2D_dX", axisType::yAxis, 5.0);
        drawPlanarPlot(combine_2D_dY, j, "Smoothed #DeltaY", "smooth_2D_dY", axisType::yAxis);
        drawPlanarPlot(combine_2D_dZ, j, "Smoothed #DeltaZ", "smooth_2D_dZ", axisType::yAxis);
    }
    
    TCanvas *can = new TCanvas(Form("can_orig"),"",600,600);
    
    gStyle->SetTitleW(0.9);
    gStyle->SetOptStat(0);
    
    h_weightsX->SetLineWidth(2.0);
    h_weightsX->SetLineColor(kBlack);
    h_weightsX->GetXaxis()->SetTitle("Weight Applied");
    std::cout << "X Weights: " << h_weightsX->Integral() << std::endl;
    h_weightsX->GetXaxis()->SetRangeUser(-1.0, 3.0);
    h_weightsX->Draw("hist");
    can->SaveAs("WeightsOnX.png");
    can->Update();
    can->Clear();
    
    h_weightsY->SetLineWidth(2.0);
    h_weightsY->SetLineColor(kBlack);
    h_weightsY->GetXaxis()->SetTitle("Weight Applied");
    std::cout << "Y Weights: " << h_weightsY->Integral() << std::endl;
    //h_weightsY->Rebin(10);
    h_weightsY->GetXaxis()->SetRangeUser(-1.0, 3.0);
    h_weightsY->Draw("hist");
    can->SaveAs("WeightsOnY.png");
    can->Update();
    can->Clear();
    
    h_weightsZ->SetLineWidth(2.0);
    h_weightsZ->SetLineColor(kBlack);
    h_weightsZ->GetXaxis()->SetTitle("Weight Applied");
    std::cout << "Z Weights: " << h_weightsZ->Integral() << std::endl;
    //h_weightsZ->Rebin(10);
    h_weightsZ->GetXaxis()->SetRangeUser(-1.0, 3.0);
    h_weightsZ->Draw("hist");
    can->SaveAs("WeightsOnZ.png");
    can->Update();
    can->Clear();
    
    
    h_diffX->SetLineWidth(2.0);
    h_diffX->SetLineColor(kBlack);
    h_diffX->GetXaxis()->SetTitle("Difference Across X Boundary (cm)");
    //std::cout << "Z Weights: " << h_weightsZ->Integral() << std::endl;
    //h_weightsZ->Rebin(10);
    h_diffX->GetXaxis()->SetRangeUser(-0.2, 0.2);
    h_diffX->Draw("hist");
    can->SaveAs("DiffX.png");
    can->Update();
    can->Clear();
    

    
    outputHistos->Write();
    outputHistos->Close();
    std::cout << "Fit Failures: " << fitFailures << std::endl;
}

void eFieldCalculator::makeSmoothMapPlots(){
    double stops[5] = {0.00,0.34,0.61,0.84,1.00};
    double red[5] = {0.00,0.00,0.87,1.00,0.51};
    double green[5] = {0.00,0.81,1.00,0.20,0.00};
    double blue[5] = {0.51,1.00,0.12,0.00,0.00};
    int   zCuts[2] = {9, 17};
    TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
    gStyle->SetNumberContours(255);
    
    TFile* fileData = new TFile("MergedMapsSmoothCosmicAndLaserV1098Median.root");
    TH3F* fwd_data_dX = (TH3F*) fileData->Get("combined_fwd_dX");
    TH3F* fwd_data_dY = (TH3F*) fileData->Get("combined_fwd_dY");
    TH3F* fwd_data_dZ = (TH3F*) fileData->Get("combined_fwd_dZ");
    
    for(Int_t k = 1; k <= fwd_data_dX->GetNbinsZ(); k++) {
        TH2F data_2D_dX(Form("dataFwd_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F data_2D_dY(Form("dataFwd_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F data_2D_dZ(Form("dataFwd_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        for(Int_t i = 1; i <= fwd_data_dX->GetNbinsX()-1; i++)
        {
            for(Int_t j = 1; j <= fwd_data_dX->GetNbinsY(); j++)
            {
                
               
               data_2D_dX.SetBinContent(i,j,fwd_data_dX->GetBinContent(i,j,k));
               data_2D_dY.SetBinContent(i,j,fwd_data_dY->GetBinContent(i,j,k));
               data_2D_dZ.SetBinContent(i,j,fwd_data_dZ->GetBinContent(i,j,k));
                
                
            }
        }
        
        drawPlanarPlot(data_2D_dX, k, "Data dX (smooth)", "data_smooth_dX", axisType::zAxis, 5.0);
        drawPlanarPlot(data_2D_dY, k, "Data dY (smooth)", "data_smooth_dY", axisType::zAxis, 5.0);
        drawPlanarPlot(data_2D_dZ, k, "Data dZ (smooth)", "data_smooth_dZ", axisType::zAxis, 5.0);
        
        
    }
    
    for(Int_t k = 1; k <= fwd_data_dX->GetNbinsY(); k++)
    {
        TH2F data_2D_dX(Form("dataFwd_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F data_2D_dY(Form("dataFwd_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F data_2D_dZ(Form("dataFwd_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        
        for(Int_t i = 1; i <= fwd_data_dX->GetNbinsX()-1; i++)
        {
            for(Int_t j = 1; j <= fwd_data_dX->GetNbinsZ(); j++)
            {
                
                
                
              data_2D_dX.SetBinContent(j,i,fwd_data_dX->GetBinContent(i,k,j));
              data_2D_dY.SetBinContent(j,i,fwd_data_dY->GetBinContent(i,k,j));
              data_2D_dZ.SetBinContent(j,i,fwd_data_dZ->GetBinContent(i,k,j));
                
                
            }
        }
        
        drawPlanarPlot(data_2D_dX, k, "Data dX (smooth)", "data_smooth_dX", axisType::yAxis, 5.0);
        drawPlanarPlot(data_2D_dY, k, "Data dY (smooth)", "data_smooth_dY", axisType::yAxis, 5.0);
        drawPlanarPlot(data_2D_dZ, k, "Data dZ (smooth)", "data_smooth_dZ", axisType::yAxis, 5.0);
        
        
    }
    
}

void eFieldCalculator::makeFwdMapPlots(){
    TFile* fileData = new TFile("/uboone/data/users/joelam/SCEInputFiles/TrueDist-N3-S50_laserdata_v1098_fwd.root");
    TH3F* fwd_data_dX = (TH3F*) fileData->Get("Reco_Displacement_X");
    TH3F* fwd_data_dY = (TH3F*) fileData->Get("Reco_Displacement_Y");
    TH3F* fwd_data_dZ = (TH3F*) fileData->Get("Reco_Displacement_Z");
    
    TFile* fileInv = new TFile("/uboone/data/users/joelam/SCEDistortionMaps/MergedFwdMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.root");
    TH3F* fwd_inv_dX = (TH3F*) fileInv->Get("combined_fwd_dX");
    TH3F* fwd_inv_dY = (TH3F*) fileInv->Get("combined_fwd_dY");
    TH3F* fwd_inv_dZ = (TH3F*) fileInv->Get("combined_fwd_dZ");
    
    TH1F *h_dX_diff = new TH1F("dX_diff", "Difference Between Inverted - Measured Fwd. (#Delta X) Map", 1000, -5.0, 5.0);
    TH1F *h_dY_diff = new TH1F("dY_diff", "Difference Between Inverted - Measured Fwd. (#Delta Y) Map", 1000, -5.0, 5.0);
    TH1F *h_dZ_diff = new TH1F("dZ_diff", "Difference Between Inverted - Measured Fwd. (#Delta Z) Map", 1000, -5.0, 5.0);
    
    /*const int startX = 6;
    const int endX   = 21;
    const int startY = 8;
    const int endY   = 19;
    const int startZ = 21;
    const int endZ   = 71;*/
    
    const int startX = 1;
    const int endX   = fwd_inv_dX->GetNbinsX();
    const int startY = 1;
    const int endY   = fwd_inv_dX->GetNbinsY();
    const int startZ = 1;
    const int endZ   = fwd_inv_dX->GetNbinsZ();
    
    for(Int_t k = startZ; k <= endZ; k++) {
        TH2F data_2D_dX(Form("dataFwd_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F data_2D_dY(Form("dataFwd_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F data_2D_dZ(Form("dataFwd_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        TH2F inv_2D_dX(Form("invFwd_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F inv_2D_dY(Form("invFwd_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F inv_2D_dZ(Form("invFwd_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        for(Int_t i = startX; i <= endX; i++)
        {
            for(Int_t j = startY; j <= endY; j++)
            {
                
                
                data_2D_dX.SetBinContent(i,j,1.0*fwd_data_dX->GetBinContent(i,j,k));
                data_2D_dY.SetBinContent(i,j,1.0*fwd_data_dY->GetBinContent(i,j,k));
                data_2D_dZ.SetBinContent(i,j,1.0*fwd_data_dZ->GetBinContent(i,j,k));
                
                inv_2D_dX.SetBinContent(i,j,1.0*fwd_inv_dX->GetBinContent(i,j,k));
                inv_2D_dY.SetBinContent(i,j,1.0*fwd_inv_dY->GetBinContent(i,j,k));
                inv_2D_dZ.SetBinContent(i,j,1.0*fwd_inv_dZ->GetBinContent(i,j,k));
                
                h_dX_diff->Fill(fwd_inv_dX->GetBinContent(i,j,k) - fwd_data_dX->GetBinContent(i,j,k));
                h_dY_diff->Fill(fwd_inv_dY->GetBinContent(i,j,k) - fwd_data_dY->GetBinContent(i,j,k));
                h_dZ_diff->Fill(fwd_inv_dZ->GetBinContent(i,j,k) - fwd_data_dZ->GetBinContent(i,j,k));
                
                
            }
        }
        
        drawPlanarPlot(data_2D_dX, k, "Data dX (fwd)", "data_Fwd_dX", axisType::zAxis, 5.0);
        drawPlanarPlot(data_2D_dY, k, "Data dY (fwd)", "data_Fwd_dY", axisType::zAxis, 5.0);
        drawPlanarPlot(data_2D_dZ, k, "Data dZ (fwd)", "data_Fwd_dZ", axisType::zAxis, 5.0);
        
        drawPlanarPlot(inv_2D_dX, k, "Inverted dX (fwd)", "inv_Fwd_dX", axisType::zAxis, 5.0);
        drawPlanarPlot(inv_2D_dY, k, "Inverted dY (fwd)", "inv_Fwd_dY", axisType::zAxis, 5.0);
        drawPlanarPlot(inv_2D_dZ, k, "Inverted dZ (fwd)", "inv_Fwd_dZ", axisType::zAxis, 5.0);
        
        
    }
    
    for(Int_t k = startY; k <= endY; k++)
    {
        TH2F data_2D_dX(Form("dataFwd_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F data_2D_dY(Form("dataFwd_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F data_2D_dZ(Form("dataFwd_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        TH2F inv_2D_dX(Form("invFwd_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F inv_2D_dY(Form("invFwd_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F inv_2D_dZ(Form("invFwd_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        
        for(Int_t i = startX; i <= endX; i++)
        {
            for(Int_t j = startZ; j <= endZ; j++)
            {
                
                
                
                data_2D_dX.SetBinContent(j,i,1.0*fwd_data_dX->GetBinContent(i,k,j));
                data_2D_dY.SetBinContent(j,i,1.0*fwd_data_dY->GetBinContent(i,k,j));
                data_2D_dZ.SetBinContent(j,i,1.0*fwd_data_dZ->GetBinContent(i,k,j));
                
                inv_2D_dX.SetBinContent(j,i,1.0*fwd_inv_dX->GetBinContent(i,k,j));
                inv_2D_dY.SetBinContent(j,i,1.0*fwd_inv_dY->GetBinContent(i,k,j));
                inv_2D_dZ.SetBinContent(j,i,1.0*fwd_inv_dZ->GetBinContent(i,k,j));
                
                
            }
        }
        
        drawPlanarPlot(data_2D_dX, k, "Data dX (fwd)", "data_Fwd_dX", axisType::yAxis, 5.0);
        drawPlanarPlot(data_2D_dY, k, "Data dY (fwd)", "data_Fwd_dY", axisType::yAxis);
        drawPlanarPlot(data_2D_dZ, k, "Data dZ (fwd)", "data_Fwd_dZ", axisType::yAxis);
        
        drawPlanarPlot(inv_2D_dX, k, "Inverted dX (fwd)", "inv_Fwd_dX", axisType::yAxis, 5.0);
        drawPlanarPlot(inv_2D_dY, k, "Inverted dY (fwd)", "inv_Fwd_dY", axisType::yAxis);
        drawPlanarPlot(inv_2D_dZ, k, "Inverted dZ (fwd)", "inv_Fwd_dZ", axisType::yAxis);
        
        
    }
    TCanvas *can = new TCanvas(Form("can_orig"),"",600,600);
    
    gStyle->SetTitleW(0.9);
    gStyle->SetOptStat("mr");
    
    h_dX_diff->SetLineWidth(2.0);
    h_dX_diff->SetLineColor(kBlack);
    h_dX_diff->GetXaxis()->SetTitle("Difference (cm)");
    //std::cout << "X Weights: " << h_weightsX->Integral() << std::endl;
    h_dX_diff->GetXaxis()->SetRangeUser(-1.0, 1.0);
    h_dX_diff->Draw("hist");
    can->SaveAs("FwdDataInvDiffdX.png");
    can->Update();
    can->Clear();
    
    h_dY_diff->SetLineWidth(2.0);
    h_dY_diff->SetLineColor(kBlack);
    h_dY_diff->GetXaxis()->SetTitle("Difference (cm)");
    //std::cout << "X Weights: " << h_weightsX->Integral() << std::endl;
    h_dY_diff->GetXaxis()->SetRangeUser(-1.0, 1.0);
    h_dY_diff->Draw("hist");
    can->SaveAs("FwdDataInvDiffdY.png");
    can->Update();
    can->Clear();
    
    h_dZ_diff->SetLineWidth(2.0);
    h_dZ_diff->SetLineColor(kBlack);
    h_dZ_diff->GetXaxis()->SetTitle("Difference (cm)");
    //std::cout << "X Weights: " << h_weightsX->Integral() << std::endl;
    h_dZ_diff->GetXaxis()->SetRangeUser(-1.0, 1.0);
    h_dZ_diff->Draw("hist");
    can->SaveAs("FwdDataInvDiffdZ.png");
    can->Update();
    can->Clear();
    

}

void eFieldCalculator::compareCalib(bool isData)
{
  
  const double xScale = TPC_X/Lx;
  const double yScale = TPC_Y/Ly;
  const double zScale = TPC_Z/Lz;
  const int   zRegions = 3; //upstream, middle, downstream
  int       zRegion = 0;
  const double zMaximum = 3.0;

  /*double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};*/
  int   zCuts[2] = {35, 65};
  /*TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);*/
  
  std::string inputLaser;
  std::string inputCosmic;
  if(isData){
     inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50_laserdata_v1098_bkwd.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_data_200k_Aug3_smoothed.root";
     
  }
  
  else{
     inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_MC_200k_Aug3_smoothed.root";
  }  

  TFile* fileLaser = new TFile(inputLaser.c_str());
  TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
  TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
  TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
  TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
  TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
  TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
  
  TFile* fileCosmic = new TFile(inputCosmic.c_str());
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
  
  TFile* fileGoodVoxels = new TFile("GoodVoxels.root", "RECREATE");

  TH3F* diff_dX = new TH3F("diff_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dY = new TH3F("diff_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dZ = new TH3F("diff_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* good_dX = new TH3F("good_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* good_dY = new TH3F("good_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* good_dZ = new TH3F("good_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* cosine_angles  = new TH3F("cosine_angles","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
  TH3F* magnitude_cosmic  = new TH3F("magnitude_cosmic","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* magnitude_laser  = new TH3F("magnitude_laser","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
  TH3F* magnitude_diff  = new TH3F("magnitude_diff","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);

    
  bool plotX = false;
  bool plotY = false;
  bool plotZ = false;
  double driftV = 0.0;
  if(isData)
    driftV = 0.001;
  driftV = driftSign*driftV;
    
    
  for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
  {
    for(Int_t j = 1; j <= diff_dX->GetNbinsY(); j++)
    {
      for(Int_t k = 1; k <= diff_dX->GetNbinsZ(); k++)
      {
	
	plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
	
	std::vector<double> cosmic_Ds;
	std::vector<double> laser_Ds;
	
	
	double cosmic_x_correction = -i*driftV*cosmicDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
    double laser_x_correction =  i*driftV*laserDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
	
	cosmic_Ds.push_back(cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction);
	cosmic_Ds.push_back(cosmic_dY->GetBinContent(i,j,k));
	cosmic_Ds.push_back(cosmic_dZ->GetBinContent(i,j,k));
	
	laser_Ds.push_back(laser_dX->GetBinContent(i,j,k) + laser_x_correction);
	laser_Ds.push_back(laser_dY->GetBinContent(i,j,k));
	laser_Ds.push_back(laser_dZ->GetBinContent(i,j,k));
	
	double cosine    = getAngle(cosmic_Ds, laser_Ds);
    double laserMag  = getVectorMagnitude(laser_Ds);
    double cosmicMag = getVectorMagnitude(cosmic_Ds);
    double magDiff   = laserMag - cosmicMag;
    if(goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)))
          magnitude_laser->SetBinContent(i,j,k,laserMag);
	
    if(goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) )
          magnitude_cosmic->SetBinContent(i,j,k,cosmicMag);
          
	if(plotX)
	{
          diff_dX->SetBinContent(i,j,k,(laser_dX->GetBinContent(i,j,k) + laser_x_correction)-(cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction) );
	}

	if(plotY)
	{
          diff_dY->SetBinContent(i,j,k,laser_dY->GetBinContent(i,j,k)-cosmic_dY->GetBinContent(i,j,k));
	}

	if(plotZ)
	{
          diff_dZ->SetBinContent(i,j,k,laser_dZ->GetBinContent(i,j,k)-cosmic_dZ->GetBinContent(i,j,k));
	}
	
	if(plotX && plotY && plotZ){
	  cosine_angles->SetBinContent(i,j,k,cosine);
      
	  
      magnitude_diff->SetBinContent(i,j,k,magDiff);
        //std::cout << cosine << std::endl;
	}  
      }
    }
  }

  gStyle->SetTitleW(0.9);
  
  TH1F *h_diffXByRegion[zRegions];
  TH1F *h_diffYByRegion[zRegions];
  TH1F *h_diffZByRegion[zRegions];
  
  plotX = false;
  plotY = false;
  plotZ = false;
  
  for(int k = 0; k < zRegions; ++k){
      h_diffXByRegion[k] = new TH1F(Form("diff_1D_dX_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
      h_diffYByRegion[k] = new TH1F(Form("diff_1D_dY_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
      h_diffZByRegion[k] = new TH1F(Form("diff_1D_dZ_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
  
  }

  for(int k = 1; k <= diff_dX->GetNbinsZ(); k++)
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
    
    TH2F cosine_2D(Form("cosine_2D_%d", k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F magnitude_cosmic_2D(Form("magnitude_cosmic_2D_%d", k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F magnitude_laser_2D(Form("magnitude_laser_2D_%d", k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F magnitude_diff_2D(Form("magnitude_diff_2D_%d", k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
 
    if( k < zCuts[0])
      zRegion = 0;
    else if( k >= zCuts[0] && k < zCuts[1])
      zRegion = 1;
    else if( k >= zCuts[1])
      zRegion = 2;
    else
      std::cout << "Cannot determine zRegion" << std::endl;  
          

    for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
    {
      for(Int_t j = 1; j <= diff_dX->GetNbinsY(); j++)
      {
          double cosmic_x_correction = -i*driftV*cosmicDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
          double laser_x_correction  = -i*driftV*laserDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);

          plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
    
    if(goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) )
              magnitude_laser_2D.SetBinContent(i,j, magnitude_laser->GetBinContent(i,j,k));
          
    if(goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) )
              magnitude_cosmic_2D.SetBinContent(i,j, magnitude_cosmic->GetBinContent(i,j,k));
    laser_2D_dX.SetBinContent(i,j,laser_dX->GetBinContent(i,j,k)   + laser_x_correction);
    cosmic_2D_dX.SetBinContent(i,j,cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction);
	
    if(plotX){
	  
	  
	  h_diffXByRegion[zRegion]->Fill(diff_dX->GetBinContent(i,j,k));
        
      if(zRegion == 1 && fabs(diff_dX->GetBinContent(i,j,k)) <= maxDiff)
	    good_dX->SetBinContent(i,j,k,1.0);
	}
    laser_2D_dY.SetBinContent(i,j,laser_dY->GetBinContent(i,j,k));
    cosmic_2D_dY.SetBinContent(i,j,cosmic_dY->GetBinContent(i,j,k));
	
          if(plotY){
	  
      
	  h_diffYByRegion[zRegion]->Fill(diff_dY->GetBinContent(i,j,k));
	  if(zRegion == 1 && fabs(diff_dY->GetBinContent(i,j,k)) <= maxDiff)
	    good_dY->SetBinContent(i,j,k,1.0);
	}
    laser_2D_dZ.SetBinContent(i,j,laser_dZ->GetBinContent(i,j,k));
    cosmic_2D_dZ.SetBinContent(i,j,cosmic_dZ->GetBinContent(i,j,k));
          
          
	if(plotZ){
	  
	  h_diffZByRegion[zRegion]->Fill(diff_dZ->GetBinContent(i,j,k));
	  if(zRegion == 1 && fabs(diff_dZ->GetBinContent(i,j,k)) <= maxDiff)
	    good_dZ->SetBinContent(i,j,k,1.0);
	
	
	}
	
	if(plotX && plotY && plotZ){
	   cosine_2D.SetBinContent(i,j, cosine_angles->GetBinContent(i,j,k));
       magnitude_diff_2D.SetBinContent(i,j, magnitude_diff->GetBinContent(i,j,k));
        
    }
	
        diff_2D_dX.SetBinContent(i,j,diff_dX->GetBinContent(i,j,k));
        diff_2D_dY.SetBinContent(i,j,diff_dY->GetBinContent(i,j,k));
        diff_2D_dZ.SetBinContent(i,j,diff_dZ->GetBinContent(i,j,k));
          
	  
      }
    }

    
    
   drawPlanarPlot(laser_2D_dX, k, "Laser #DeltaX", "laser_2D_dX", axisType::zAxis, 5.0);
   drawPlanarPlot(laser_2D_dY, k, "Laser #DeltaY", "laser_2D_dY", axisType::zAxis);
   drawPlanarPlot(laser_2D_dZ, k, "Laser #DeltaZ", "laser_2D_dZ", axisType::zAxis);
   drawPlanarPlot(cosmic_2D_dX, k, "Cosmic #DeltaX", "cosmic_2D_dX", axisType::zAxis, 5.0);
   drawPlanarPlot(cosmic_2D_dY, k, "Cosmic #DeltaY", "cosmic_2D_dY", axisType::zAxis);
   drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic #DeltaZ", "cosmic_2D_dZ", axisType::zAxis);
   drawPlanarPlot(diff_2D_dX, k, "Laser-Cosmic #DeltaX", "diff_2D_dX", axisType::zAxis, zMaximum);
   drawPlanarPlot(diff_2D_dY, k, "Laser-Cosmic #DeltaY", "diff_2D_dY", axisType::zAxis, zMaximum);
   drawPlanarPlot(diff_2D_dZ, k, "Laser-Cosmic #DeltaZ", "diff_2D_dZ", axisType::zAxis, zMaximum);
   drawPlanarPlot(cosine_2D,  k, "Cosine of Laser - Cosmic", "cos_2D", axisType::zAxis, 1.0);
   drawPlanarPlot(magnitude_cosmic_2D,  k, "Magnitude of Cosmic Deflection", "cosmic_2D_mag", axisType::zAxis, 10.0);
   drawPlanarPlot(magnitude_laser_2D,  k, "Magnitude of Laser Deflection", "laser_2D_mag", axisType::zAxis, 10.0);
   drawPlanarPlot(magnitude_diff_2D,  k, "Differece in Magnitude", "diff_2D_mag", axisType::zAxis, 8.0);
      
   

    

    
   

    
    
    
  }
  double zOne = (((double) zCuts[0] - 1)/100.0)*TPC_Z;
  double zTwo = (((double) zCuts[1] - 1)/100.0)*TPC_Z;
  double maximum = -999.9;
  //std::vector<int> binRange = calculateFWHMBins(h_diffXByRegion[0]);
  //if(binRange.size() != 0)
  //  std::cout << h_diffXByRegion[0]->GetBinCenter(binRange[0])  << " " << h_diffXByRegion[0]->GetBinCenter(binRange[1]) << std::endl;
  h_diffXByRegion[0]->Scale( TPC_Z / zOne);
  h_diffXByRegion[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffXByRegion[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffXByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[0]->GetMaximum();
  if(h_diffXByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[1]->GetMaximum();
  if(h_diffXByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[2]->GetMaximum();      
  
  h_diffXByRegion[0]->SetMaximum(1.1*maximum);
  h_diffXByRegion[1]->SetMaximum(1.1*maximum);
  h_diffXByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dX("c_1D_diff_dX", "", 600, 600);
  c_1D_diff_dX.cd();
  TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.55);
  legend->SetTextSize(0.03);
  h_diffXByRegion[0]->SetLineColor(kBlack);  
  legend->AddEntry(h_diffXByRegion[0], Form("Z < %.2f", zOne), "l");
  h_diffXByRegion[1]->SetLineColor(kRed+2);  
  legend->AddEntry(h_diffXByRegion[1], Form("%.2f #leq Z < %.2f", zOne, zTwo), "l");
  h_diffXByRegion[2]->SetLineColor(kBlue+2);  
  legend->AddEntry(h_diffXByRegion[2], Form("%.2f #leq Z", zTwo), "l");
  
  for(int i=0; i < zRegions; ++i){
    h_diffXByRegion[i]->SetTitle("Laser-Cosmic #DeltaX [cm]");
    h_diffXByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffXByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffXByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffXByRegion[i]->Draw("hist");
    else  
      h_diffXByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dX.SaveAs("1D_plots/XY_diff_1D_dX.png");
  
  maximum = -999.9;
  h_diffYByRegion[0]->Scale( TPC_Z / zOne);
  h_diffYByRegion[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffYByRegion[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffYByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[0]->GetMaximum();
  if(h_diffYByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[1]->GetMaximum();
  if(h_diffYByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[2]->GetMaximum();      
  
  h_diffYByRegion[0]->SetMaximum(1.1*maximum);
  h_diffYByRegion[1]->SetMaximum(1.1*maximum);
  h_diffYByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dY("c_1D_diff_dY", "", 600, 600);
  c_1D_diff_dY.cd();
  h_diffYByRegion[0]->SetLineColor(kBlack);  
  h_diffYByRegion[1]->SetLineColor(kRed+2);  
  h_diffYByRegion[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffYByRegion[i]->SetTitle("Laser-Cosmic #DeltaY [cm]");
    h_diffYByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffYByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffYByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffYByRegion[i]->Draw("hist");
    else  
      h_diffYByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dY.SaveAs("1D_plots/XY_diff_1D_dY.png");
  
  maximum = -999.9;
  h_diffZByRegion[0]->Scale( TPC_Z / zOne);
  h_diffZByRegion[1]->Scale( TPC_Z / (zTwo - zOne)  );
  h_diffZByRegion[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffZByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[0]->GetMaximum();
  if(h_diffZByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[1]->GetMaximum();
  if(h_diffZByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[2]->GetMaximum();
   
  h_diffZByRegion[0]->SetMaximum(1.1*maximum);
  h_diffZByRegion[1]->SetMaximum(1.1*maximum);
  h_diffZByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dZ("c_1D_diff_dZ", "", 600, 600);
  c_1D_diff_dZ.cd();
  h_diffZByRegion[0]->SetLineColor(kBlack);  
  h_diffZByRegion[1]->SetLineColor(kRed+2);  
  h_diffZByRegion[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffZByRegion[i]->SetTitle("Laser-Cosmic #DeltaZ [cm]");
    h_diffZByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffZByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffZByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffZByRegion[i]->Draw("hist");
    else  
      h_diffZByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dZ.SaveAs("1D_plots/XY_diff_1D_dZ.png");
  
  TCanvas c_1D_goodVoxels_dX("c_godVoxels_dX", "", 1300, 900);
  c_1D_goodVoxels_dX.cd();
  good_dX->GetXaxis()->SetTitle("X [cm]");
  good_dX->GetXaxis()->SetTitleOffset(1.7);
  good_dX->GetYaxis()->SetTitle("Y [cm]");
  good_dX->GetYaxis()->SetTitleOffset(1.7);
  good_dX->GetZaxis()->SetTitle("Z [cm]");
  good_dX->GetZaxis()->SetTitleOffset(1.7);
  
  good_dX-> Draw("lego");
  c_1D_goodVoxels_dX.SaveAs("goodVoxels_dX.C");
  
  TCanvas c_1D_goodVoxels_dY("c_godVoxels_dY", "", 1300, 900);
  good_dY->GetXaxis()->SetTitle("X [cm]");
  good_dY->GetXaxis()->SetTitleOffset(1.7);
  good_dY->GetYaxis()->SetTitle("Y [cm]");
  good_dY->GetYaxis()->SetTitleOffset(1.7);
  good_dY->GetZaxis()->SetTitle("Z [cm]");
  good_dY->GetZaxis()->SetTitleOffset(1.7);
  c_1D_goodVoxels_dY.cd();
  good_dY-> Draw("lego");
  c_1D_goodVoxels_dY.SaveAs("goodVoxels_dY.C");
  
  TCanvas c_1D_goodVoxels_dZ("c_godVoxels_dZ", "", 1300, 900);
  good_dZ->GetXaxis()->SetTitle("X [cm]");
  good_dZ->GetXaxis()->SetTitleOffset(1.7);
  good_dZ->GetYaxis()->SetTitle("Y [cm]");
  good_dZ->GetYaxis()->SetTitleOffset(1.7);
  good_dZ->GetZaxis()->SetTitle("Z [cm]");
  good_dZ->GetZaxis()->SetTitleOffset(1.7);
  c_1D_goodVoxels_dZ.cd();
  good_dZ-> Draw("lego");
  c_1D_goodVoxels_dZ.SaveAs("goodVoxels_dZ.C");
  
    fileGoodVoxels->Write();
    fileGoodVoxels->Close();
    
  
}

void eFieldCalculator::compareCalibZXPlane(bool isData)
{
  
  const double xScale = TPC_X/Lx;
  const double yScale = TPC_Y/Ly;
  const double zScale = TPC_Z/Lz;
  const int   zRegions = 3; //top, middle, down
  int       zRegion = 0;
  const double zMaximum = 3.0;

  /*double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};*/
  int   zCuts[2] = {9, 17};
  /*TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);*/
  
  std::string inputLaser;
  std::string inputCosmic;
  if(isData){
     inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50_laserdata_v1098_bkwd.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_data_200k_Aug3_smoothed.root";
  }
  
  else{
     inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_MC_200k_Aug3_smoothed.root";
  }  

  TFile* fileLaser = new TFile(inputLaser.c_str());
  TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
  TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
  TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
  TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
  TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
  TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
  
  TFile* fileCosmic = new TFile(inputCosmic.c_str());
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
  


  //Convert from cosmic to laser coordiantes
  TH3F* diff_dX = new TH3F("diff_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dY = new TH3F("diff_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_dZ = new TH3F("diff_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
  TH3F* cosine_angles  = new TH3F("cosine_angles","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH3F* magnitude_cosmic  = new TH3F("magnitude_cosmic","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH3F* magnitude_laser  = new TH3F("magnitude_laser","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH3F* magnitude_diff  = new TH3F("magnitude_diff","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  bool plotX = false;
  bool plotY = false;
  bool plotZ = false;
  double driftV = 0.0;
  if(isData)
    driftV = 0.001;
  driftV = driftSign*driftV;

  for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
  {
    for(Int_t j = 1; j <= diff_dX->GetNbinsY(); j++)
    {
      for(Int_t k = 1; k <= diff_dX->GetNbinsZ(); k++)
      {
	plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
	
    std::vector<double> cosmic_Ds;
    std::vector<double> laser_Ds;
          
	double cosmic_x_correction = -i*driftV*cosmicDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
    double laser_x_correction  = i*driftV*laserDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
          
	/*if(j== 1 && k == 1)
	  std::cout << i*diff_dX->GetXaxis()->GetBinWidth(i) << std::endl;*/
	//x_correction = 0.0;
    
     cosmic_Ds.push_back(cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction);
     cosmic_Ds.push_back(cosmic_dY->GetBinContent(i,j,k));
     cosmic_Ds.push_back(cosmic_dZ->GetBinContent(i,j,k));
          
     laser_Ds.push_back(laser_dX->GetBinContent(i,j,k) + laser_x_correction);
     laser_Ds.push_back(laser_dY->GetBinContent(i,j,k));
     laser_Ds.push_back(laser_dZ->GetBinContent(i,j,k));
          
     double cosine    = getAngle(cosmic_Ds, laser_Ds);
     double laserMag  = getVectorMagnitude(laser_Ds);
     double cosmicMag = getVectorMagnitude(cosmic_Ds);
     double magDiff   = laserMag - cosmicMag;
     
     if(goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) ){
              magnitude_laser->SetBinContent(i,j,k,laserMag);
         
     }
          
     if(goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) )
             magnitude_cosmic->SetBinContent(i,j,k,cosmicMag);
          
			
	if(plotX)
	{
          diff_dX->SetBinContent(i,j,k,(laser_dX->GetBinContent(i,j,k) + laser_x_correction) -(cosmic_dX->GetBinContent(i,j,k)+cosmic_x_correction) );
	}

	if(plotY)
	{
          diff_dY->SetBinContent(i,j,k,laser_dY->GetBinContent(i,j,k)-cosmic_dY->GetBinContent(i,j,k));
	}

	if(plotZ)
	{
          diff_dZ->SetBinContent(i,j,k,laser_dZ->GetBinContent(i,j,k)-cosmic_dZ->GetBinContent(i,j,k));
	}
   
    if(plotX && plotY && plotZ){
      
      cosine_angles->SetBinContent(i,j,k,cosine);
      magnitude_diff->SetBinContent(i,j,k,magDiff);
       
          }
          
      }
    }
  }

  gStyle->SetTitleW(0.9);
  
  TH1F *h_diffXByRegion[zRegions];
  TH1F *h_diffYByRegion[zRegions];
  TH1F *h_diffZByRegion[zRegions];
  
  plotX = false;
  plotY = false;
  plotZ = false;
  
  for(int k = 0; k < zRegions; ++k){
      h_diffXByRegion[k] = new TH1F(Form("diff_1D_dX_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
      h_diffYByRegion[k] = new TH1F(Form("diff_1D_dY_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
      h_diffZByRegion[k] = new TH1F(Form("diff_1D_dZ_%d", k), "", ((int)zMaximum*100), -zMaximum, zMaximum);
  
  }

  for(int k = 1; k <= diff_dX->GetNbinsY(); k++)
  {
    TH2F laser_2D_dX(Form("laser_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F laser_2D_dY(Form("laser_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F laser_2D_dZ(Form("laser_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);

    TH2F cosmic_2D_dX(Form("cosmic_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F cosmic_2D_dY(Form("cosmic_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F cosmic_2D_dZ(Form("cosmic_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);

    TH2F diff_2D_dX(Form("diff_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F diff_2D_dY(Form("diff_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F diff_2D_dZ(Form("diff_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    
    TH2F cosine_2D(Form("cosine_2D_%d", k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F magnitude_cosmic_2D(Form("magnitude_cosmic_2D_%d", k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F magnitude_laser_2D(Form("magnitude_laser_2D_%d", k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
    TH2F magnitude_diff_2D(Form("magnitude_diff_2D_%d", k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
 
    if( k < zCuts[0])
      zRegion = 0;
    else if( k >= zCuts[0] && k < zCuts[1])
      zRegion = 1;
    else if( k >= zCuts[1])
      zRegion = 2;
    else
      std::cout << "Cannot determine zRegion" << std::endl; 
      
    //std::cout << k << " " << zRegion << std::endl;   
          

    for(Int_t i = 1; i <= diff_dX->GetNbinsX(); i++)
    {
      for(Int_t j = 1; j <= diff_dX->GetNbinsZ(); j++)
      {
        double cosmic_x_correction = -i*driftV*cosmicDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
        double laser_x_correction  =  i*driftV*laserDriftVScale*diff_dX->GetXaxis()->GetBinWidth(i);
          
	plotX = (goodLaser(laser_dX->GetBinContent(i,k,j), laser_dX_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dX->GetBinContent(i,k,j), cosmic_dX_err->GetBinContent(i,k,j)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,k,j), laser_dY_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dY->GetBinContent(i,k,j), cosmic_dY_err->GetBinContent(i,k,j)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,k,j), laser_dZ_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dZ->GetBinContent(i,k,j), cosmic_dZ_err->GetBinContent(i,k,j)) );
    
          if(goodLaser(laser_dX->GetBinContent(i,k,j), laser_dX_err->GetBinContent(i,k,j)) && goodLaser(laser_dY->GetBinContent(i,k,j), laser_dY_err->GetBinContent(i,k,j)) && goodLaser(laser_dZ->GetBinContent(i,k,j), laser_dZ_err->GetBinContent(i,k,j)) ){
              magnitude_laser_2D.SetBinContent(j,i, magnitude_laser->GetBinContent(i,k,j));
              if(k == 1 && i == 3)
                  magnitude_laser_2D.SetBinContent(j,i, -5.0);
          }
    if(goodCosmic(cosmic_dX->GetBinContent(i,k,j), cosmic_dX_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dY->GetBinContent(i,k,j), cosmic_dY_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dZ->GetBinContent(i,k,j), cosmic_dZ_err->GetBinContent(i,k,j)) )
              magnitude_cosmic_2D.SetBinContent(j,i, magnitude_cosmic->GetBinContent(i,k,j));
	
     laser_2D_dX.SetBinContent(j,i,laser_dX->GetBinContent(i,k,j)   + laser_x_correction);
     cosmic_2D_dX.SetBinContent(j,i,cosmic_dX->GetBinContent(i,k,j) + cosmic_x_correction);
	 
     if(plotX){
	  
	  diff_2D_dX.SetBinContent(j,i,diff_dX->GetBinContent(i,k,j));
	  h_diffXByRegion[zRegion]->Fill(diff_dX->GetBinContent(i,k,j));
	}
	
    laser_2D_dY.SetBinContent(j,i,laser_dY->GetBinContent(i,k,j));
    cosmic_2D_dY.SetBinContent(j,i,cosmic_dY->GetBinContent(i,k,j));
          
	if(plotY){
	  
      
	  diff_2D_dY.SetBinContent(j,i,diff_dY->GetBinContent(i,k,j));
	  h_diffYByRegion[zRegion]->Fill(diff_dY->GetBinContent(i,k,j));
	}
    
    laser_2D_dZ.SetBinContent(j,i,laser_dZ->GetBinContent(i,k,j));
    cosmic_2D_dZ.SetBinContent(j,i,cosmic_dZ->GetBinContent(i,k,j));
          
	if(plotZ){

	  diff_2D_dZ.SetBinContent(j,i,diff_dZ->GetBinContent(i,k,j));
	  h_diffZByRegion[zRegion]->Fill(diff_dZ->GetBinContent(i,k,j));
	
	
	}
          
    if(plotX && plotY && plotZ){
       cosine_2D.SetBinContent(j,i,cosine_angles->GetBinContent(i,k,j));
       magnitude_diff_2D.SetBinContent(j,i, magnitude_diff->GetBinContent(i,k,j));
    }
		  
      }
    }

    drawPlanarPlot(laser_2D_dX, k, "Laser #DeltaX", "laser_2D_dX", axisType::yAxis, 5.0);
    drawPlanarPlot(laser_2D_dY, k, "Laser #DeltaY", "laser_2D_dY", axisType::yAxis);
    drawPlanarPlot(laser_2D_dZ, k, "Laser #DeltaZ", "laser_2D_dZ", axisType::yAxis);
    drawPlanarPlot(cosmic_2D_dX, k, "Cosmic #DeltaX", "cosmic_2D_dX", axisType::yAxis, 5.0);
    drawPlanarPlot(cosmic_2D_dY, k, "Cosmic #DeltaY", "cosmic_2D_dY", axisType::yAxis);
    drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic #DeltaZ", "cosmic_2D_dZ", axisType::yAxis);
    drawPlanarPlot(diff_2D_dX, k, "Laser-Cosmic #DeltaX", "diff_2D_dX", axisType::yAxis, zMaximum);
    drawPlanarPlot(diff_2D_dY, k, "Laser-Cosmic #DeltaY", "diff_2D_dY", axisType::yAxis, zMaximum);
    drawPlanarPlot(diff_2D_dZ, k, "Laser-Cosmic #DeltaZ", "diff_2D_dZ", axisType::yAxis, zMaximum);
    drawPlanarPlot(cosine_2D,  k, "Cosine of Laser - Cosmic", "cos_2D", axisType::yAxis, 1.0);
    drawPlanarPlot(magnitude_cosmic_2D,  k, "Magnitude of Cosmic Deflection", "cosmic_2D_mag", axisType::yAxis, 10.0);
    drawPlanarPlot(magnitude_laser_2D,  k, "Magnitude of Laser Deflection", "laser_2D_mag", axisType::yAxis, 10.0);
    drawPlanarPlot(magnitude_diff_2D,  k, "Differnce in Magnitude", "diff_2D_mag", axisType::yAxis, 10.0);
    
  }
  double zOne = (((double) zCuts[0] - 1)/25.0)*TPC_Y;
  double zTwo = (((double) zCuts[1] - 1)/25.0)*TPC_Y;
  double maximum = -999.9;
  h_diffXByRegion[0]->Scale( TPC_Y / zOne);
  h_diffXByRegion[1]->Scale( TPC_Y / (zTwo - zOne));
  h_diffXByRegion[2]->Scale( TPC_Y / (TPC_Y - zTwo) );
  
  if(h_diffXByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[0]->GetMaximum();
  if(h_diffXByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[1]->GetMaximum();
  if(h_diffXByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffXByRegion[2]->GetMaximum();      
  
  h_diffXByRegion[0]->SetMaximum(1.1*maximum);
  h_diffXByRegion[1]->SetMaximum(1.1*maximum);
  h_diffXByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dX("c_1D_diff_dX", "", 600, 600);
  c_1D_diff_dX.cd();
  TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.55);
  legend->SetTextSize(0.03);
  h_diffXByRegion[0]->SetLineColor(kBlack);  
  legend->AddEntry(h_diffXByRegion[0], Form("Y < %.2f", zOne), "l");
  h_diffXByRegion[1]->SetLineColor(kRed+2);  
  legend->AddEntry(h_diffXByRegion[1], Form("%.2f #leq Y < %.2f", zOne, zTwo), "l");
  h_diffXByRegion[2]->SetLineColor(kBlue+2);  
  legend->AddEntry(h_diffXByRegion[2], Form("%.2f #geq Y", zTwo), "l");
  
  for(int i=0; i < zRegions; ++i){
    h_diffXByRegion[i]->SetTitle("Laser-Cosmic #DeltaX [cm]");
    h_diffXByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffXByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffXByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffXByRegion[i]->Draw("hist");
    else  
      h_diffXByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dX.SaveAs("1D_plots/ZX_diff_1D_dX.png");
  
  maximum = -999.9;
  h_diffYByRegion[0]->Scale( TPC_Y / zOne);
  h_diffYByRegion[1]->Scale( TPC_Y / (zTwo - zOne));
  h_diffYByRegion[2]->Scale( TPC_Y / (TPC_Y - zTwo) );
  
  if(h_diffYByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[0]->GetMaximum();
  if(h_diffYByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[1]->GetMaximum();
  if(h_diffYByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffYByRegion[2]->GetMaximum();      
  
  h_diffYByRegion[0]->SetMaximum(1.1*maximum);
  h_diffYByRegion[1]->SetMaximum(1.1*maximum);
  h_diffYByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dY("c_1D_diff_dY", "", 600, 600);
  c_1D_diff_dY.cd();
  h_diffYByRegion[0]->SetLineColor(kBlack);  
  h_diffYByRegion[1]->SetLineColor(kRed+2);  
  h_diffYByRegion[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffYByRegion[i]->SetTitle("Laser-Cosmic #DeltaY [cm]");
    h_diffYByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffYByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffYByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffYByRegion[i]->Draw("hist");
    else  
      h_diffYByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dY.SaveAs("1D_plots/ZX_diff_1D_dY.png");
  
  maximum = -999.9;
  h_diffZByRegion[0]->Scale( TPC_Y / zOne);
  h_diffZByRegion[1]->Scale( TPC_Y / (zTwo - zOne)  );
  h_diffZByRegion[2]->Scale( TPC_Y / (TPC_Y - zTwo) );
  
  if(h_diffZByRegion[0]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[0]->GetMaximum();
  if(h_diffZByRegion[1]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[1]->GetMaximum();
  if(h_diffZByRegion[2]->GetMaximum() > maximum)
     maximum = h_diffZByRegion[2]->GetMaximum();
   
  h_diffZByRegion[0]->SetMaximum(1.1*maximum);
  h_diffZByRegion[1]->SetMaximum(1.1*maximum);
  h_diffZByRegion[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diff_dZ("c_1D_diff_dZ", "", 600, 600);
  c_1D_diff_dZ.cd();
  h_diffZByRegion[0]->SetLineColor(kBlack);  
  h_diffZByRegion[1]->SetLineColor(kRed+2);  
  h_diffZByRegion[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffZByRegion[i]->SetTitle("Laser-Cosmic #DeltaZ [cm]");
    h_diffZByRegion[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffZByRegion[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffZByRegion[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffZByRegion[i]->Draw("hist");
    else  
      h_diffZByRegion[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diff_dZ.SaveAs("1D_plots/ZX_diff_1D_dZ.png");
  
}

void eFieldCalculator::compareTruth(bool sigmaDiff)
{
  
  const double xScale = TPC_X/Lx;
  const double yScale = TPC_Y/Ly;
  const double zScale = TPC_Z/Lz;
  const int   zRegions = 3; //upstream, middle, downstream
  int       zRegion = 0;
  double zMaximum = 5.0;

  double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};
  int   zCuts[2] = {21, 87};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TFile* fileLaser = new TFile("RecoCorr-N3-S50-LaserMC-2side-Anode.root");
  TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
  TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
  TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
  TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
  TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
  TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
  
    //Fix this at a later data
  TFile* fileCosmic = new TFile("output_hists_MC_200k_Aug3.root");
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
  
    /*
  TFile* fileCosmic = new TFile("MergedMaps.root");
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("combined_dX");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("combined_dY");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("combined_dZ");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("combined_dX_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("combined_dY_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("combined_dZ_Error");
    */
    
  
  TFile* fileTruth = new TFile("output_truth_hists.root");
  TH3F* truth_dX = (TH3F*) fileTruth->Get("True_Displacement_X");
  TH3F* truth_dY = (TH3F*) fileTruth->Get("True_Displacement_Y");
  TH3F* truth_dZ = (TH3F*) fileTruth->Get("True_Displacement_Z");
  TH3F* truth_dX_err = (TH3F*) fileTruth->Get("True_Displacement_X_Error");
  TH3F* truth_dY_err = (TH3F*) fileTruth->Get("True_Displacement_Y_Error");
  TH3F* truth_dZ_err = (TH3F*) fileTruth->Get("True_Displacement_Z_Error");    

  std::cout << xMin << " " << xMax << std::endl;
  std::cout << yMin << " " << yMax << std::endl;
  std::cout << zMin << " " << zMax << std::endl;

  
  TH3F* diff_cosmic_dX = new TH3F("diff_cosmic_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_cosmic_dY = new TH3F("diff_cosmic_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_cosmic_dZ = new TH3F("diff_cosmic_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
 
 TH3F* diff_laser_dX = new TH3F("diff_laser_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_laser_dY = new TH3F("diff_laser_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* diff_laser_dZ = new TH3F("diff_laser_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax); 
  bool plotX = false;
  bool plotY = false;
  bool plotZ = false;

  for(Int_t i = 1; i <= diff_laser_dX->GetNbinsX(); i++)
  {
    for(Int_t j = 1; j <= diff_laser_dX->GetNbinsY(); j++)
    {
      for(Int_t k = 1; k <= diff_laser_dX->GetNbinsZ(); k++)
      {
	plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
	
	
	double denom = 1.0;
	
	
	if(plotX)
	{
          (sigmaDiff) ? (denom = laser_dX_err->GetBinContent(i,j,k)) : (denom = 1.0);
	  //std::cout << denom << std::endl; 
	    
	  
	  diff_laser_dX->SetBinContent(i,j,k,(laser_dX->GetBinContent(i,j,k)-truth_dX->GetBinContent(i,j,k))/denom);
	  (sigmaDiff) ? (denom = cosmic_dX_err->GetBinContent(i,j,k)) : (denom = 1.0);	  	  
	  diff_cosmic_dX->SetBinContent(i,j,k,(cosmic_dX->GetBinContent(i,j,k)-truth_dX->GetBinContent(i,j,k))/denom);
	}

	if(plotY)
	{
	  (sigmaDiff) ? (denom = laser_dY_err->GetBinContent(i,j,k)) : (denom = 1.0);
	  diff_laser_dY->SetBinContent(i,j,k,(laser_dY->GetBinContent(i,j,k)-truth_dY->GetBinContent(i,j,k))/denom);
	  (sigmaDiff) ? (denom = cosmic_dY_err->GetBinContent(i,j,k)) : (denom = 1.0);
	  diff_cosmic_dY->SetBinContent(i,j,k,(cosmic_dY->GetBinContent(i,j,k)-truth_dY->GetBinContent(i,j,k))/denom);
	}

	if(plotZ)
	{
          (sigmaDiff) ? (denom = laser_dZ_err->GetBinContent(i,j,k)) : (denom = 1.0);
	  diff_laser_dZ->SetBinContent(i,j,k,(laser_dZ->GetBinContent(i,j,k)-truth_dZ->GetBinContent(i,j,k))/denom);
	  (sigmaDiff) ? (denom = cosmic_dZ_err->GetBinContent(i,j,k)) : (denom = 1.0);
	  diff_cosmic_dZ->SetBinContent(i,j,k,(cosmic_dZ->GetBinContent(i,j,k)-truth_dZ->GetBinContent(i,j,k))/denom);
	}
      }
    }
  }

  gStyle->SetTitleW(0.9);
  
  TH1F *h_diffXByRegionCosmic[zRegions];
  TH1F *h_diffYByRegionCosmic[zRegions];
  TH1F *h_diffZByRegionCosmic[zRegions];
  
  TH1F *h_diffXByRegionLaser[zRegions];
  TH1F *h_diffYByRegionLaser[zRegions];
  TH1F *h_diffZByRegionLaser[zRegions];
  
  plotX = false;
  plotY = false;
  plotZ = false;
  
  for(int k = 0; k < zRegions; ++k){
      h_diffXByRegionCosmic[k] = new TH1F(Form("Cdiff_1D_dX_%d", k), "", 300, -3.0, 3.0);
      h_diffYByRegionCosmic[k] = new TH1F(Form("Cdiff_1D_dY_%d", k), "", 300, -3.0, 3.0);
      h_diffZByRegionCosmic[k] = new TH1F(Form("Cdiff_1D_dZ_%d", k), "", 300, -3.0, 3.0);
      
      h_diffXByRegionLaser[k] = new TH1F(Form("Ldiff_1D_dX_%d", k), "", 300, -3.0, 3.0);
      h_diffYByRegionLaser[k] = new TH1F(Form("Ldiff_1D_dY_%d", k), "", 300, -3.0, 3.0);
      h_diffZByRegionLaser[k] = new TH1F(Form("Ldiff_1D_dZ_%d", k), "", 300, -3.0, 3.0);
  
  
  }

  for(Int_t k = 1; k <= diff_laser_dX->GetNbinsZ(); k++)
  {
    TH2F laser_2D_dX(Form("laser_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F laser_2D_dY(Form("laser_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F laser_2D_dZ(Form("laser_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);

    TH2F cosmic_2D_dX(Form("cosmic_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F cosmic_2D_dY(Form("cosmic_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
    TH2F cosmic_2D_dZ(Form("cosmic_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);


    if( k < zCuts[0])
      zRegion = 0;
    else if( k >= zCuts[0] && k < zCuts[1])
      zRegion = 1;
    else if( k >= zCuts[1])
      zRegion = 2;
    else
      std::cout << "Cannot determine zRegion" << std::endl;  
          

    for(Int_t i = 1; i <= diff_laser_dX->GetNbinsX(); i++)
    {
      for(Int_t j = 1; j <= diff_laser_dX->GetNbinsY(); j++)
      {
        
	plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
	
	
	
	if(plotX){	
	  laser_2D_dX.SetBinContent(i,j,diff_laser_dX->GetBinContent(i,j,k));
	  cosmic_2D_dX.SetBinContent(i,j,diff_cosmic_dX->GetBinContent(i,j,k));
	  h_diffXByRegionCosmic[zRegion]->Fill(diff_cosmic_dX->GetBinContent(i,j,k));
	  h_diffXByRegionLaser[zRegion]->Fill(diff_laser_dX->GetBinContent(i,j,k));
	}
	
	if(plotY){
	  
          laser_2D_dY.SetBinContent(i,j,diff_laser_dY->GetBinContent(i,j,k));
	  cosmic_2D_dY.SetBinContent(i,j,diff_cosmic_dY->GetBinContent(i,j,k));
	  h_diffYByRegionCosmic[zRegion]->Fill(diff_cosmic_dY->GetBinContent(i,j,k));
	  h_diffYByRegionLaser[zRegion]->Fill(diff_laser_dY->GetBinContent(i,j,k));
	}
	
	if(plotZ){
	  laser_2D_dZ.SetBinContent(i,j,diff_laser_dZ->GetBinContent(i,j,k));
	  cosmic_2D_dZ.SetBinContent(i,j,diff_cosmic_dZ->GetBinContent(i,j,k));
	  h_diffZByRegionCosmic[zRegion]->Fill(diff_cosmic_dZ->GetBinContent(i,j,k));
	  h_diffZByRegionLaser[zRegion]->Fill(diff_laser_dZ->GetBinContent(i,j,k));
	
	
	}  

      }
    }

    drawPlanarPlot(laser_2D_dX, k, "Laser - Truth dX", "LaserTruth_dX", axisType::zAxis, 2.5);
    drawPlanarPlot(laser_2D_dY, k, "Laser - Truth dY", "LaserTruth_dY", axisType::zAxis, 2.5);
    drawPlanarPlot(laser_2D_dZ, k, "Laser - Truth dZ", "LaserTruth_dZ", axisType::zAxis, 2.5);
      
    drawPlanarPlot(cosmic_2D_dX, k, "Cosmic - Truth dX", "CosmicTruth_dX", axisType::zAxis, 2.5);
    drawPlanarPlot(cosmic_2D_dY, k, "Cosmic - Truth dY", "CosmicTruth_dY", axisType::zAxis, 2.5);
    drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic - Truth dZ", "CosmicTruth_dZ", axisType::zAxis, 2.5);

    
  }
    
    for(Int_t k = 1; k <= diff_laser_dX->GetNbinsY(); k++)
    {
        TH2F laser_2D_dX(Form("laser_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F laser_2D_dY(Form("laser_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F laser_2D_dZ(Form("laser_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        TH2F cosmic_2D_dX(Form("cosmic_2D_dX_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F cosmic_2D_dY(Form("cosmic_2D_dY_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F cosmic_2D_dZ(Form("cosmic_2D_dZ_%d",k),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        /*
        if( k < zCuts[0])
            zRegion = 0;
        else if( k >= zCuts[0] && k < zCuts[1])
            zRegion = 1;
        else if( k >= zCuts[1])
            zRegion = 2;
        else
            std::cout << "Cannot determine zRegion" << std::endl;
        */
        
        for(Int_t i = 1; i <= diff_laser_dX->GetNbinsX(); i++)
        {
            for(Int_t j = 1; j <= diff_laser_dX->GetNbinsZ(); j++)
            {
                
                plotX = (goodLaser(laser_dX->GetBinContent(i,k,j), laser_dX_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dX->GetBinContent(i,k,j), cosmic_dX_err->GetBinContent(i,k,j)) );
                plotY = (goodLaser(laser_dY->GetBinContent(i,k,j), laser_dY_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dY->GetBinContent(i,k,j), cosmic_dY_err->GetBinContent(i,k,j)) );
                plotZ = (goodLaser(laser_dZ->GetBinContent(i,k,j), laser_dZ_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dZ->GetBinContent(i,k,j), cosmic_dZ_err->GetBinContent(i,k,j)) );
                
                
                
                if(plotX){
                    laser_2D_dX.SetBinContent(j,i,diff_laser_dX->GetBinContent(i,k,j));
                    cosmic_2D_dX.SetBinContent(j,i,diff_cosmic_dX->GetBinContent(i,k,j));
                   
                }
                
                if(plotY){
                    
                    laser_2D_dY.SetBinContent(j,i,diff_laser_dY->GetBinContent(i,k,j));
                    cosmic_2D_dY.SetBinContent(j,i,diff_cosmic_dY->GetBinContent(i,k,j));
                    
                }
                
                if(plotZ){
                    laser_2D_dZ.SetBinContent(j,i,diff_laser_dZ->GetBinContent(i,k,j));
                    cosmic_2D_dZ.SetBinContent(k,i,diff_cosmic_dZ->GetBinContent(i,k,j));
                    
                    
                }
                
            }
        }
        
        drawPlanarPlot(laser_2D_dX, k, "Laser - Truth dX", "LaserTruth_dX", axisType::yAxis, 2.5);
        drawPlanarPlot(laser_2D_dY, k, "Laser - Truth dY", "LaserTruth_dY", axisType::yAxis, 2.5);
        drawPlanarPlot(laser_2D_dZ, k, "Laser - Truth dZ", "LaserTruth_dZ", axisType::yAxis, 2.5);
        
        drawPlanarPlot(cosmic_2D_dX, k, "Cosmic - Truth dX", "CosmicTruth_dX", axisType::yAxis, 2.5);
        drawPlanarPlot(cosmic_2D_dY, k, "Cosmic - Truth dY", "CosmicTruth_dY", axisType::yAxis, 2.5);
        drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic - Truth dZ", "CosmicTruth_dZ", axisType::yAxis, 2.5);
        
        
    }
    
  double zOne = (((double) zCuts[0] - 1)/100.0)*TPC_Z;
  double zTwo = (((double) zCuts[1] - 1)/100.0)*TPC_Z;
  double maximum = -999.9;
  
  std::cout << calculateFWHM(h_diffXByRegionCosmic[0]) << std::endl;
  
  h_diffXByRegionCosmic[0]->Scale( TPC_Z / zOne);
  h_diffXByRegionCosmic[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffXByRegionCosmic[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffXByRegionCosmic[0]->GetMaximum() > maximum)
     maximum = h_diffXByRegionCosmic[0]->GetMaximum();
  if(h_diffXByRegionCosmic[1]->GetMaximum() > maximum)
     maximum = h_diffXByRegionCosmic[1]->GetMaximum();
  if(h_diffXByRegionCosmic[2]->GetMaximum() > maximum)
     maximum = h_diffXByRegionCosmic[2]->GetMaximum();      
  
  h_diffXByRegionCosmic[0]->SetMaximum(1.1*maximum);
  h_diffXByRegionCosmic[1]->SetMaximum(1.1*maximum);
  h_diffXByRegionCosmic[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffC_dX("c_1D_diffC_dX", "", 600, 600);
  c_1D_diffC_dX.cd();
  TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.55);
  legend->SetTextSize(0.03);
  h_diffXByRegionCosmic[0]->SetLineColor(kBlack);  
  legend->AddEntry(h_diffXByRegionCosmic[0], Form("Z < %.2f", zOne), "l");
  h_diffXByRegionCosmic[1]->SetLineColor(kRed+2);  
  legend->AddEntry(h_diffXByRegionCosmic[1], Form("%.2f #leq Z < %.2f", zOne, zTwo), "l");
  h_diffXByRegionCosmic[2]->SetLineColor(kBlue+2);  
  legend->AddEntry(h_diffXByRegionCosmic[2], Form("%.2f #leq Z", zTwo), "l");
  
  for(int i=0; i < zRegions; ++i){
    h_diffXByRegionCosmic[i]->SetTitle("Cosmic - Truth #DeltaX [cm]");
    h_diffXByRegionCosmic[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffXByRegionCosmic[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffXByRegionCosmic[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffXByRegionCosmic[i]->Draw("hist");
    else  
      h_diffXByRegionCosmic[i]->Draw("hist same");
  }
  legend->Draw("same");
  
  c_1D_diffC_dX.SaveAs("Truth_plots/diff_Cosmic_1D_dX.png");
  
  h_diffXByRegionLaser[0]->Scale( TPC_Z / zOne);
  h_diffXByRegionLaser[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffXByRegionLaser[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffXByRegionLaser[0]->GetMaximum() > maximum)
     maximum = h_diffXByRegionLaser[0]->GetMaximum();
  if(h_diffXByRegionLaser[1]->GetMaximum() > maximum)
     maximum = h_diffXByRegionLaser[1]->GetMaximum();
  if(h_diffXByRegionLaser[2]->GetMaximum() > maximum)
     maximum = h_diffXByRegionLaser[2]->GetMaximum();      
  
  h_diffXByRegionLaser[0]->SetMaximum(1.1*maximum);
  h_diffXByRegionLaser[1]->SetMaximum(1.1*maximum);
  h_diffXByRegionLaser[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffL_dX("c_1D_diffL_dX", "", 600, 600);
  c_1D_diffL_dX.cd();
  h_diffXByRegionLaser[0]->SetLineColor(kBlack);  
  h_diffXByRegionLaser[1]->SetLineColor(kRed+2);  
  h_diffXByRegionLaser[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffXByRegionLaser[i]->SetTitle("Laser - Truth #DeltaX [cm]");
    h_diffXByRegionLaser[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffXByRegionLaser[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffXByRegionLaser[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffXByRegionLaser[i]->Draw("hist");
    else  
      h_diffXByRegionLaser[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diffL_dX.SaveAs("Truth_plots/diff_Laser_1D_dX.png");
  
  maximum = -999.9;
  
  h_diffYByRegionCosmic[0]->Scale( TPC_Z / zOne);
  h_diffYByRegionCosmic[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffYByRegionCosmic[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffYByRegionCosmic[0]->GetMaximum() > maximum)
     maximum = h_diffYByRegionCosmic[0]->GetMaximum();
  if(h_diffYByRegionCosmic[1]->GetMaximum() > maximum)
     maximum = h_diffYByRegionCosmic[1]->GetMaximum();
  if(h_diffYByRegionCosmic[2]->GetMaximum() > maximum)
     maximum = h_diffYByRegionCosmic[2]->GetMaximum();      
  
  h_diffYByRegionCosmic[0]->SetMaximum(1.1*maximum);
  h_diffYByRegionCosmic[1]->SetMaximum(1.1*maximum);
  h_diffYByRegionCosmic[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffC_dY("c_1D_diffC_dY", "", 600, 600);
  c_1D_diffC_dY.cd();
  h_diffYByRegionCosmic[0]->SetLineColor(kBlack);  
  
  h_diffYByRegionCosmic[1]->SetLineColor(kRed+2);  
  h_diffYByRegionCosmic[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffYByRegionCosmic[i]->SetTitle("Cosmic - Truth #DeltaY [cm]");
    h_diffYByRegionCosmic[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffYByRegionCosmic[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffYByRegionCosmic[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffYByRegionCosmic[i]->Draw("hist");
    else  
      h_diffYByRegionCosmic[i]->Draw("hist same");
  }
  legend->Draw("same");
  
  c_1D_diffC_dY.SaveAs("Truth_plots/diff_Cosmic_1D_dY.png");
  
  h_diffYByRegionLaser[0]->Scale( TPC_Z / zOne);
  h_diffYByRegionLaser[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffYByRegionLaser[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffYByRegionLaser[0]->GetMaximum() > maximum)
     maximum = h_diffYByRegionLaser[0]->GetMaximum();
  if(h_diffYByRegionLaser[1]->GetMaximum() > maximum)
     maximum = h_diffYByRegionLaser[1]->GetMaximum();
  if(h_diffYByRegionLaser[2]->GetMaximum() > maximum)
     maximum = h_diffYByRegionLaser[2]->GetMaximum();      
  
  h_diffYByRegionLaser[0]->SetMaximum(1.1*maximum);
  h_diffYByRegionLaser[1]->SetMaximum(1.1*maximum);
  h_diffYByRegionLaser[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffL_dY("c_1D_diffL_dY", "", 600, 600);
  c_1D_diffL_dY.cd();
  h_diffYByRegionLaser[0]->SetLineColor(kBlack);  
  h_diffYByRegionLaser[1]->SetLineColor(kRed+2);  
  h_diffYByRegionLaser[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffYByRegionLaser[i]->SetTitle("Laser - Truth #DeltaY [cm]");
    h_diffYByRegionLaser[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffYByRegionLaser[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffYByRegionLaser[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffYByRegionLaser[i]->Draw("hist");
    else  
      h_diffYByRegionLaser[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diffL_dY.SaveAs("Truth_plots/diff_Laser_1D_dY.png");
  
  maximum = -999.9;
  
  h_diffZByRegionCosmic[0]->Scale( TPC_Z / zOne);
  h_diffZByRegionCosmic[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffZByRegionCosmic[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffZByRegionCosmic[0]->GetMaximum() > maximum)
     maximum = h_diffZByRegionCosmic[0]->GetMaximum();
  if(h_diffZByRegionCosmic[1]->GetMaximum() > maximum)
     maximum = h_diffZByRegionCosmic[1]->GetMaximum();
  if(h_diffZByRegionCosmic[2]->GetMaximum() > maximum)
     maximum = h_diffZByRegionCosmic[2]->GetMaximum();      
  
  h_diffZByRegionCosmic[0]->SetMaximum(1.1*maximum);
  h_diffZByRegionCosmic[1]->SetMaximum(1.1*maximum);
  h_diffZByRegionCosmic[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffC_dZ("c_1D_diffC_dZ", "", 600, 600);
  c_1D_diffC_dZ.cd();
  h_diffZByRegionCosmic[0]->SetLineColor(kBlack);    
  h_diffZByRegionCosmic[1]->SetLineColor(kRed+2);  
  h_diffZByRegionCosmic[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffZByRegionCosmic[i]->SetTitle("Cosmic - Truth #DeltaZ [cm]");
    h_diffZByRegionCosmic[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffZByRegionCosmic[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffZByRegionCosmic[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffZByRegionCosmic[i]->Draw("hist");
    else  
      h_diffZByRegionCosmic[i]->Draw("hist same");
  }
  legend->Draw("same");
  
  c_1D_diffC_dZ.SaveAs("Truth_plots/diff_Cosmic_1D_dZ.png");
  
  h_diffZByRegionLaser[0]->Scale( TPC_Z / zOne);
  h_diffZByRegionLaser[1]->Scale( TPC_Z / (zTwo - zOne));
  h_diffZByRegionLaser[2]->Scale( TPC_Z / (TPC_Z - zTwo) );
  
  if(h_diffZByRegionLaser[0]->GetMaximum() > maximum)
     maximum = h_diffZByRegionLaser[0]->GetMaximum();
  if(h_diffZByRegionLaser[1]->GetMaximum() > maximum)
     maximum = h_diffZByRegionLaser[1]->GetMaximum();
  if(h_diffZByRegionLaser[2]->GetMaximum() > maximum)
     maximum = h_diffZByRegionLaser[2]->GetMaximum();      
  
  h_diffZByRegionLaser[0]->SetMaximum(1.1*maximum);
  h_diffZByRegionLaser[1]->SetMaximum(1.1*maximum);
  h_diffZByRegionLaser[2]->SetMaximum(1.1*maximum);
  
  TCanvas c_1D_diffL_dZ("c_1D_diffL_dZ", "", 600, 600);
  c_1D_diffL_dZ.cd();
  h_diffZByRegionLaser[0]->SetLineColor(kBlack);  
  h_diffZByRegionLaser[1]->SetLineColor(kRed+2);  
  h_diffZByRegionLaser[2]->SetLineColor(kBlue+2);  
  
  for(int i=0; i < zRegions; ++i){
    h_diffZByRegionLaser[i]->SetTitle("Laser - Truth #DeltaZ [cm]");
    h_diffZByRegionLaser[i]->GetXaxis()->SetTitle("Difference (cm)");
    h_diffZByRegionLaser[i]->GetXaxis()->SetTitleOffset(1.0);
    h_diffZByRegionLaser[i]->GetXaxis()->SetTitleSize(0.04);
    if(i == 0)
      h_diffZByRegionLaser[i]->Draw("hist");
    else  
      h_diffZByRegionLaser[i]->Draw("hist same");
  }
  legend->Draw("same");
  c_1D_diffL_dZ.SaveAs("Truth_plots/diff_Laser_1D_dZ.png");
  
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
    //std::cout << inputVec[1] << " " << return_vec[1] << std::endl;
    return return_vec;
}


std::vector<double> eFieldCalculator::laserToCosmic(std::vector<double> inputVec){
    std::vector<double> return_vec;
    const double cathodeX = 2.5604;
    const double bottomY  = 1.1625;
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

double eFieldCalculator::calculateFWHM(TH1F *hist){
  const double epsilon = 5.0;
  int maxBin    = hist->GetMaximumBin();
  double maximum = hist->GetBinContent(maxBin); 
  int leftBin  = -1;
  int rightBin = -1;
  for(int i =1; i < maxBin; ++i){
     if( (hist->GetBinContent(i)/2.0 - maximum) < epsilon){
         leftBin = i;
	 break;
     }
        
  
  }
  
  for(int i = maxBin; i < hist->GetNbinsX(); ++i){
     if( (hist->GetBinContent(i)/2.0 - maximum) < epsilon){
         rightBin = i;
	 break;
     }
  
  }
  
  if(leftBin != -1 && rightBin != -1){
    return (hist->GetBinCenter(rightBin) - hist->GetBinCenter(leftBin) );
  
  }
  
  else
     return 0.0;
  


}

double eFieldCalculator::getDotProduct(std::vector<double> vecOne, std::vector<double> vecTwo){
  return (vecOne[0]*vecTwo[0] + vecOne[1]*vecTwo[1] + vecOne[2]*vecTwo[2]);

}

double eFieldCalculator::getVectorMagnitude(std::vector<double> vec){
    double magnitude = (pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2) );
    return sqrt(magnitude);

}

double eFieldCalculator::getAngle(std::vector<double> vecOne, std::vector<double> vecTwo){
   double dotProduct = getDotProduct(vecOne, vecTwo);
   double magOne = getVectorMagnitude(vecOne);
   double magTwo = getVectorMagnitude(vecTwo);
   
   return (dotProduct / (magOne*magTwo) );

}

std::vector<int> eFieldCalculator::calculateFWHMBins(TH1F *hist){
  const double epsilon = 25.0;
  int maxBin    = hist->GetMaximumBin();
  double maximum = hist->GetBinContent(maxBin); 
  int leftBin  = -1;
  int rightBin = -1;
  std::vector<int> return_vec;
  std::cout << maxBin << " " << maximum << std::endl;
  for(int i =1; i < maxBin; ++i){
     
     if( fabs(hist->GetBinContent(i) - (maximum/2.0) ) < epsilon){
         leftBin = i;
	 return_vec.push_back(leftBin);
	 break;
     }
        
  
  }
  
  for(int i = maxBin; i < hist->GetNbinsX(); ++i){
     if( fabs(hist->GetBinContent(i) - (maximum/2.0) ) < epsilon){
         
	 rightBin = i;
	 return_vec.push_back(rightBin);
	 break;
     }
  
  }
  
  return return_vec;
  


}

void eFieldCalculator::drawPlanarPlot(TH2F hist, int planeNum, const char* label, const char* filename, axisType axis, double zMax){
    int canWidth = 600;
    int canHeight = 600;
    const char *xAxisLabel;
    const char *yAxisLabel;
    const char *dir;
    char plotTitle[100];
    
    //XY Plane
    if(axis == axisType::zAxis){
      xAxisLabel = "X_{reco} [cm]";
      yAxisLabel = "Y_{reco} [cm]";
      dir        = "2DPlots_XY";
      sprintf(plotTitle, "%s [cm]:  ^{}Z_{reco} = %.2f cm", label,(((double) planeNum-1)/100.0)*TPC_Z);
    }
    
    //ZX Plane
    if(axis == axisType::yAxis){
      canWidth = 900;
      xAxisLabel = "Z_{reco} [cm]";
      yAxisLabel = "X_{reco} [cm]";
      dir        = "2DPlots_ZX";
      sprintf(plotTitle, "%s [cm]:  ^{}Y_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_Y) - TPC_Y/2.0) );
    }
    
    //YZ Plane
    if(axis == axisType::xAxis){
        canWidth = 900;
        xAxisLabel = "Z_{reco} [cm]";
        yAxisLabel = "y_{reco} [cm]";
        dir        = "2DPlots_YZ";
        sprintf(plotTitle, "%s [cm]:  ^{}X_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_X)));
    }
    
    
    TCanvas can(Form("can"),"",canWidth,canHeight);
    can.cd();
    hist.SetTitle(plotTitle);
    hist.GetXaxis()->SetTitle(xAxisLabel);
    hist.GetXaxis()->SetTitleOffset(1.0);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitle(yAxisLabel);
    hist.GetYaxis()->SetTitleOffset(1.2);
    hist.GetYaxis()->SetTitleSize(0.04);
    hist.GetZaxis()->SetNoExponent(kTRUE);
    hist.GetZaxis()->SetLabelSize(0.025);
    hist.SetStats(0);
    hist.Draw("COLZ");
    if(zMax != 0.0)
      hist.GetZaxis()->SetRangeUser(-zMax, zMax);
    can.SaveAs(Form("%s/%s_%d.png",dir,filename,planeNum));

}

void eFieldCalculator::drawPlanarPlot(TH2F *hist, int planeNum, const char* label, const char* filename, axisType axis, double zMax){
    int canWidth = 600;
    int canHeight = 600;
    const char *xAxisLabel;
    const char *yAxisLabel;
    const char *dir;
    char plotTitle[100];
    
    //XY Plane
    if(axis == axisType::zAxis){
        xAxisLabel = "X_{reco} [cm]";
        yAxisLabel = "Y_{reco} [cm]";
        dir        = "2DPlots_XY";
        if(planeNum > 0)
            sprintf(plotTitle, "%s [cm]:  ^{}Z_{reco} = %.2f cm", label,(((double) planeNum-1)/100.0)*TPC_Z);
        else
            sprintf(plotTitle, "%s", label);
        
    }
    
    //ZX Plane
    if(axis == axisType::yAxis){
        canWidth = 900;
        xAxisLabel = "Z_{reco} [cm]";
        yAxisLabel = "X_{reco} [cm]";
        dir        = "2DPlots_ZX";
        if(planeNum > 0)
            sprintf(plotTitle, "%s [cm]:  ^{}Y_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_Y) - TPC_Y/2.0) );
        else
            sprintf(plotTitle, "%s", label);
    }
    
    //YZ Plane
    if(axis == axisType::xAxis){
        canWidth = 900;
        xAxisLabel = "Z_{reco} [cm]";
        yAxisLabel = "y_{reco} [cm]";
        dir        = "2DPlots_YZ";
        if(planeNum > 0)
            sprintf(plotTitle, "%s [cm]:  ^{}X_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_X)));
        else
            sprintf(plotTitle, "%s", label);
    }
    
    
    TCanvas can(Form("can"),"",canWidth,canHeight);
    can.cd();
    hist->SetTitle(plotTitle);
    hist->GetXaxis()->SetTitle(xAxisLabel);
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitle(yAxisLabel);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetZaxis()->SetNoExponent(kTRUE);
    hist->GetZaxis()->SetLabelSize(0.025);
    hist->SetStats(0);
    hist->Draw("COLZ");
    if(zMax != 0.0)
        hist->GetZaxis()->SetRangeUser(-zMax, zMax);
    can.SaveAs(Form("%s/%s_%d.png",dir,filename,planeNum));
    
}

void eFieldCalculator::draw1DPlot(TH1F *histOne, TH1F *histTwo, int planeNum, const char* label, const char* filename, axisType axis, double zMax){
    int canWidth = 600;
    int canHeight = 600;
    const char *xAxisLabel;
    const char *yAxisLabel;
    const char *dir;
    char plotTitle[100];
    
    //XY Plane
    if(axis == axisType::zAxis){
        xAxisLabel = "#Delta #theta [degrees]";
        yAxisLabel = "Fraction of Tracks";
        dir        = "AnglePlots";
        sprintf(plotTitle, "%s ^{}Z_{reco} = %.2f cm", label,(((double) planeNum-1)/100.0)*TPC_Z);
    }
    
    //ZX Plane
    if(axis == axisType::yAxis){
        xAxisLabel = "#Delta #theta [degrees]";
        yAxisLabel = "Fraction of Tracks";
        dir        = "AnglePlots";
        sprintf(plotTitle, "%s ^{}Y_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_Y) - TPC_Y/2.0) );
    }
    
    //YZ Plane
    if(axis == axisType::xAxis){
        xAxisLabel = "#Delta #theta [degrees]";
        yAxisLabel = "Fraction of Tracks";
        dir        = "AnglePlots";
        sprintf(plotTitle, "%s  ^{}X_{reco} = %.2f cm", label, (double)(TPC_X*((25-planeNum)/25.0))/*(((((double) planeNum-1)/25.0)*TPC_X))*/);
    }
    
    //histOne->Scale(1.0/histOne->Integral());
    //histTwo->Scale(1.0/histTwo->Integral());
    
    histOne->SetLineWidth(2.0);
    histOne->SetLineColor(kRed);
    histTwo->SetLineWidth(2.0);
    histTwo->SetLineColor(kGreen+2);
    
    gStyle->SetTitleW(0.9);
    gStyle->SetOptStat(0);
    
    TCanvas *c_orig = new TCanvas(Form("can_orig"),"",canWidth,canHeight);
    TLegend *leg_combined = new TLegend(0.55,0.70,0.88,0.85);
    leg_combined->SetLineColor(kWhite);
    leg_combined->AddEntry(histOne,"Before Calibration","L");
    leg_combined->AddEntry(histTwo,"After Calibration","L");
    histOne->GetXaxis()->SetTitle(xAxisLabel);
    histOne->GetYaxis()->SetTitle(yAxisLabel);
    histOne->SetTitle(plotTitle);
    histOne->Draw("HIST");
    histTwo->Draw("HISTsame");
    leg_combined->Draw("same");
    histOne->Draw("AXISsame");
    histOne->SetMaximum(1.1*std::max(histOne->GetMaximum(),histTwo->GetMaximum()));
    histOne->SetMinimum(0.0001);
    c_orig->SaveAs(Form("%s_%d.png",filename, planeNum ));
    
}

void eFieldCalculator::compareFaces(bool isData){
    std::string inputLaser;
    std::string inputCosmic;
    std::vector<std::string> inputFace;
    if(isData){
        inputLaser  = "RecoCorr-N3-S50-Data-2side-Anode.root";
        inputCosmic = "output_hists_data_200k_Aug3.root";
        inputFace.push_back("data_top.root");
        inputFace.push_back("data_cathode.root");
    }
    
    else{
        inputLaser  = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
        inputCosmic = "output_hists_MC_200k_Aug3.root";
        inputFace.push_back("MC_top.root");
        inputFace.push_back("MC_cathode.root");
    }
    
    TFile *fileFacesTop = new TFile(inputFace[0].c_str());
    TFile *fileFacesCat = new TFile(inputFace[1].c_str());
    //TH3F* face_dX = (TH3F*) fileLaser->Get("th2_for_looking_at_reco_offset_values_top");
    TH2F* face_dY = (TH2F*) fileFacesTop->Get("th2_for_looking_at_reco_offset_values_top");
    TH2F* face_dX = (TH2F*) fileFacesCat->Get("th2_for_looking_at_reco_offset_values_cathode");
    //TH3F* face_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
    
    int nZbins  = face_dY->GetNbinsY();
    int nXbins  = face_dY->GetNbinsX();
    int nYbins  = 0;
    
    double driftV = 0.0;
    if(isData)
        driftV = 0.01;

    
    TFile* fileLaser = new TFile(inputLaser.c_str());
    TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
    TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
    TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
    TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
    TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
    TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
    
    TFile* fileCosmic = new TFile(inputCosmic.c_str());
    TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
    TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
    TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
    TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
    TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
    TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
    
    TH2F laser_2D_dY(Form("laser_2D_dY"),"",nZbins,zMin,zMax,nXbins,xMin,xMax);
    TH2F cosmic_2D_dY(Form("cosmic_2D_dY"),"",nZbins,zMin,zMax,nXbins,xMin,xMax);

    for(int i = 1; i <= nXbins; i++){
        for(int j = 1; j <= nZbins; j++){
            double xVal = face_dY->GetXaxis()->GetBinCenter(i);
            double zVal = face_dY->GetYaxis()->GetBinCenter(j);
            double laserVal  =  laser_dY->Interpolate(xVal,  105.0, zVal);
            double cosmicVal =  cosmic_dY->Interpolate(xVal, 105.0, zVal);
            double faceVal = -1.0*face_dY->GetBinContent(i,j);
            
            //bool  plotLaser  =  goodLaser(laser_dY->GetBinContent(i,25,j), laser_dY_err->GetBinContent(i,25,j));
            //bool  plotCosmic =  goodCosmic(cosmic_dY->GetBinContent(i,25,j), cosmic_dY_err->GetBinContent(i,25,j));
            bool plotLaser  = true;
            bool plotCosmic = true;
            /*
            if(plotLaser && plotCosmic)
              std::cout << laserVal << " " << cosmicVal << " " << fabs(face_dY->GetBinContent(i,j)) << std::endl;
            */
            if(plotLaser)
              laser_2D_dY.SetBinContent(j,i,  (laserVal  - faceVal ) );
            if(plotCosmic)
              cosmic_2D_dY.SetBinContent(j,i, (cosmicVal - faceVal ) );
            
            
            
        }
    }
    
    nZbins  = face_dX->GetNbinsX();
    nYbins  = face_dX->GetNbinsY();
    
    TH2F laser_2D_dX(Form("laser_2D_dX"),"",nZbins,zMin,zMax,nYbins,yMin,yMax);
    TH2F cosmic_2D_dX(Form("cosmic_2D_dX"),"",nZbins,zMin,zMax,nYbins,yMin,yMax);
    
    for(int i = 1; i <= nZbins; i++){
        for(int j = 1; j <= nYbins; j++){
            double zVal = face_dX->GetXaxis()->GetBinCenter(i);
            double yVal = face_dX->GetYaxis()->GetBinCenter(j);
            double xVal = 230.0;
            double x_correction = driftV*xVal;
            //double x_correction = 0.0;
            double laserVal  =   laser_dX->Interpolate(xVal, yVal, zVal);
            double cosmicVal =   cosmic_dX->Interpolate(xVal, yVal, zVal) - x_correction;
            double faceVal   =   -1.0*(face_dX->GetBinContent(i,j) + x_correction);
            
            //bool  plotLaser  =  goodLaser(laser_dY->GetBinContent(i,25,j), laser_dY_err->GetBinContent(i,25,j));
            //bool  plotCosmic =  goodCosmic(cosmic_dY->GetBinContent(i,25,j), cosmic_dY_err->GetBinContent(i,25,j));
            bool plotLaser  = true;
            bool plotCosmic = true;
            
            if(plotLaser && plotCosmic && j == 11)
                 std::cout << zVal << ", " << yVal << ": " <<  laserVal << " " << cosmicVal << " " <<  faceVal << std::endl;
             
            if(plotLaser)
               laser_2D_dX.SetBinContent(j,i,  (laserVal  - faceVal ) ) ;
            if(plotCosmic)
               cosmic_2D_dX.SetBinContent(j,i, (cosmicVal - faceVal ) ) ;
            
            
            
        }
    }
    
    drawPlanarPlot(laser_2D_dY, 26, "Laser - Face #DeltaY (top)", "laser_2D_top_dY", axisType::yAxis, 6.0);
    drawPlanarPlot(cosmic_2D_dY, 26, "Cosmic - Face #DeltaY (top)", "cosmic_2D_top_dY", axisType::yAxis, 6.0);
    
    drawPlanarPlot(laser_2D_dX, 26,  "Laser - Face #DeltaX (cathode)", "laser_2D_cathode_dX", axisType::xAxis, 4.0);
    drawPlanarPlot(cosmic_2D_dX, 26, "Cosmic - Face #DeltaX (cathode)", "cosmic_2D_cathode_dX", axisType::xAxis, 4.0);
}

void eFieldCalculator::doFits(){
    std::string inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
    std::string inputCosmic = "output_hists_MC_200k_Aug3.root";


   TFile* fileLaser = new TFile(inputLaser.c_str());
   TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
   TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
   TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
   TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
   TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
   TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");

   TFile* fileCosmic = new TFile(inputCosmic.c_str());
   TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
   TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
   TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
   TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
   TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
   TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
    
   TFile* fileTruth = new TFile("output_truth_hists.root");
   TH3F* truth_dX = (TH3F*) fileTruth->Get("True_Displacement_X");
   TH3F* truth_dY = (TH3F*) fileTruth->Get("True_Displacement_Y");
   TH3F* truth_dZ = (TH3F*) fileTruth->Get("True_Displacement_Z");
    
   TFile* fileWeights = new TFile("LaserCosmicWeights.root", "RECREATE");
   TH1F *weights_x = new TH1F("weights_x", "", 100, 0, 1.0);
   TH1F *chisq_x = new TH1F("chisq_x", "", 100, 0, 1.0);
   TH3F *weights_x_3D = new TH3F("weights_3D_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH1F *weights_y = new TH1F("weights_y", "", 100, 0, 1.0);
    TH1F *chisq_y = new TH1F("chisq_y", "", 100, 0, 1.0);
    TH3F *weights_y_3D = new TH3F("weights_3D_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);

    TH1F *weights_z = new TH1F("weights_z", "", 100, 0, 1.0);
    TH1F *chisq_z = new TH1F("chisq_z", "", 100, 0, 1.0);
    TH3F *weights_z_3D = new TH3F("weights_3D_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    //TF1 *f1 = new TF1("chi2", calcChi2, 0, 100.0, 2);
    
    
    for(int xBin = 1; xBin < cosmic_dX->GetNbinsX(); ++xBin){
        for(int yBin = 1; yBin < cosmic_dX->GetNbinsY(); ++yBin){
            for(int zBin = 1; zBin < cosmic_dX->GetNbinsZ(); ++zBin){
                bool useCosmic = goodCosmic(cosmic_dX->GetBinContent(xBin, yBin, zBin), cosmic_dX_err->GetBinContent(xBin, yBin, zBin));
                bool useLaser  = goodLaser(laser_dX->GetBinContent(xBin, yBin, zBin), laser_dX_err->GetBinContent(xBin, yBin, zBin));
                
                if(useLaser && useCosmic){
                
                  ROOT::Math::Minimizer* minX =  ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
                  minX->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
                  minX->SetMaxIterations(10000);  // for GSL
                  minX->SetTolerance(0.001);
                  minX->SetPrintLevel(kPrint);
                  ROOT::Math::Functor f1(&calcChi2,5);
                  minX->SetFunction(f1);
                
                  minX->SetVariable(0, "x", 0.5,  0.01);
                  minX->SetVariableLimits(0, 0.0, 1.0);
                
                  minX->SetFixedVariable(1, "cosmicVal", cosmic_dX->GetBinContent(xBin, yBin, zBin));
                  minX->SetFixedVariable(2, "laserVal", laser_dX->GetBinContent(xBin, yBin, zBin));
                  minX->SetFixedVariable(3, "truthVal", truth_dX->GetBinContent(xBin, yBin, zBin));
                  minX->SetFixedVariable(4, "cosmicErr",cosmic_dX_err->GetBinContent(xBin, yBin, zBin) );
                  minX->SetFixedVariable(5, "laserErr", laser_dX_err->GetBinContent(xBin, yBin, zBin) );
                  minX->Minimize();
                
                  const double *xs = minX->X();
                  weights_x->Fill(xs[0]);
                  chisq_x->Fill(minX->MinValue());
                  weights_x_3D->SetBinContent(xBin, yBin, zBin, xs[0]);
                }
                
                else{
                    weights_x->Fill(-1.0);
                    chisq_x->Fill(-1.0);
                    weights_x_3D->SetBinContent(xBin, yBin, zBin, -1.0);
                }
                
                useCosmic = goodCosmic(cosmic_dY->GetBinContent(xBin, yBin, zBin), cosmic_dY_err->GetBinContent(xBin, yBin, zBin));
                useLaser  = goodLaser(laser_dY->GetBinContent(xBin, yBin, zBin), laser_dY_err->GetBinContent(xBin, yBin, zBin));
                if(useLaser && useCosmic){
                    
                    ROOT::Math::Minimizer* minX =  ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
                    minX->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
                    minX->SetMaxIterations(10000);  // for GSL
                    minX->SetTolerance(0.001);
                    minX->SetPrintLevel(kPrint);
                    ROOT::Math::Functor f1(&calcChi2,5);
                    minX->SetFunction(f1);
                    
                    minX->SetVariable(0, "x", 0.5,  0.01);
                    minX->SetVariableLimits(0, 0.0, 1.0);
                    
                    minX->SetFixedVariable(1, "cosmicVal", cosmic_dY->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(2, "laserVal", laser_dY->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(3, "truthVal", truth_dY->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(4, "cosmicErr",cosmic_dY_err->GetBinContent(xBin, yBin, zBin) );
                    minX->SetFixedVariable(5, "laserErr", laser_dY_err->GetBinContent(xBin, yBin, zBin) );
                    minX->Minimize();
                    
                    const double *xs = minX->X();
                    weights_y->Fill(xs[0]);
                    chisq_y->Fill(minX->MinValue());
                    weights_y_3D->SetBinContent(xBin, yBin, zBin, xs[0]);
                }
                
                else{
                    weights_y->Fill(-1.0);
                    chisq_y->Fill(-1.0);
                    weights_y_3D->SetBinContent(xBin, yBin, zBin, -1.0);
                }
                
                useCosmic = goodCosmic(cosmic_dZ->GetBinContent(xBin, yBin, zBin), cosmic_dZ_err->GetBinContent(xBin, yBin, zBin));
                useLaser  = goodLaser(laser_dZ->GetBinContent(xBin, yBin, zBin), laser_dZ_err->GetBinContent(xBin, yBin, zBin));
                if(useLaser && useCosmic){
                    
                    ROOT::Math::Minimizer* minX =  ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
                    minX->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
                    minX->SetMaxIterations(10000);  // for GSL
                    minX->SetTolerance(0.001);
                    minX->SetPrintLevel(kPrint);
                    ROOT::Math::Functor f1(&calcChi2,5);
                    minX->SetFunction(f1);
                    
                    minX->SetVariable(0, "x", 0.5,  0.01);
                    minX->SetVariableLimits(0, 0.0, 1.0);
                    
                    minX->SetFixedVariable(1, "cosmicVal", cosmic_dZ->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(2, "laserVal", laser_dZ->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(3, "truthVal", truth_dZ->GetBinContent(xBin, yBin, zBin));
                    minX->SetFixedVariable(4, "cosmicErr",cosmic_dZ_err->GetBinContent(xBin, yBin, zBin) );
                    minX->SetFixedVariable(5, "laserErr", laser_dZ_err->GetBinContent(xBin, yBin, zBin) );
                    minX->Minimize();
                    
                    const double *xs = minX->X();
                    weights_z->Fill(xs[0]);
                    chisq_z->Fill(minX->MinValue());
                    weights_z_3D->SetBinContent(xBin, yBin, zBin, xs[0]);
                }
                
                else{
                    weights_z->Fill(-1.0);
                    chisq_z->Fill(-1.0);
                    weights_z_3D->SetBinContent(xBin, yBin, zBin, -1.0);
                }

            }
        }
    }
    
    for(int zBin = 1; zBin < weights_z_3D->GetNbinsZ(); ++zBin){
        TH2F weights_2D_dX(Form("weights_2D_dX_%d",zBin),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax) ;
        TH2F weights_2D_dY(Form("weights_2D_dY_%d",zBin),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax) ;
        TH2F weights_2D_dZ(Form("weights_2D_dZ_%d",zBin),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        for(int xBin = 1; xBin < weights_z_3D->GetNbinsX(); ++xBin){
          for(int yBin = 1; yBin < weights_z_3D->GetNbinsY(); ++yBin){
              weights_2D_dX.SetBinContent(xBin, yBin, weights_x_3D->GetBinContent(xBin, yBin, zBin) );
              weights_2D_dY.SetBinContent(xBin, yBin, weights_y_3D->GetBinContent(xBin, yBin, zBin) );
              weights_2D_dZ.SetBinContent(xBin, yBin, weights_z_3D->GetBinContent(xBin, yBin, zBin) );
            }
        }
         drawPlanarPlot(weights_2D_dX,  zBin, "Cosmic Weight dX", "weight_2D_dX", axisType::zAxis, 1.0);
         drawPlanarPlot(weights_2D_dY,  zBin, "Cosmic Weight dY", "weight_2D_dY", axisType::zAxis, 1.0);
         drawPlanarPlot(weights_2D_dZ,  zBin, "Cosmic Weight dZ", "weight_2D_dZ", axisType::zAxis, 1.0);
    }
    
    for(int yBin = 1; yBin < weights_z_3D->GetNbinsY(); ++yBin){
        TH2F weights_2D_dX(Form("weights_2D_dX_%d",yBin),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F weights_2D_dY(Form("weights_2D_dY_%d",yBin),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F weights_2D_dZ(Form("weights_2D_dZ_%d",yBin),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        for(int xBin = 1; xBin < weights_z_3D->GetNbinsX(); ++xBin){
            for(int zBin = 1; zBin < weights_z_3D->GetNbinsZ(); ++zBin){
                weights_2D_dX.SetBinContent(zBin, xBin, weights_x_3D->GetBinContent(xBin, yBin, zBin) );
                weights_2D_dY.SetBinContent(zBin, xBin, weights_y_3D->GetBinContent(xBin, yBin, zBin) );
                weights_2D_dZ.SetBinContent(zBin, xBin, weights_z_3D->GetBinContent(xBin, yBin, zBin) );
            }
        }
        drawPlanarPlot(weights_2D_dX,  yBin, "Cosmic Weight dX", "weight_2D_dX", axisType::yAxis, 1.0);
        drawPlanarPlot(weights_2D_dY,  yBin, "Cosmic Weight dY", "weight_2D_dY", axisType::yAxis, 1.0);
        drawPlanarPlot(weights_2D_dZ,  yBin, "Cosmic Weight dZ", "weight_2D_dZ", axisType::yAxis, 1.0);
    }
    
    fileWeights->Write();
    fileWeights->Close();
    
    
    
}

void eFieldCalculator::combineMaps(bool isData, bool skipLaser, bool skipCosmic){
  const int lowZ    = 0;
  const int highZ   = 1000;
  //const int lowZ    = 0;
  //const int highZ   = 1000;
  const int lowX    = 0;
  const int highX   = 100;
  //const int lowY    = 4;
  const int lowY    = 0;
  const int highY   = 210;
  double driftV = 0.0;
  
  
  if(isData)
      driftV = 0.007;
  
  if(skipLaser && skipCosmic){
      std::cout << "Skipping Laser and Cosmic, i.e. do nothing!" << std::endl;
      return;
  }
    
  std::string inputLaser;
  std::string inputCosmic;
  std::string inputGoodVoxels;
  std::string inputTruth;
  std::string outputName;
  if(skipLaser)
      outputName = "MergedMapsCosmicOnlySmoothed.root";
  else if(skipCosmic)
      outputName = "MergedMapsLaserOnly.root";
  else
      outputName = "MergedMapsCosmicAndLaserNoDriftV.root";
  
  if(isData){
     inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50_laserdata_v1098_bkwd.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_data_200k_Aug3_smoothed.root";
     inputGoodVoxels = "LaserCosmicWeights.root";
     inputTruth      = "output_truth_hists.root";
  }
  
  else{
     inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_MC_200k_Aug3_smoothed.root";
     inputGoodVoxels = "GoodVoxels.root";
     inputTruth      = "output_truth_hists.root";
  }
  
  TFile *outputHistos = new TFile(outputName.c_str(), "RECREATE");

  TFile* fileLaser = new TFile(inputLaser.c_str());
  TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
  TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
  TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
  TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
  TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
  TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
  
  TFile* fileCosmic = new TFile(inputCosmic.c_str());
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
  
  TFile* fileTruth = new TFile(inputTruth.c_str());
  TH3F* truth_dX = (TH3F*) fileTruth->Get("True_Displacement_X");
  TH3F* truth_dY = (TH3F*) fileTruth->Get("True_Displacement_Y");
  TH3F* truth_dZ = (TH3F*) fileTruth->Get("True_Displacement_Z");
  
  outputHistos->cd();
  
  TH3F* combine_dX = new TH3F("Reco_Displacement_X","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  TH3F* combine_dY = new TH3F("Reco_Displacement_Y","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  TH3F* combine_dZ = new TH3F("Reco_Displacement_Z","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  TH3F* combine_dX_err = new TH3F("Reco_Displacement_X_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  TH3F* combine_dY_err = new TH3F("Reco_Displacement_Y_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  TH3F* combine_dZ_err = new TH3F("Reco_Displacement_Z_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
  
  for(int i = 1; i <= combine_dX->GetNbinsX(); i++){

    for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
      for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
         
	 double numerator   = 0.0;
     double denominator = 0.0;
     double error       = 0.0;
	 
	 bool useCosmic = false;
	 bool useLaser  = false;
     bool goodAgreement = false;
     
     double x_correction = -i*driftV*cosmicDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i);
	 
	 useCosmic = (goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)  ) && !skipCosmic);
	 useLaser  = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && i < highX && i >= lowX && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
     //std:: cout << "Use Cosmic: " << useCosmic << " " << "Use Laser: " << useLaser << std::endl;
    if(useCosmic && useLaser){
       numerator = ( (cosmic_dX->GetBinContent(i,j,k) +x_correction)/cosmic_dX_err->GetBinContent(i,j,k) + (laser_dX->GetBinContent(i,j,k))/laser_dX_err->GetBinContent(i,j,k) );
	   denominator = 1/cosmic_dX_err->GetBinContent(i,j,k) + 1/laser_dX_err->GetBinContent(i,j,k);
	   error = sqrt(pow(cosmic_dX_err->GetBinContent(i,j,k), 2) + pow(laser_dX_err->GetBinContent(i,j,k), 2));
        
         }
	 
	 else if(!useCosmic && useLaser){
       numerator = laser_dX->GetBinContent(i,j,k);
	   denominator = 1.0;
	   error = laser_dX_err->GetBinContent(i,j,k);
	 }
	 
	 else if(useCosmic && !useLaser){
       numerator = cosmic_dX->GetBinContent(i,j,k) + x_correction;
	   denominator = 1.0;
	   error = cosmic_dX_err->GetBinContent(i,j,k);
	 }
	 
	 else{
       //std::cout << laser_dX->GetBinContent(i,j,k) << " " << laser_dX_err->GetBinContent(i,j,k) << std::endl;
       numerator = 0.0;
	   denominator = 1.0;
	   error = 0.0;
	 }
	 
	 if(denominator != 0){
	    combine_dX->SetBinContent(i,j,k, numerator/denominator);
	    combine_dX_err->SetBinContent(i,j,k, error);
	 
	 }
	 
	 else
	   std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
	 
	 useCosmic = (goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) && !skipCosmic);
	 useLaser  = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
	 
	if(useCosmic && useLaser){
           numerator = (cosmic_dY->GetBinContent(i,j,k)/cosmic_dY_err->GetBinContent(i,j,k) + laser_dY->GetBinContent(i,j,k)/laser_dY_err->GetBinContent(i,j,k) );
	   denominator = 1/cosmic_dY_err->GetBinContent(i,j,k) + 1/laser_dY_err->GetBinContent(i,j,k);
	   error = sqrt(pow(cosmic_dY_err->GetBinContent(i,j,k), 2) + pow(laser_dY_err->GetBinContent(i,j,k), 2));
         }
	 
	 else if(!useCosmic && useLaser){
           numerator = laser_dY->GetBinContent(i,j,k);
	   denominator = 1.0;
	   error = laser_dY_err->GetBinContent(i,j,k);
	 }
	 
	 else if(useCosmic && !useLaser){
           numerator = cosmic_dY->GetBinContent(i,j,k);
	   denominator = 1.0;
	   error = cosmic_dY_err->GetBinContent(i,j,k);
	 }
	 
	 else{
       //std::cout << laser_dY->GetBinContent(i,j,k) << " " << laser_dY_err->GetBinContent(i,j,k) << std::endl;
       numerator = 0.0;
	   denominator = 1.0;
	   error = 0.0;
	 }
	 
	 if(denominator != 0){
	    combine_dY->SetBinContent(i,j,k, numerator/denominator);
	    combine_dY_err->SetBinContent(i,j,k, error);
	 
	 }
	 
	 else
	   std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
	 
	 useCosmic = (goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) && !skipCosmic);
	 useLaser  = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
	 
	 if(useCosmic && useLaser){
           numerator = (cosmic_dZ->GetBinContent(i,j,k)/cosmic_dZ_err->GetBinContent(i,j,k) + laser_dZ->GetBinContent(i,j,k)/laser_dZ_err->GetBinContent(i,j,k) );
	   denominator = 1/cosmic_dZ_err->GetBinContent(i,j,k) + 1/laser_dZ_err->GetBinContent(i,j,k);
	   error = sqrt(pow(cosmic_dZ_err->GetBinContent(i,j,k), 2) + pow(laser_dZ_err->GetBinContent(i,j,k), 2));
         }
	 
	 else if(!useCosmic && useLaser){
           numerator = laser_dZ->GetBinContent(i,j,k);
	   denominator = 1.0;
	   error = laser_dZ_err->GetBinContent(i,j,k);
	 }
	 
	 else if(useCosmic && !useLaser){
           numerator = cosmic_dZ->GetBinContent(i,j,k);
	   denominator = 1.0;
	   error = cosmic_dZ_err->GetBinContent(i,j,k);
	 }
	 
	 else{
       //std::cout << laser_dZ->GetBinContent(i,j,k) << " " << laser_dZ_err->GetBinContent(i,j,k) << std::endl;
       numerator = 0.0;
	   denominator = 1.0;
	   error = 0.0;
	 }
	 
	 if(denominator != 0){
	    combine_dZ->SetBinContent(i,j,k, numerator/denominator);
	    combine_dZ_err->SetBinContent(i,j,k, error);
	 
	 }
	 
	 else
	   std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;    
	 
	    
      
      
      }//end of loop over y bins
    }//end of loop over x bins
  
  }//end of loop over z bins
    
    
  
  
    
  for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
      
      for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
          for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
              combine_2D_dX.SetBinContent(i,j,combine_dX->GetBinContent(i,j,k));
              combine_2D_dY.SetBinContent(i,j,combine_dY->GetBinContent(i,j,k));
              combine_2D_dZ.SetBinContent(i,j,combine_dZ->GetBinContent(i,j,k));
    
          }
      }
      drawPlanarPlot(combine_2D_dX, k, "Combined #DeltaX", "combined_2D_dX", axisType::zAxis);
      drawPlanarPlot(combine_2D_dY, k, "Combined #DeltaY", "combined_2D_dY", axisType::zAxis);
      drawPlanarPlot(combine_2D_dZ, k, "Combined #DeltaZ", "combined_2D_dZ", axisType::zAxis);
  }
    
    outputHistos->Write();
    outputHistos->Close();

}

void eFieldCalculator::combineWeightedMaps(){
    const int lowZ    = 21;
    const int highZ   = 87;
    const int lowX    = 0;
    const int highX   = 40;
    const int lowY    = 4;
    const int highY   = 210;
    double driftV = 0.01;

    std::string outputName = "MergedMapsWeighted.root";
    std::string  inputLaser = "RecoCorr-N3-S50-Data-2side-Anode.root";
    std::string  inputCosmic = "output_hists_data_200k_Aug3.root";
    std::string  inputGoodVoxels = "LaserCosmicWeights.root";
    std::string  inputTruth      = "output_truth_hists.root";
    
    
    TFile *outputHistos = new TFile(outputName.c_str(), "RECREATE");
    
    TFile* fileLaser = new TFile(inputLaser.c_str());
    TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
    TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
    TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
    TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
    TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
    TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
    
    TFile* fileCosmic = new TFile(inputCosmic.c_str());
    TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
    TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
    TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
    TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
    TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
    TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
    
    TFile* fileTruth = new TFile(inputTruth.c_str());
    TH3F* truth_dX = (TH3F*) fileTruth->Get("True_Displacement_X");
    TH3F* truth_dY = (TH3F*) fileTruth->Get("True_Displacement_Y");
    TH3F* truth_dZ = (TH3F*) fileTruth->Get("True_Displacement_Z");
    
    TFile* fileGoodVoxels = new TFile(inputGoodVoxels.c_str());
    TH3F*  weights_dX = (TH3F*) fileGoodVoxels->Get("weights_3D_dX");
    TH3F*  weights_dY = (TH3F*) fileGoodVoxels->Get("weights_3D_dY");
    TH3F*  weights_dZ = (TH3F*) fileGoodVoxels->Get("weights_3D_dZ");
    
    outputHistos->cd();
    
    TH3F* combine_dX = new TH3F("combined_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dY = new TH3F("combined_dY","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dZ = new TH3F("combined_dZ","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH3F* combine_dX_err = new TH3F("combined_dX_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dY_err = new TH3F("combined_dY_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dZ_err = new TH3F("combined_dZ_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH1F *chi2_plot = new TH1F("chi2_plot", "Chi2", 1000, 0, 100.0);
    
    for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
        
        for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                
                double numerator   = 0.0;
                double denominator = 0.0;
                double error       = 0.0;
                
                bool useCosmic = false;
                bool useLaser  = false;
                bool goodAgreement = false;
                
                double x_correction = -i*driftV*cosmic_dX->GetXaxis()->GetBinWidth(i);
                double cosmicWeight = weights_dX->GetBinContent(i,j,k);
                double laserWeight  = 1.0 - cosmicWeight;
                
                useCosmic = goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k));
                useLaser  = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY);
                //useLaser = false;
                if(useCosmic && useLaser){
                    numerator = ( cosmicWeight*(cosmic_dX->GetBinContent(i,j,k)+x_correction)  + laserWeight*laser_dX->GetBinContent(i,j,k) );
                    denominator = 1.0;
                    error = sqrt(pow(cosmic_dX_err->GetBinContent(i,j,k), 2) + pow(laser_dX_err->GetBinContent(i,j,k), 2));
                    
                }
                
                else if(!useCosmic && useLaser){
                    numerator = laser_dX->GetBinContent(i,j,k);
                    denominator = 1.0;
                    error = laser_dX_err->GetBinContent(i,j,k);
                }
                
                else if(useCosmic && !useLaser){
                    numerator = cosmic_dX->GetBinContent(i,j,k) + x_correction;
                    denominator = 1.0;
                    error = cosmic_dX_err->GetBinContent(i,j,k);
                }
                
                else{
                    numerator = 0.0;
                    denominator = 1.0;
                    error = 0.0;
                }
                
                if(denominator != 0){
                    combine_dX->SetBinContent(i,j,k, numerator/denominator);
                    combine_dX_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
                
                useCosmic = goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k) );
                useLaser  = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY);
                 cosmicWeight = weights_dY->GetBinContent(i,j,k);
                laserWeight  = 1.0 - cosmicWeight;
                
                if(useCosmic && useLaser){
                    numerator = (cosmicWeight*cosmic_dY->GetBinContent(i,j,k) + laserWeight*laser_dY->GetBinContent(i,j,k) );
                    denominator = 1.0;
                    error = sqrt(pow(cosmic_dY_err->GetBinContent(i,j,k), 2) + pow(laser_dY_err->GetBinContent(i,j,k), 2));
                }
                
                else if(!useCosmic && useLaser){
                    numerator = laser_dY->GetBinContent(i,j,k);
                    denominator = 1.0;
                    error = laser_dY_err->GetBinContent(i,j,k);
                }
                
                else if(useCosmic && !useLaser){
                    numerator = cosmic_dY->GetBinContent(i,j,k);
                    denominator = 1.0;
                    error = cosmic_dY_err->GetBinContent(i,j,k);
                }
                
                else{
                    numerator = 0.0;
                    denominator = 1.0;
                    error = 0.0;
                }
                
                if(denominator != 0){
                    combine_dY->SetBinContent(i,j,k, numerator/denominator);
                    combine_dY_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
                
                useCosmic = goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k));
                useLaser  = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY);
                cosmicWeight = weights_dZ->GetBinContent(i,j,k);
                laserWeight  = 1.0 - cosmicWeight;
                
                if(useCosmic && useLaser){
                    numerator = (cosmicWeight*cosmic_dZ->GetBinContent(i,j,k) + laserWeight*laser_dZ->GetBinContent(i,j,k) );
                    denominator = 1.0;
                    error = sqrt(pow(cosmic_dZ_err->GetBinContent(i,j,k), 2) + pow(laser_dZ_err->GetBinContent(i,j,k), 2));
                }
                
                else if(!useCosmic && useLaser){
                    numerator = laser_dZ->GetBinContent(i,j,k);
                    denominator = 1.0;
                    error = laser_dZ_err->GetBinContent(i,j,k);
                }
                
                else if(useCosmic && !useLaser){
                    numerator = cosmic_dZ->GetBinContent(i,j,k);
                    denominator = 1.0;
                    error = cosmic_dZ_err->GetBinContent(i,j,k);
                }
                
                else{
                    numerator = 0.0;
                    denominator = 1.0;
                    error = 0.0;
                }
                
                if(denominator != 0){
                    combine_dZ->SetBinContent(i,j,k, numerator/denominator);
                    combine_dZ_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
                
                
                
                
            }//end of loop over y bins
        }//end of loop over x bins
        
    }//end of loop over z bins
    
    
    
    
    
    for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
            for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
                combine_2D_dX.SetBinContent(i,j,combine_dX->GetBinContent(i,j,k));
                combine_2D_dY.SetBinContent(i,j,combine_dY->GetBinContent(i,j,k));
                combine_2D_dZ.SetBinContent(i,j,combine_dZ->GetBinContent(i,j,k));
                
            }
        }
        drawPlanarPlot(combine_2D_dX, k, "Combined #DeltaX", "combined_2D_dX", axisType::zAxis);
        drawPlanarPlot(combine_2D_dY, k, "Combined #DeltaY", "combined_2D_dY", axisType::zAxis);
        drawPlanarPlot(combine_2D_dZ, k, "Combined #DeltaZ", "combined_2D_dZ", axisType::zAxis);
    }
    
    outputHistos->Write();
    outputHistos->Close();
    
}

void eFieldCalculator::combineMaps(int lowX, int highX, int lowY, int highY, int lowZ, int highZ, bool isData){
    /*const int lowZ    = 21;
    const int highZ   = 87;
    const int lowX    = 0;
    const int highX   = 100;
    const int lowY    = 0;
    const int highY   = 210;*/
    double driftV = 0.0;
    const bool SmoothBoundary = true;
    
    
    //if(isData)
    //    driftV = 0.008;
    //For tuning drift V scale
    if(isData)
        driftV = 0.001;
    driftV = driftSign*driftV;
    std::cout << "Laser Scale: " << laserDriftVScale << " Cosmic scale: " <<cosmicDriftVScale << std::endl;
    std::string inputLaser;
    std::string inputCosmic;
    std::string inputGoodVoxels;
    std::string inputTruth;
    std::string outputName;
    if(isData)
        outputName = "MergedMapsSmoothCosmicAndLaserNoDriftVSkip.root";
    else
        outputName = "MergedMapsSmoothCosmicAndLaserMC.root";
    
    if(isData){
        inputLaser = "/uboone/data/users/joelam/SCEInputFiles/RecoCorr-N3-S50_laserdata_v1098_bkwd.root";
        inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_data_200k_Aug3_smoothed.root";
        inputGoodVoxels = "LaserCosmicWeights.root";
        inputTruth      = "output_truth_hists.root";
    }
    
    else{
        inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
        inputCosmic = "/uboone/data/users/joelam/SCEInputFiles/output_hists_MC_200k_Aug3_smoothed.root";
        inputGoodVoxels = "GoodVoxels.root";
        inputTruth      = "output_truth_hists.root";
    }
    
    TFile *outputHistos = new TFile(outputName.c_str(), "RECREATE");
    
    TFile* fileLaser = new TFile(inputLaser.c_str());
    TH3F* laser_dX = (TH3F*) fileLaser->Get("Reco_Displacement_X");
    TH3F* laser_dY = (TH3F*) fileLaser->Get("Reco_Displacement_Y");
    TH3F* laser_dZ = (TH3F*) fileLaser->Get("Reco_Displacement_Z");
    TH3F* laser_dX_err = (TH3F*) fileLaser->Get("Reco_Displacement_X_Error");
    TH3F* laser_dY_err = (TH3F*) fileLaser->Get("Reco_Displacement_Y_Error");
    TH3F* laser_dZ_err = (TH3F*) fileLaser->Get("Reco_Displacement_Z_Error");
    
    TFile* fileCosmic = new TFile(inputCosmic.c_str());
    TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
    TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
    TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
    TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
    TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
    TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");
    
    TFile* fileTruth = new TFile(inputTruth.c_str());
    TH3F* truth_dX = (TH3F*) fileTruth->Get("True_Displacement_X");
    TH3F* truth_dY = (TH3F*) fileTruth->Get("True_Displacement_Y");
    TH3F* truth_dZ = (TH3F*) fileTruth->Get("True_Displacement_Z");
    
    outputHistos->cd();
    
    TH3F* combine_dX = new TH3F("Reco_Displacement_X","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dY = new TH3F("Reco_Displacement_Y","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dZ = new TH3F("Reco_Displacement_Z","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    TH3F* combine_dX_err = new TH3F("Reco_Displacement_X_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dY_err = new TH3F("Reco_Displacement_Y_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    TH3F* combine_dZ_err = new TH3F("Reco_Displacement_Z_Error","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
        
        for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                
                double numerator   = 0.0;
                double denominator = 0.0;
                double error       = 0.0;
                
                bool useCosmic = false;
                bool useLaser  = false;
                bool goodAgreement = false;
                bool atXLowerBoundary = false;
                bool atYLowerBoundary = false;
                bool atZLowerBoundary = false;
                bool atXUpperBoundary = false;
                bool atYUpperBoundary = false;
                bool atZUpperBoundary = false;
                double cosmic_x_correction = -i*driftV*cosmicDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i);
                double laser_x_correction  = i*driftV*laserDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i);
                
                /*
                if(driftSign > 0.0)
                    laser_x_correction = 0.0;
                else
                    cosmic_x_correction = 0.0;
                */
                //std::cout << "Cosmic: " << cosmic_x_correction << " Laser: " << laser_x_correction << std::endl;
                useCosmic = goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k) );
                useLaser  = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && i < highX && i >= lowX && k >= lowZ && k < highZ && j >= lowY && j < highY);
                atXLowerBoundary = (i == lowX && j > lowY && j < highY && k > lowZ && k < highZ);
                atYLowerBoundary = (j == lowY && i > lowX && i < highX && k > lowZ && k < highZ);
                atZLowerBoundary = (k == lowZ && i > lowX && i < highX && j > lowY && j < highY);
                
                atXUpperBoundary = (i == highX && j > lowY && j < highY && k > lowZ && k < highZ);
                atYUpperBoundary = (j == highY && i > lowX && i < highX && k > lowZ && k < highZ);
                atZUpperBoundary = (k == highZ && i > lowX && i < highX && j > lowY && j < highY);

                if(SmoothBoundary && (atXLowerBoundary || atXUpperBoundary) ){
                    numerator = ( cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction + laser_dX->GetBinContent(i,j,k) + laser_x_correction + cosmic_dX->GetBinContent(i+1,j,k) -(i+1)*driftV*cosmicDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i) +  laser_dX->GetBinContent(i+1,j,k)  - (i+1)*driftV*laserDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i) + cosmic_dX->GetBinContent(i-1,j,k) - (i-1)*driftV*cosmicDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i) + laser_dX->GetBinContent(i-1,j,k) - (i-1)*driftV*laserDriftVScale*cosmic_dX->GetXaxis()->GetBinWidth(i) );
                    denominator = 6.0;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(SmoothBoundary && (atYLowerBoundary || atYUpperBoundary) ){
                    numerator = ( cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction + laser_dX->GetBinContent(i,j,k) + laser_x_correction + cosmic_dX->GetBinContent(i,j+1,k) + cosmic_x_correction + laser_dX->GetBinContent(i,j+1,k) + laser_x_correction + cosmic_dX->GetBinContent(i,j-1,k) + cosmic_x_correction + laser_dX->GetBinContent(i,j-1,k) + laser_x_correction );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(SmoothBoundary && (atZLowerBoundary || atZUpperBoundary) ){
                    numerator = ( cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction + laser_dX->GetBinContent(i,j,k) + laser_x_correction + cosmic_dX->GetBinContent(i,j,k+1) + cosmic_x_correction + laser_dX->GetBinContent(i,j,k+1) + laser_x_correction + cosmic_dX->GetBinContent(i,j,k-1) + cosmic_x_correction + laser_dX->GetBinContent(i,j,k-1) + laser_x_correction );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                
                else if(useLaser){
                    numerator = (laser_dX->GetBinContent(i,j,k) + laser_x_correction);
                    denominator = 1.0;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                else if(useCosmic){
                    numerator = ( (cosmic_dX->GetBinContent(i,j,k) + cosmic_x_correction)  );
                    denominator = 1.0;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                //use the MC.
                else{
                    //std::cout << laser_dX->GetBinContent(i,j,k) << " " << laser_dX_err->GetBinContent(i,j,k) << std::endl;
                    std::cout << "Using MC! " << i << " " << " " << j << " " << k << std::endl;
                    numerator = truth_dX->GetBinContent(i,j,k) + laser_x_correction;
                    denominator = 1.0;
                    //in this case, error is ginormous
                    error = 0.0;
                }
                
                if(denominator != 0){
                    combine_dX->SetBinContent(i,j,k, numerator/denominator);
                    combine_dX_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
                
                useCosmic = goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k) );
                useLaser  = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && i < highX && i >= lowX && k >= lowZ && k < highZ && j >= lowY && j < highY);
                
                if(SmoothBoundary && (atXLowerBoundary || atXUpperBoundary) ){
                    numerator = ( cosmic_dY->GetBinContent(i,j,k) + laser_dY->GetBinContent(i,j,k) + cosmic_dY->GetBinContent(i+1,j,k) + laser_dY->GetBinContent(i+1,j,k) + cosmic_dY->GetBinContent(i-1,j,k) + laser_dY->GetBinContent(i-1,j,k) );
                    denominator = 6.0;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                    
                    
                    
                }
                
                else if(SmoothBoundary && (atYLowerBoundary || atYUpperBoundary) ){
                    numerator = ( cosmic_dY->GetBinContent(i,j,k) + laser_dY->GetBinContent(i,j,k) + cosmic_dY->GetBinContent(i,j+1,k) + laser_dY->GetBinContent(i,j+1,k) + cosmic_dY->GetBinContent(i,j-1,k) + laser_dY->GetBinContent(i,j-1,k) );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(SmoothBoundary && (atZLowerBoundary || atZUpperBoundary) ){
                    numerator = ( cosmic_dY->GetBinContent(i,j,k) + laser_dY->GetBinContent(i,j,k) + cosmic_dY->GetBinContent(i,j,k+1) + laser_dY->GetBinContent(i,j,k+1) + cosmic_dY->GetBinContent(i,j,k-1) + laser_dY->GetBinContent(i,j,k-1) );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                    
                    
                    
                }
                
                else if(useLaser){
                    numerator = laser_dY->GetBinContent(i,j,k);
                    denominator = 1.0;
                    if(!isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                else if(useCosmic){
                    numerator = ( cosmic_dY->GetBinContent(i,j,k)  );
                    denominator = 1.0;
                    if(!isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                //use the MC.
                else{
                    //std::cout << laser_dX->GetBinContent(i,j,k) << " " << laser_dX_err->GetBinContent(i,j,k) << std::endl;
                    numerator = truth_dY->GetBinContent(i,j,k);
                    denominator = 1.0;
                    //in this case, error is ginormous
                    error = 0.0;
                }
                
                if(denominator != 0){
                    combine_dY->SetBinContent(i,j,k, numerator/denominator);
                    combine_dY_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;

                useCosmic = goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k) );
                useLaser  = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && i < highX && i >= lowX && k >= lowZ && k < highZ && j >= lowY && j < highY);
                if(SmoothBoundary && (atXLowerBoundary || atXUpperBoundary) ){
                    numerator = ( cosmic_dZ->GetBinContent(i,j,k) + laser_dZ->GetBinContent(i,j,k) + cosmic_dZ->GetBinContent(i+1,j,k) + laser_dZ->GetBinContent(i+1,j,k) + cosmic_dZ->GetBinContent(i-1,j,k) + laser_dZ->GetBinContent(i-1,j,k) );
                    denominator = 6.0;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(SmoothBoundary && (atYLowerBoundary || atYUpperBoundary) ){
                    numerator = ( cosmic_dZ->GetBinContent(i,j,k) + laser_dZ->GetBinContent(i,j,k) + cosmic_dZ->GetBinContent(i,j+1,k) + laser_dZ->GetBinContent(i,j+1,k) + cosmic_dZ->GetBinContent(i,j-1,k) + laser_dZ->GetBinContent(i,j-1,k) );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(SmoothBoundary && (atZLowerBoundary || atZUpperBoundary) ){
                    numerator = ( cosmic_dZ->GetBinContent(i,j,k) + laser_dZ->GetBinContent(i,j,k) + cosmic_dZ->GetBinContent(i,j,k+1) + laser_dZ->GetBinContent(i,j,k+1) + cosmic_dZ->GetBinContent(i,j,k-1) + laser_dZ->GetBinContent(i,j,k-1) );
                    denominator = 6.0;
                    //std::cout << "At Z Boundary" << std::endl;
                    if(isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                 
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                 
                 
                 
                }
                
                else if(useLaser){
                    numerator = laser_dZ->GetBinContent(i,j,k);
                    denominator = 1.0;
                    if(!isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                else if(useCosmic){
                    numerator = ( cosmic_dZ->GetBinContent(i,j,k)  );
                    denominator = 1.0;
                    if(!isData){
                        //should be a cool calculation
                        error = -1.0;
                    }
                    
                    else{
                        //in this case, error is meaningless
                        error = -1.0;
                    }
                }
                
                //I think this needs to be flipped
                
                
                //use the MC.
                else{
                    //std::cout << laser_dX->GetBinContent(i,j,k) << " " << laser_dX_err->GetBinContent(i,j,k) << std::endl;
                    numerator = truth_dZ->GetBinContent(i,j,k);
                    denominator = 1.0;
                    //in this case, error is ginormous
                    error = -1.0;
                }
                
                if(denominator != 0){
                    combine_dZ->SetBinContent(i,j,k, numerator/denominator);
                    combine_dZ_err->SetBinContent(i,j,k, error);
                    
                }
                
                else
                    std::cout << "Denomiantor for combination == 0! This should not happen!" << std::endl;
               
                
                
                
                
            }//end of loop over y bins
        }//end of loop over x bins
        
    }//end of loop over z bins
    
    
    
    
    
    for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",k),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
            for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
                combine_2D_dX.SetBinContent(i,j,combine_dX->GetBinContent(i,j,k));
                combine_2D_dY.SetBinContent(i,j,combine_dY->GetBinContent(i,j,k));
                combine_2D_dZ.SetBinContent(i,j,combine_dZ->GetBinContent(i,j,k));
                
            }
        }
        drawPlanarPlot(combine_2D_dX, k, "Combined #DeltaX", "combined_2D_dX", axisType::zAxis, 5.0);
        drawPlanarPlot(combine_2D_dY, k, "Combined #DeltaY", "combined_2D_dY", axisType::zAxis);
        drawPlanarPlot(combine_2D_dZ, k, "Combined #DeltaZ", "combined_2D_dZ", axisType::zAxis);
    }
    
    for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
        TH2F combine_2D_dX(Form("combined_2D_dX_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F combine_2D_dY(Form("combined_2D_dY_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        TH2F combine_2D_dZ(Form("combined_2D_dZ_%d",j),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                combine_2D_dX.SetBinContent(k,i,combine_dX->GetBinContent(i,j,k));
                combine_2D_dY.SetBinContent(k,i,combine_dY->GetBinContent(i,j,k));
                combine_2D_dZ.SetBinContent(k,i,combine_dZ->GetBinContent(i,j,k));
                
            }
        }
        drawPlanarPlot(combine_2D_dX, j, "Combined #DeltaX", "combined_2D_dX", axisType::yAxis, 5.0);
        drawPlanarPlot(combine_2D_dY, j, "Combined #DeltaY", "combined_2D_dY", axisType::yAxis);
        drawPlanarPlot(combine_2D_dZ, j, "Combined #DeltaZ", "combined_2D_dZ", axisType::yAxis);
    }
    
    outputHistos->Write();
    outputHistos->Close();
    
    delete fileLaser;
    delete outputHistos;
    delete fileCosmic;
    delete fileTruth;
    
}

void eFieldCalculator::MakeDistorionTree(){
    std::string inputHistograms = "/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsSmoothCosmicAndLaser.root";
    std::string inputEField     = "/uboone/data/users/joelam/SCEInputFiles/Emap-NTT-1-MergedMapsCosmicAndLaserData_.root";
    std::string outputFileName      = "SCEDistortionsMerged.root";
    TFile* fileInput = new TFile(inputHistograms.c_str());
    TFile* fileInputEField = new TFile(inputEField.c_str());
    
    TH3F* combine_dX = (TH3F*) fileInput->Get("combined_dX");
    TH3F* combine_dY = (TH3F*) fileInput->Get("combined_dY");
    TH3F* combine_dZ = (TH3F*) fileInput->Get("combined_dZ");
    
    TH3F* E_field_X = (TH3F*) fileInputEField->Get("Distorted_EField_X");
    TH3F* E_field_Y = (TH3F*) fileInputEField->Get("Distorted_EField_Y");
    TH3F* E_field_Z = (TH3F*) fileInputEField->Get("Distorted_EField_Z");
    
    const int eFieldVoxeslX = 20;
    const int eFieldVoxeslY = 20;
    const int eFieldVoxeslZ = 80;
    
    const double eFieldVoexlWdith = 0.125;
    const double kVtoVolts = 1000.0;
    const double metersTocm = 100.0;
    const double cmToMeters = 0.01;
    const double E0 = 0.273; // kV / cm
    
    
    double x_true, y_true, z_true;
    double x_reco, y_reco, z_reco;
    double x_fwd, y_fwd, z_fwd;
    double x_bkwd, y_bkwd, z_bkwd;
    double Dx, Dy, Dz;
    double Ex, Ey, Ez;
    int elecFate;
    
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *T_calib = new TTree("SpaCEtree_bkwdDisp","SpaCEtree_bkwdDisp");
    T_calib->Branch("x_true",&x_true,"data_bkwdDisp/D");
    T_calib->Branch("y_true",&y_true,"data_bkwdDisp/D");
    T_calib->Branch("z_true",&z_true,"data_bkwdDisp/D");
    T_calib->Branch("x_reco",&x_reco,"data_bkwdDisp/D");
    T_calib->Branch("y_reco",&y_reco,"data_bkwdDisp/D");
    T_calib->Branch("z_reco",&z_reco,"data_bkwdDisp/D");
    
    T_calib->Branch("Dx",&Dx,"data_bkwdDisp/D");
    T_calib->Branch("Dy",&Dy,"data_bkwdDisp/D");
    T_calib->Branch("Dz",&Dz,"data_bkwdDisp/D");
    T_calib->Branch("elecFate",&elecFate,"data_bkwdDisp/I");
    
    TTree *T_field_calib = new TTree("SpaCEtree","SpaCEtree");
    T_field_calib->Branch("xpoint",&x_reco,"data/D");
    T_field_calib->Branch("ypoint",&y_reco,"data/D");
    T_field_calib->Branch("zpoint",&z_reco,"data/D");
    T_field_calib->Branch("Ex",&Ex,"data/D");
    T_field_calib->Branch("Ey",&Ey,"data/D");
    T_field_calib->Branch("Ez",&Ez,"data/D");
    T_field_calib->SetDirectory(outputFile);
    
    

    double xBins[combine_dX->GetNbinsX()];
    combine_dX->GetXaxis()->GetLowEdge(xBins);
    double yBins[combine_dX->GetNbinsY()];
    combine_dX->GetYaxis()->GetLowEdge(yBins);
    double zBins[combine_dX->GetNbinsZ()];
    combine_dX->GetZaxis()->GetLowEdge(zBins);
    
    for(int i = 0; i <= eFieldVoxeslX; ++i){
        x_reco = i*eFieldVoexlWdith;
        double x_interp = x_reco;
        if(i == 0)
            x_interp = (i+1)*eFieldVoexlWdith;
        else if( i == eFieldVoxeslX)
            x_interp = (i-1)*eFieldVoexlWdith;
        
        for(int j = 0; j <= eFieldVoxeslY; ++j){
            y_reco = j*eFieldVoexlWdith;
            double y_interp = y_reco;
            if(j == 0)
                y_interp = (j+1)*eFieldVoexlWdith;
            else if( j == eFieldVoxeslY)
                y_interp = (j-1)*eFieldVoexlWdith;
            for(int k = 0; k <= eFieldVoxeslZ; ++k){
                z_reco = k*eFieldVoexlWdith;
                double z_interp = z_reco;
                if(k == 0)
                    z_interp = (k+1)*eFieldVoexlWdith;
                else if(k == eFieldVoxeslZ)
                    z_interp = (k-1)*eFieldVoexlWdith;
                if(E_field_X->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)) > 10e10){
                    std::cout << "Bogus Point: " << x_interp << " " << y_interp << " " << z_interp << std::endl;
                    Ex = 0.0;
                }
                
                else
                    Ex = kVtoVolts*metersTocm*(E0 - E_field_X->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)));
                if(E_field_Y->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)) > 10e10){
                    std::cout << "Bogus Point: " << x_interp << " " << y_interp << " " << z_interp << std::endl;
                    Ey = 0.0;
                }
                
                else
                    Ey = kVtoVolts*metersTocm*(E_field_Y->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)));
                if(E_field_Z->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)) > 10e10){
                    std::cout << "Bogus Point: " << x_interp << " " << y_interp << " " << z_interp << std::endl;
                    Ez = 0.0;
                }
                
                else
                    Ez = kVtoVolts*metersTocm*(E_field_Z->Interpolate(doInvCoordTransformX(x_interp), doInvCoordTransformY(y_interp), doInvCoordTransformZ(z_interp)));
                //std::cout << Ex << " " << x_reco << " " << doInvCoordTransformX(x_reco) << " " << Ey << " " << y_reco << " " << doInvCoordTransformY(y_reco) << " " << Ez << " " << z_reco << " " << doInvCoordTransformZ(z_reco) << std::endl;
                T_field_calib->Fill();
            }
        }
    }
    
    
    for(int i = combine_dX->GetNbinsX(); i > 0; i--){
        
        if(i!=combine_dX->GetNbinsX())
          x_reco = doCoordTransformX((xBins[i-1] + xBins[i])/2.0);
        else
            x_reco = doCoordTransformX(256);
        //  for(int i = 1; i <= combine_dX->GetNbinsX() ; i++){
        for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
            if(j != combine_dX->GetNbinsY())
               y_reco = doCoordTransformY((yBins[j-1] + yBins[j])/2.0);
            else
                y_reco = doCoordTransformY(116.5);
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                if(k != combine_dX->GetNbinsZ())
                    z_reco = doCoordTransformZ((zBins[k-1] + zBins[k])/2.0);
                else
                    z_reco = doCoordTransformZ(1037);
            
                Dx = -1.0*cmToMeters*combine_dX->GetBinContent(i, j, k);
                Dy = cmToMeters*combine_dY->GetBinContent(i, j, k);
                Dz = cmToMeters*combine_dZ->GetBinContent(i, j, k);

                x_true = x_reco + Dx;
                y_true = y_reco + Dy;
                z_true = z_reco + Dz;
                
                elecFate = 1;
                T_calib->Fill();
                //T_field_calib->Fill();
    
            }
        }
    }
    
    
    
    /*
    TTree *T_true = new TTree("SpaCEtree_fwdDisp","SpaCEtree_fwdDisp");
    T_true->Branch("x_true",&x_true,"data_calib/D");
    T_true->Branch("y_true",&y_true,"data_calib/D");
    T_true->Branch("z_true",&z_true,"data_calib/D");
    T_true->Branch("x_reco",&x_reco,"data_calib/D");
    T_true->Branch("y_reco",&y_reco,"data_calib/D");
    T_true->Branch("z_reco",&z_reco,"data_calib/D");
    T_true->Branch("Dx",&Dx,"data_calib/D");
    T_true->Branch("Dy",&Dy,"data_calib/D");
    T_true->Branch("Dz",&Dz,"data_calib/D");
    T_true->Branch("elecFate",&elecFate,"data_calib/I");
    T_true->SetDirectory(outputFile);
    
    for(int i = 1; i <= combine_dX->GetNbinsX(); i++){
        
        for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                x_reco = (xBins[i-1] + xBins[i])/2.0;
                y_reco = (yBins[j-1] + yBins[j])/2.0;
                z_reco = (zBins[k-1] + zBins[k])/2.0;
                
                Dx = -1.0*combine_dX->GetBinContent(i, j, k);
                Dy = -1.0*combine_dY->GetBinContent(i, j, k);
                Dz = -1.0*combine_dZ->GetBinContent(i, j, k);
                
                x_true = x_reco + Dx;
                y_true = y_reco + Dy;
                z_true = z_reco + Dz;
                
                elecFate = 1;
                T_calib->Fill();
                
            }
        }
    }*/
    T_field_calib->Write();
    T_calib->Write();
    outputFile->Close();
    
    
}

void eFieldCalculator::MakeDistortionHistograms(bool isFwd){
    std::string inputHistograms;
    if(isFwd)
      inputHistograms = "/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsFwdCosmicOnly.root";
    else
      inputHistograms = "/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsCosmicOnlySmoothed.root";
    
    std::string inputEField     = "/uboone/data/users/joelam/SCEInputFiles/Emap-NTT-1-MergedMapsCosmicOnlySmoothed_newCurve.root";
    std::string outputFileName;
    if(isFwd)
      outputFileName = "/uboone/data/users/joelam/SCEDistortionMaps/SCEDistortionsMergedFwdCosmicOnlyE273.root";
    else
      outputFileName = "/uboone/data/users/joelam/SCEDistortionMaps/SCEDistortionsCosmicOnlyE273.root";
    TFile* fileInput = new TFile(inputHistograms.c_str());
    TFile* fileInputEField = new TFile(inputEField.c_str());
    
    TH3F* combine_dX;
    TH3F* combine_dY;
    TH3F* combine_dZ;
    
    if(isFwd){
       combine_dX = (TH3F*) fileInput->Get("combined_fwd_dX");
       combine_dY = (TH3F*) fileInput->Get("combined_fwd_dY");
       combine_dZ = (TH3F*) fileInput->Get("combined_fwd_dZ");
    }
    
    else{
        combine_dX = (TH3F*) fileInput->Get("Reco_Displacement_X");
        combine_dY = (TH3F*) fileInput->Get("Reco_Displacement_Y");
        combine_dZ = (TH3F*) fileInput->Get("Reco_Displacement_Z");
    }
    
    
    TH3F* E_field_X = (TH3F*) fileInputEField->Get("Distorted_EField_X");
    TH3F* E_field_Y = (TH3F*) fileInputEField->Get("Distorted_EField_Y");
    TH3F* E_field_Z = (TH3F*) fileInputEField->Get("Distorted_EField_Z");
    
    std::cout << xMin << ", " << xMax << std::endl;
    std::cout << yMin << ", " << yMax << std::endl;
    std::cout << zMin << ", " << zMax << std::endl;

    
    const double kVtoVolts = 1000.0;
    const double metersTocm = 100.0;
    const double cmToMeters = 0.01;
    const double E0 = 0.2739; // kV / cm
    
    const int nXbins = 26;
    const int nYbins = 26;
    const int nZbins = 101;
    
    const double minX = -0.05;
    const double maxX = 2.55;
    const double minY = -0.05;
    const double maxY = 2.55;
    const double minZ = -0.05;
    const double maxZ = 10.05;
    
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    TH3F *hDx = new TH3F("hDx", "hDx", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    TH3F *hDy = new TH3F("hDy", "hDy", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    TH3F *hDz = new TH3F("hDz", "hDz", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    
    TH3F *hEx = new TH3F("hEx", "hEx", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    TH3F *hEy = new TH3F("hEy", "hEy", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    TH3F *hEz = new TH3F("hEz", "hEz", nXbins, minX, maxX, nYbins, minY, maxY, nZbins, minZ, maxZ);
    
    if(1){
       for(int i = combine_dX->GetNbinsX(); i > 0; --i){
           int xBin = (nXbins - i)+1;
           int xEBin = i;
        if(i == combine_dX->GetNbinsX())
            xEBin = xEBin - 1;
        else if(i == 1)
            xEBin = xEBin + 1;
        
        for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
            int yEBin = j;
            if(j == 1)
                yEBin = yEBin + 1;
            else if(j == combine_dX->GetNbinsY())
                yEBin = yEBin - 1;
            for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                int zEBin = k;
                if(k == 1)
                    zEBin = zEBin + 1;
                else if(k == combine_dX->GetNbinsZ())
                    zEBin = zEBin - 1;
               double Dx = -1.0*combine_dX->GetBinContent(xEBin, yEBin, zEBin);
               double Dy = combine_dY->GetBinContent(xEBin, yEBin, zEBin);
               double Dz = combine_dZ->GetBinContent(xEBin, yEBin, zEBin);
               
               double Ex = (E0 - E_field_X->GetBinContent(xEBin, yEBin, zEBin));
               double Ey = E_field_Y->GetBinContent(xEBin, yEBin, zEBin);
               double Ez = E_field_Z->GetBinContent(xEBin, yEBin, zEBin);
                
                std::cout << i << ", " << j << ", " << k << ", " << Dx << ", " << Dy << ", " << Dz <<", " << Ex << ", " << Ey << ", " << Ez << std::endl;
                
               hDx->SetBinContent(xBin, j, k, Dx);
               hDy->SetBinContent(xBin, j, k, Dy);
               hDz->SetBinContent(xBin, j, k, Dz);
                
               hEx->SetBinContent(xBin, j, k, Ex);
               hEy->SetBinContent(xBin, j, k, Ey);
               hEz->SetBinContent(xBin, j, k, Ez);
                

                //T_field_calib->Fill();
                
            }
        }
    }
    }
    else{
        for(int i = 1; i <= combine_dX->GetNbinsX(); ++i){
            int xEBin = i;
            if(i == combine_dX->GetNbinsX())
                xEBin = xEBin - 1;
            else if(i == 1)
                xEBin = xEBin + 1;
            
            for(int j = 1; j <= combine_dX->GetNbinsY(); j++){
                int yEBin = j;
                if(j == 1)
                    yEBin = yEBin + 1;
                else if(j == combine_dX->GetNbinsY())
                    yEBin = yEBin - 1;
                for(int k = 1; k <= combine_dX->GetNbinsZ(); k++){
                    int zEBin = k;
                    if(k == 1)
                        zEBin = zEBin + 1;
                    else if(k == combine_dX->GetNbinsZ())
                        zEBin = zEBin - 1;
                    double Dx = 1.0*combine_dX->GetBinContent(xEBin, yEBin, zEBin);
                    double Dy = combine_dY->GetBinContent(xEBin, yEBin, zEBin);
                    double Dz = combine_dZ->GetBinContent(xEBin, yEBin, zEBin);
                    
                    double Ex = (E0 - E_field_X->GetBinContent(xEBin, yEBin, zEBin));
                    double Ey = E_field_Y->GetBinContent(xEBin, yEBin, zEBin);
                    double Ez = E_field_Z->GetBinContent(xEBin, yEBin, zEBin);
                    
                    std::cout << i << ", " << j << ", " << k << ", " << Dx << ", " << Dy << ", " << Dz <<", " << Ex << ", " << Ey << ", " << Ez << std::endl;
                    
                    hDx->SetBinContent(i, j, k, Dx);
                    hDy->SetBinContent(i, j, k, Dy);
                    hDz->SetBinContent(i, j, k, Dz);
                    
                    hEx->SetBinContent(i, j, k, Ex);
                    hEy->SetBinContent(i, j, k, Ey);
                    hEz->SetBinContent(i, j, k, Ez);
                    
                    
                    //T_field_calib->Fill();
                    
                }
            }
        }
    }
    outputFile->Write();
    outputFile->Close();
    
    
}

float eFieldCalculator::LinInterp(float x, float x1, float x2, float q00, float q01) {
    return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
}



float eFieldCalculator::TrilinInterp(float x, float y, float z, float q000, float q001, float q010, float q011, float q100, float q101, float q110, float q111, float x1, float x2, float y1, float y2, float z1, float z2) {
    float x00 = LinInterp(x, x1, x2, q000, q100);
    float x10 = LinInterp(x, x1, x2, q010, q110);
    float x01 = LinInterp(x, x1, x2, q001, q101);
    float x11 = LinInterp(x, x1, x2, q011, q111);
    float r0 = LinInterp(y, y1, y2, x00, x01);
    float r1 = LinInterp(y, y1, y2, x10, x11);
    
    return LinInterp(z, z1, z2, r0, r1);
}

PCAResults DoPCA(const PointCloud &points) {
    
    TVector3 outputCentroid;
    std::pair<TVector3,TVector3> outputEndPoints;
    float outputLength;
    TVector3 outputEigenValues;
    std::vector<TVector3> outputEigenVecs;
    
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

std::vector<double> eFieldCalculator::studyResults2(std::string inputMapFileName, std::string plotName)
{
    std::vector<double> return_vec;
    bool containsLaser = false;
    bool isMC = false;
    
    if(inputMapFileName.find("MC") != std::string::npos)
        isMC = true;
    if(inputMapFileName.find("Laser") != std::string::npos)
        containsLaser = true;
    
    std::cout << "Is MC: " << isMC << " Is laser: " << containsLaser << std::endl;
    
    //Bin along track segments.
    //Compare to the far edge.
    
    
    
    //Look at larger track slices
    //const double bufferLength = 0.05*100;
    double driftV = 0.001;
    //double driftV = 0.01; //for optimization studies
    if(isMC)
       driftV   = 0.15;
    const double bufferLength = 0.0;
    const double epsilon = 1e-18;
    const bool skipZeroes = false;
    const bool doTH3Int = true;
    
    const int minTrackPoints = 75;
    //const int minTrackPoints = 50; //unit test
    const int numTrackSegPoints = 15;
    const int lastSegment       = 2;
    const int firstSegment      = 0;
    //const int firstSegment      = 15; //30 for anode piercing
    //const int maxTrackSegments  = 2;
    const int maxTrackSegments    = 1;
    const int stepsBack           = 1; //unit test
    //const int stepsBack           = 15;
    const int maxBadPoints        = 5; //normally 5
    
    TGaxis::SetMaxDigits(3);
    
    double numDivisions_x = nCalibDivisions_x;
    double numDivisions_y = nCalibDivisions_y;
    double numDivisions_z = nCalibDivisions_z;
    /*
    const double zLow = 0.0;
    const double zHigh = Lz;
    */
    
    const double xLow  = 0.0;
    //const double xLow  = 0.5;
    const double xHigh = Lx;
    //const double xHigh  = 1.0;
    
    //const double yLow = 0.0;
    //const double yHigh = Ly;
    
    //Remove most of the tail: low == 1.3 high == Ly
    const double yLow = 0.0;
    const double yHigh = Ly;
    
    const double zLow = 2.1;
    const double zHigh = 8.0;
    
   // const double zLow = 0.0;
   // const double zHigh = Lz;
    /*
    Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
    Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
    Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
    Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
    TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
    gStyle->SetNumberContours(255);
    */
    
    
    
    TH1F *origAngHist = new TH1F("origAngHist","",50,0.0,10.0);
    TH1F *corrAngHist = new TH1F("corrAngHist","",50,0.0,10.0);
    
    TH1F *h_xBins[maxTrackSegments];
    TH1F *h_yBins[maxTrackSegments];
    TH1F *h_zBins[maxTrackSegments];
    
    TH1F *h_last_xBins[maxTrackSegments];
    TH1F *h_last_yBins[maxTrackSegments];
    TH1F *h_last_zBins[maxTrackSegments];
    
    
    for(int i=0; i<maxTrackSegments; ++i){
        h_xBins[i] = new TH1F(Form("xBins_%d", i),"",nCalibDivisions_x+1,0.0,25.0);
        h_yBins[i] = new TH1F(Form("yBins_%d", i),"",nCalibDivisions_y+1,0.0,25.0);
        h_zBins[i] = new TH1F(Form("zBins_%d", i),"",nCalibDivisions_z+1,0.0,100.0);
        
        h_last_xBins[i] = new TH1F(Form("last_xBins_%d", i),"",nCalibDivisions_x+1,0.0,25.0);
        h_last_yBins[i] = new TH1F(Form("last_yBins_%d", i),"",nCalibDivisions_y+1,0.0,25.0);
        h_last_zBins[i] = new TH1F(Form("last_zBins_%d", i),"",nCalibDivisions_z+1,0.0,100.0);
        
        
    }
    
    TH1F *origAngHistByZ[nCalibDivisions_z+1];
    TH1F *corrAngHistByZ[nCalibDivisions_z+1];
    TH1F *origAngHistByX[nCalibDivisions_x+1];
    TH1F *corrAngHistByX[nCalibDivisions_x+1];
    TH1F *origAngHistByY[nCalibDivisions_y+1];
    TH1F *corrAngHistByY[nCalibDivisions_y+1];
    
    TH2F *PointsXY[nCalibDivisions_z+1];
    TH2F *distPointsXY[nCalibDivisions_z+1];
    
    
    TH2F *PointsYZ[nCalibDivisions_x+1];
    TH2F *distPointsYZ[nCalibDivisions_x+1];
    
    
    TH2F *PointsZX[nCalibDivisions_y+1];
    TH2F *distPointsZX[nCalibDivisions_y+1];
    
    
    
    for(int i =0; i < nCalibDivisions_z+1; ++i){
        origAngHistByZ[i] = new TH1F(Form("origAngHistZ_%d", i),"",50,0.0,10.0);
        corrAngHistByZ[i] = new TH1F(Form("corrAngHistZ_%d", i),"",50,0.0,10.0);
        PointsXY[i]      = new TH2F(Form("PointsXY_%d", i),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        distPointsXY[i]      = new TH2F(Form("distPointsXY_%d", i),"",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax);
        
        
    }
    
    for(int i =0; i < nCalibDivisions_x+1; ++i){
        origAngHistByX[i] = new TH1F(Form("origAngHistX_%d", i),"",50,0.0,10.0);
        corrAngHistByX[i] = new TH1F(Form("corrAngHistX_%d", i),"",50,0.0,10.0);
        PointsYZ[i]       = new TH2F(Form("PointsYZ_%d", i),"",nCalibDivisions_z+1,zMin,zMax,25,-125.0,125.0);
        //std::cout << zMin << " " << zMax << " " << yMin << " " << yMax <<std::endl;
        distPointsYZ[i]   = new TH2F(Form("distPointsYZ_%d", i),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_y+1,yMin,yMax);
        
        
    }
    
    for(int i =0; i < nCalibDivisions_y+1; ++i){
        origAngHistByY[i] = new TH1F(Form("origAngHistY_%d", i),"",50,0.0,10.0);
        corrAngHistByY[i] = new TH1F(Form("corrAngHistY_%d", i),"",50,0.0,10.0);
        PointsZX[i]      = new TH2F(Form("PointsZX_%d", i),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        distPointsZX[i]      = new TH2F(Form("distPointsZX_%d", i),"",nCalibDivisions_z+1,zMin,zMax,nCalibDivisions_x+1,xMin,xMax);
        
        
    }
    
    std::string meanFileName = "3DTracks.root";
    /*
    if(skipLaser){
       inputMapFileName = "/uboone/data/users/joelam/SCEInputFiles/MergedMapsCosmicOnly.root";
       plotName = "AnglePlots/CombinedAngCosmicOnlyNoZeroes.png";
       meanFileName = "MeansCosmicOnly.root";
    }
    
    else{
        std::cout << "Using Laser..." << std::endl;
        inputMapFileName = "MergedMapsSmoothCosmicAndLaserMC.root";
        //inputMapFileName = "CombinedTmp.root";
        plotName = "AnglePlots/MergedMapsSmoothCosmicAndLaserMC.png";
        meanFileName
    }*/

    
    TFile *distortionMapInput = new TFile(inputMapFileName.c_str());
    /*TH3F* dist_dX = (TH3F*) distortionMapInput->Get("combined_dX");
    TH3F* dist_dY = (TH3F*) distortionMapInput->Get("combined_dY");
    TH3F* dist_dZ = (TH3F*) distortionMapInput->Get("combined_dZ");
    
    TH3F* dist_dX_err = (TH3F*) distortionMapInput->Get("combined_dX_Error");
    TH3F* dist_dY_err = (TH3F*) distortionMapInput->Get("combined_dY_Error");
    TH3F* dist_dZ_err = (TH3F*) distortionMapInput->Get("combined_dZ_Error");
    */
    
    TH3F* dist_dX  = (TH3F*) distortionMapInput->Get("Reco_Displacement_X");
    TH3F* dist_dY  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Y");
    TH3F* dist_dZ  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Z");
    
    TH3F* dist_dX_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_X_Error");
    TH3F* dist_dY_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_Y_Error");
    TH3F* dist_dZ_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_Z_Error");
    
    
    double corr_Dx[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    double corr_Dy[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    double corr_Dz[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    
    double origAngMean[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    double origAngNum[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    
    double corrAngMean[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    double corrAngNum[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    
    
    std::string TrackFileName = "/uboone/data/users/mrmooney/ForJoel/output_data_tracks.root";
   // if(isMC)
   //     TrackFileName = "/uboone/data/users/mrmooney/ForJoel/output_MC_tracks.root";
    
    //TFile* inputFile = new TFile("output.root");
    TFile* inputFile = new TFile(TrackFileName.c_str());
    //TFile* inputFile = new TFile("/uboone/data/users/mrmooney/ForJoel/output_MC_tracks.root");
   // TFile* inputFile = new TFile("Naive_straight_tracks_w_hits.root");
   //TFile* inputFile = new TFile("straight_tracks_w_distortion_v1.root");
   //TFile* inputFile = new TFile("/uboone/data/users/joelam/SCEInputFiles/LaserTracksData.root");
    

    TTreeReader readerTracks("SpaCEtree_tracks", inputFile);
    TTreeReaderValue<Int_t> nElec_tracks(readerTracks, "nElec_tracks");
    TTreeReaderArray<Double_t> elecX_tracks_orig(readerTracks, "elecX_tracks");
    TTreeReaderArray<Double_t> elecY_tracks_orig(readerTracks, "elecY_tracks");
    TTreeReaderArray<Double_t> elecZ_tracks_orig(readerTracks, "elecZ_tracks");
    
    //std::cout << dist_dX->GetNbinsX() << " " << dist_dX->GetNbinsY() << " " << dist_dX->GetNbinsZ() << std::endl;
    
    for (int x = 0; x <dist_dX->GetNbinsX(); x++ ) {
        for (int y = 0; y <dist_dX->GetNbinsY(); y++ ) {
            for (int z = 0; z <dist_dX->GetNbinsZ(); z++ ) {
                corr_Dx[x][y][z] = 0.0;
                corr_Dy[x][y][z] = 0.0;
                corr_Dz[x][y][z] = 0.0;
                
            }
        }
    }
    
    for (int x = 0; x <nCalibDivisions_x+1; x++ ) {
        for (int y = 0; y <nCalibDivisions_y+1; y++ ) {
            for (int z = 0; z <nCalibDivisions_z+1; z++ ) {
                origAngMean[x][y][z] = 0.0;
                origAngNum[x][y][z]  = 0.0;
                
                corrAngMean[x][y][z] = 0.0;
                corrAngNum[x][y][z]  = 0.0;
               
            }
        }
    }
    
    //Change this to loop over TH3 Bins
    for(int i = 0; i < dist_dX->GetNbinsX(); i++){
        for(int j = 0; j < dist_dX->GetNbinsY(); j++){
            for(int k = 0; k < dist_dX->GetNbinsZ(); k++){
                corr_Dx[i][j][k] = -0.01*dist_dX->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1);
                /*if(!goodLaser(dist_dX->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1), dist_dX_err->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1))){
                    corr_Dx[i][j][k] = 0.0;
                    //std::cout << "Bad point Dx: " <<  i << " " << " " << j << " " << k << std::endl;
                }*/
                corr_Dy[i][j][k] = 0.01*dist_dY->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1);
                /*if(!goodLaser(dist_dY->GetBinContent(dist_dY->GetNbinsX()-i,j+1,k+1), dist_dY_err->GetBinContent(dist_dY->GetNbinsX()-i,j+1,k+1))){
                    corr_Dy[i][j][k] = 0.0;
                    //std::cout << "Bad point Dy: " <<  i << " " << " " << j << " " << k << std::endl;
                }*/
                corr_Dz[i][j][k] = 0.01*dist_dZ->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1);
                /*if(!goodLaser(dist_dZ->GetBinContent(dist_dZ->GetNbinsX()-i,j+1,k+1), dist_dZ_err->GetBinContent(dist_dZ->GetNbinsX()-i,j+1,k+1))){
                    corr_Dz[i][j][k] = 0.0;
                    //std::cout << "Bad point Dz: " <<  i << " " << " " << j << " " << k << std::endl;
                }*/
                
            }
        }
    }
    
    TH3F* diff_dX = new TH3F("diff_dX","",nCalibDivisions_x+1,xMin,xMax,nCalibDivisions_y+1,yMin,yMax,nCalibDivisions_z+1,zMin,zMax);
    
    //TFile* fileMeans = new TFile(meanFileName.c_str(), "RECREATE");
    
    const int nTracks = 10;
    int trackNo       = 0;
    TH3F* tracksOrig[nTracks];
    TH3F* tracksCorr[nTracks];
    
    /*for(int i = 0; i < nTracks; ++i){
       tracksOrig[i] = new TH3F(Form("tracksOrig%d", i),"",1000,0,10.0,250,0,2.5,250,0,2.5);
       tracksCorr[i] = new TH3F(Form("tracksCorr%d", i),"",1000,0,10.0,250,0,2.5,250,0,2.5);
       tracksOrig[i]->GetXaxis()->SetTitle("Z");
       tracksOrig[i]->GetYaxis()->SetTitle("X");
       tracksOrig[i]->GetZaxis()->SetTitle("Y");
        
       tracksCorr[i]->GetXaxis()->SetTitle("Z");
       tracksCorr[i]->GetYaxis()->SetTitle("X");
       tracksCorr[i]->GetZaxis()->SetTitle("Y");
        
       tracksOrig[i]->SetMarkerStyle(kFullDotLarge);
       tracksCorr[i]->SetMarkerStyle(kFullDotLarge);
       tracksOrig[i]->SetMarkerColor(kRed);
       tracksCorr[i]->SetMarkerColor(kGreen+2);
        
        
    }*/
    int counter = 0;
    //while (readerTracks.Next())
    //while ((readerTracks.Next()) && (counter < 100000))
    while ((readerTracks.Next()) && (counter < 213000))
    {
        //std::cout << "Go go go" << std::endl;
        if(counter % 10000 == 0) {std::cout << counter << std::endl;}
        counter++;
        
        if (*nElec_tracks >= minTrackPoints) {
          for(int trackSeg = 0; trackSeg < maxTrackSegments; ++trackSeg){
            PointCloud startPoints;
            PointCloud endPoints;
            
            std::vector<double> elecX_tracks;
            std::vector<double> elecY_tracks;
            std::vector<double> elecZ_tracks;
   
            for(int i = 0; i < *nElec_tracks; i++){
                //std::vector<double> points_laser;
                //std::vector<double> points_cosmic;
                double driftVFactor = cosmicDriftVScale;
                /*if(isInLaserRegion(doInvCoordTransformX(elecX_tracks_orig[i] ), doInvCoordTransformY(elecY_tracks_orig[i]), doInvCoordTransformZ(elecZ_tracks_orig[i]) ) )
                    driftVFactor = laserDriftVScale;
                else
                    driftVFactor = cosmicDriftVScale;
                */
                const double pointOffset = driftV*driftVFactor*(Lx-elecX_tracks_orig[i]);
              
                elecX_tracks.push_back(elecX_tracks_orig[i]  + pointOffset);
                elecY_tracks.push_back(elecY_tracks_orig[i]);
                elecZ_tracks.push_back(elecZ_tracks_orig[i]);
                
                /*
                elecX_tracks.push_back(doCoordTransformX(elecX_tracks_orig[i]) );
                elecY_tracks.push_back(doCoordTransformY(elecY_tracks_orig[i]) );
                elecZ_tracks.push_back(doCoordTransformZ(elecZ_tracks_orig[i]) );
                */
                //std::cout << i << " " << elecX_tracks.at(i) << " " << elecY_tracks.at(i) << " " << elecZ_tracks.at(i) << std::endl;
            }
              
              int xBin = TMath::Nint(elecX_tracks[firstSegment]*numDivisions_x/Lx);
              int yBin = TMath::Nint(elecY_tracks[firstSegment]*numDivisions_y/Ly);
              int zBin = TMath::Nint(elecZ_tracks[firstSegment]*numDivisions_z/Lz);
              
              int last_xBin = TMath::Nint(elecX_tracks[*nElec_tracks-stepsBack]*numDivisions_x/Lx);
              int last_yBin = TMath::Nint(elecY_tracks[*nElec_tracks-stepsBack]*numDivisions_y/Ly);
              int last_zBin = TMath::Nint(elecZ_tracks[*nElec_tracks-stepsBack]*numDivisions_z/Lz);
              
              if(xBin > nCalibDivisions_x)
                  xBin = nCalibDivisions_x;
              if(yBin > nCalibDivisions_y)
                  yBin = nCalibDivisions_y;
              if(zBin > nCalibDivisions_z)
                  zBin = nCalibDivisions_z;
              
              if(xBin < 0 )
                  xBin = 0;
              if(yBin < 0)
                  yBin = 0.0;
              if(zBin < 0)
                  zBin = 0.0;
              
              if(last_xBin > nCalibDivisions_x)
                  last_xBin = nCalibDivisions_x;
              if(last_yBin > nCalibDivisions_y)
                  last_yBin = nCalibDivisions_y;
              if(last_zBin > nCalibDivisions_z)
                  last_zBin = nCalibDivisions_z;
              
              if(last_xBin < 0 )
                  last_xBin = 0;
              if(last_yBin < 0)
                  last_yBin = 0.0;
              if(last_zBin < 0)
                  last_zBin = 0.0;

            
            
            Int_t numBadPoints_start = 0;
              
              for (int i = firstSegment; i < numTrackSegPoints+firstSegment; i++) {
                //Calculate points for tri-linear interp here
                  //Calculate points for tri-linear interp here
                  float x_center = elecX_tracks[i];
                  float y_center = elecY_tracks[i];
                  float z_center = elecZ_tracks[i];
                  
                  float dXddd = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXddu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXdud = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXduu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXudd = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXudu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXuud = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXuuu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float dYddd = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYddu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYdud = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYduu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYudd = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYudu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYuud = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYuuu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float dZddd = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZddu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZdud = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZduu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZudd = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZudu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZuud = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZuuu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float x_up = TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                  float x_down = TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                  float y_up   = TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                  float y_down = TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                  float z_up = TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                  float z_down = TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                  
                  float xDist = TrilinInterp(x_center, y_center, z_center, dXddd, dXddu, dXdud, dXduu, dXudd, dXudu, dXuud, dXuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  float yDist = TrilinInterp(x_center, y_center, z_center, dYddd, dYddu, dYdud, dYduu, dYudd, dYudu, dYuud, dYuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  float zDist = TrilinInterp(x_center, y_center, z_center, dZddd, dZddu, dZdud, dZduu, dZudd, dZudu, dZuud, dZuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  
                  float xDistInt = -0.01*dist_dX->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      xDist = xDistInt;
                  float yDistInt = 0.01*dist_dY->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      yDist = yDistInt;
                  float zDistInt = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      zDist = zDistInt;
                  
                  
                if ((elecX_tracks[i] <= xLow) || (elecX_tracks[i] >= xHigh) || (elecY_tracks[i] <= yLow) || (elecY_tracks[i] >= yHigh) || (elecZ_tracks[i] <= zLow) || (elecZ_tracks[i] >= zHigh)) {
                    //std::cout << "Bad point: " << elecX_tracks[i] << " " << elecY_tracks[i] << " " << elecZ_tracks[i] << std::endl;
                    numBadPoints_start++;
                    continue;
                }
                  
                if(skipZeroes && (fabs(xDist) < epsilon || fabs(yDist) < epsilon || fabs(zDist) < epsilon ) ){
                      //std::cout << "Bad point: " << xDist << " " << yDist << " " << zDist << std::endl;
                      numBadPoints_start++;
                      continue;
                    
                  }
                  
                if(xDist>1E10 || yDist>1E10 || zDist>1E10){
                    //std::cout << "Bad distortion: " << xDist << " " << yDist << " " << zDist << std::endl;
                    numBadPoints_start++;
                    continue;
                
                  }
                  Point tempPoint;
                tempPoint.x = elecX_tracks[i];
                tempPoint.y = elecY_tracks[i];
                tempPoint.z = elecZ_tracks[i];
                
                int xBinI = TMath::Nint(elecX_tracks[i]*numDivisions_x/Lx);
                int yBinI = TMath::Nint(elecY_tracks[i]*numDivisions_y/Ly);
                int zBinI = TMath::Nint(elecZ_tracks[i]*numDivisions_z/Lz);
                
                std::vector<double> points_cosmic;
                std::vector<double> points_laser;
                points_laser.push_back(tempPoint.x);
                points_laser.push_back(tempPoint.y);
                points_laser.push_back(tempPoint.z);
                points_cosmic = cosmicToLaser(points_laser);
                //if(points_cosmic[1] < 0.0) std::cout << xBinI<< std::endl;
                PointsXY[zBinI]->Fill(points_cosmic[0], points_cosmic[1]);
                PointsYZ[xBinI]->Fill(points_cosmic[2], points_cosmic[1]);
                PointsZX[yBinI]->Fill(points_cosmic[2], points_cosmic[0]);
                
                startPoints.push_back(tempPoint);
            }
            if (numBadPoints_start > maxBadPoints) continue;
            
            Int_t numBadPoints_end = 0;
            for (int i = *nElec_tracks-stepsBack; i > *nElec_tracks-stepsBack-numTrackSegPoints; i--) {
                float x_center = elecX_tracks[i];
                float y_center = elecY_tracks[i];
                float z_center = elecZ_tracks[i];
                
                float dXddd = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXddu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXdud = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXduu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXudd = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXudu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXuud = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXuuu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float dYddd = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYddu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYdud = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYduu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYudd = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYudu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYuud = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYuuu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float dZddd = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZddu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZdud = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZduu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZudd = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZudu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZuud = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZuuu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float x_up = TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                float x_down = TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                float y_up   = TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                float y_down = TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                float z_up = TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                float z_down = TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                
                float xDist = TrilinInterp(x_center, y_center, z_center, dXddd, dXddu, dXdud, dXduu, dXudd, dXudu, dXuud, dXuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                float yDist = TrilinInterp(x_center, y_center, z_center, dYddd, dYddu, dYdud, dYduu, dYudd, dYudu, dYuud, dYuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                float zDist = TrilinInterp(x_center, y_center, z_center, dZddd, dZddu, dZdud, dZduu, dZudd, dZudu, dZuud, dZuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                
                float xDistInt = -0.01*dist_dX->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    xDist = xDistInt;
                float yDistInt = 0.01*dist_dY->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    yDist = yDistInt;
                float zDistInt = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    zDist = zDistInt;
                
                if ((elecX_tracks[i] <= xLow) || (elecX_tracks[i] >= xHigh) || (elecY_tracks[i] <= yLow) || (elecY_tracks[i] >= yHigh) || (elecZ_tracks[i] <= zLow) || (elecZ_tracks[i] >= zHigh)) {
                    numBadPoints_end++;
                    continue;
                }
                
                if(skipZeroes && (fabs(xDist) < epsilon || fabs(yDist) < epsilon || fabs(zDist) < epsilon ) ){
                    //std::cout << "Bad point: " << xDist << " " << yDist << " " << zDist << std::endl;
                    numBadPoints_end++;
                    continue;
                    
                }
                
                if(xDist>1E10 || yDist>1E10 || zDist>1E10){
                    //std::cout << "Bad distortion: " << xDist << " " << yDist << " " << zDist << std::endl;
                    numBadPoints_end++;
                    continue;
                }
                
                
                Point tempPoint;
                tempPoint.x = elecX_tracks[i];
                tempPoint.y = elecY_tracks[i];
                tempPoint.z = elecZ_tracks[i];
                
                endPoints.push_back(tempPoint);
            }
            if (numBadPoints_end > maxBadPoints) continue;
            
              
            PCAResults results_start = DoPCA(startPoints);
            PCAResults results_end = DoPCA(endPoints);
            
            Double_t dotProd = results_start.eVecs[0](0)*results_end.eVecs[0](0) + results_start.eVecs[0](1)*results_end.eVecs[0](1) + results_start.eVecs[0](2)*results_end.eVecs[0](2);
            Double_t startMag = sqrt(pow(results_start.eVecs[0](0),2) + pow(results_start.eVecs[0](1),2) + pow(results_start.eVecs[0](2),2));
            Double_t endMag = sqrt(pow(results_end.eVecs[0](0),2) + pow(results_end.eVecs[0](1),2) + pow(results_end.eVecs[0](2),2));
            
            Double_t dTheta = TMath::ACos(dotProd/(startMag*endMag));
            //std::cout << dTheta << " " << std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)) << std::endl;
            if(trackSeg == 0)
              origAngHist->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByX[xBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByY[yBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByZ[zBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngMean[xBin][yBin][zBin] += std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265));
            origAngNum[xBin][yBin][zBin]  += 1.0;
            /*if(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)) < 2.5){
                  h_xBins[trackSeg]->Fill(xBin);
                  h_yBins[trackSeg]->Fill(yBin);
                  h_zBins[trackSeg]->Fill(zBin);
                  
                  

              }
              
              else{
                  h_last_xBins[trackSeg]->Fill(xBin);
                  h_last_yBins[trackSeg]->Fill(yBin);
                  h_last_zBins[trackSeg]->Fill(zBin);
                  
              }*/
            
            PointCloud startPointsCorr;
            PointCloud endPointsCorr;
            
            Int_t numBadPoints_start_corr = 0;
              
              
              for (int i = firstSegment; i < numTrackSegPoints+firstSegment; i++) {
                  //Calculate points for tri-linear interp here
                  float pointOffset = 0.0;
                  float x_center = elecX_tracks[i];
                  //Do not correct the corrections, only the input tracks (maybe)
                  if(0){
                      //std::cout << "Correcting x point: " << x_center << " " << elecX_tracks[i] - pointOffset << std::endl;
                      double driftVFactor = 1.0;
                      if(isInLaserRegion(doInvCoordTransformX(elecX_tracks_orig[i] ), doInvCoordTransformY(elecY_tracks_orig[i]), doInvCoordTransformZ(elecZ_tracks_orig[i]) ) )
                          driftVFactor = laserDriftVScale;
                      else
                          driftVFactor = cosmicDriftVScale;
                      pointOffset = driftV*driftVFactor*(Lx-elecX_tracks[i]);
                  }
                  float y_center = elecY_tracks[i];
                  float z_center = elecZ_tracks[i];
                  
                  float dXddd = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXddu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXdud = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXduu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXudd = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXudu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXuud = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dXuuu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float dYddd = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYddu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYdud = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYduu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYudd = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYudu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYuud = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dYuuu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float dZddd = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZddu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZdud = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZduu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZudd = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZudu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZuud = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                  float dZuuu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                  
                  float x_up = TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                  float x_down = TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                  float y_up   = TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                  float y_down = TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                  float z_up = TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                  float z_down = TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                  
                  float xDist = TrilinInterp(x_center, y_center, z_center, dXddd, dXddu, dXdud, dXduu, dXudd, dXudu, dXuud, dXuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  float yDist = TrilinInterp(x_center, y_center, z_center, dYddd, dYddu, dYdud, dYduu, dYudd, dYudu, dYuud, dYuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  float zDist = TrilinInterp(x_center, y_center, z_center, dZddd, dZddu, dZdud, dZduu, dZudd, dZudu, dZuud, dZuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                  
                  float xDistInt = -0.01*dist_dX->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      xDist = xDistInt;
                  float yDistInt = 0.01*dist_dY->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      yDist = yDistInt;
                  float zDistInt = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                  if(doTH3Int)
                      zDist = zDistInt;
                  
                  
                  if ((elecX_tracks[i] <= xLow) || (elecX_tracks[i] >= xHigh) || (elecY_tracks[i] <= yLow) || (elecY_tracks[i] >= yHigh ) || (elecZ_tracks[i] <= zLow) || (elecZ_tracks[i] >= zHigh)) {
                    numBadPoints_start_corr++;
                    continue;
                  }
                  
                  if(skipZeroes && (fabs(xDist) < epsilon || fabs(yDist) < epsilon || fabs(zDist) < epsilon ) ){
                      //std::cout << "Bad point: " << xDist << " " << yDist << " " << zDist << std::endl;
                      numBadPoints_start_corr++;
                      continue;
                      
                  }
                  
                  if(xDist>1E10 || yDist>1E10 || zDist>1E10){
                      //std::cout << "Bad distortion: " << xDist << " " << yDist << " " << zDist << std::endl;
                      numBadPoints_start_corr++;
                      continue;
                      
                  }
                
                
               /*
                std::cout << counter << " " << i << "   ";
                std::cout << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << std::endl;
                */
                
                Point tempPoint;
               
               
                tempPoint.x = x_center + xDist - pointOffset;
                tempPoint.y = y_center + yDist;
                tempPoint.z = z_center + zDist;
                
                
                std::vector<double> points_laser;
                std::vector<double> points_cosmic;
                
                points_laser.push_back(tempPoint.x);
                points_laser.push_back(tempPoint.y);
                points_laser.push_back(tempPoint.z);
                points_cosmic = cosmicToLaser(points_laser);
                
                distPointsXY[zBin]->Fill(points_cosmic[0], points_cosmic[1]);
                distPointsYZ[xBin]->Fill(points_cosmic[2], points_cosmic[1]);
                distPointsZX[yBin]->Fill(points_cosmic[2], points_cosmic[0]);
                
                startPointsCorr.push_back(tempPoint);
            }
            if (numBadPoints_start_corr > maxBadPoints) continue;
            
            Int_t numBadPoints_end_corr = 0;
            for (int i = *nElec_tracks-stepsBack; i > *nElec_tracks-stepsBack-numTrackSegPoints; i--) {
                //Calculate points for tri-linear interp here
                double pointOffset = 0.0;
                float x_center = elecX_tracks[i];
                //Do not correct the corrections, only the input tracks (maybe)
                if(0){
                    //std::cout << "Correcting x point: " << x_center << " " << elecX_tracks[i] - pointOffset << std::endl;
                    double driftVFactor = 1.0;
                    if(isInLaserRegion(doInvCoordTransformX(elecX_tracks_orig[i] ), doInvCoordTransformY(elecY_tracks_orig[i]), doInvCoordTransformZ(elecZ_tracks_orig[i]) ) )
                        driftVFactor = laserDriftVScale;
                    else
                        driftVFactor = cosmicDriftVScale;
                    pointOffset = driftV*driftVFactor*(Lx-elecX_tracks[i]);
                    //pointOffset = driftV*(Lx-elecX_tracks_orig[i]);
                }
                float y_center = elecY_tracks[i];
                float z_center = elecZ_tracks[i];
                
                float dXddd = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXddu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXdud = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXduu = corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXudd = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXudu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXuud = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dXuuu = corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float dYddd = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYddu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYdud = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYduu = corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYudd = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYudu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYuud = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dYuuu = corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float dZddd = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZddu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZdud = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZduu = corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZudd = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZudu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZuud = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)];
                float dZuuu = corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)];
                
                float x_up = TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                float x_down = TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x);
                float y_up   = TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                float y_down = TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y);
                float z_up = TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                float z_down = TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z);
                
                float xDist = TrilinInterp(x_center, y_center, z_center, dXddd, dXddu, dXdud, dXduu, dXudd, dXudu, dXuud, dXuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                float yDist = TrilinInterp(x_center, y_center, z_center, dYddd, dYddu, dYdud, dYduu, dYudd, dYudu, dYuud, dYuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                float zDist = TrilinInterp(x_center, y_center, z_center, dZddd, dZddu, dZdud, dZduu, dZudd, dZudu, dZuud, dZuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                
                float xDistInt = -0.01*dist_dX->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    xDist = xDistInt;
                float yDistInt = 0.01*dist_dY->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    yDist = yDistInt;
                float zDistInt = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(x_center), doInvCoordTransformY(y_center), doInvCoordTransformZ(z_center));
                if(doTH3Int)
                    zDist = zDistInt;
                
                if ((elecX_tracks[i] <= xLow) || (elecX_tracks[i] >= xHigh) || (elecY_tracks[i] <=yLow) || (elecY_tracks[i] >= yHigh) || (elecZ_tracks[i] <= zLow) || (elecZ_tracks[i] >= zHigh) ) {
                    numBadPoints_end_corr++;
                    continue;
                }
                
                if(skipZeroes && (fabs(xDist) < epsilon || fabs(yDist) < epsilon || fabs(zDist) < epsilon ) ){
                    //std::cout << "Bad point: " << xDist << " " << yDist << " " << zDist << std::endl;
                    numBadPoints_end_corr++;
                    continue;
                    
                }
                
                if(xDist>1E10 || yDist>1E10 || zDist>1E10){
                    //std::cout << "Bad distortion: " << xDist << " " << yDist << " " << zDist << std::endl;
                    numBadPoints_end_corr++;
                    continue;
                
                }
               /* std::cout << counter << " " << i << "   ";
                std::cout << elecX_tracks.at(i)<< " " << elecY_tracks.at(i) << " " << elecZ_tracks.at(i) << " " << std::endl;
                std::cout << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << std::endl;
                */
                
                Point tempPoint;
                
                tempPoint.x = x_center + xDist - pointOffset;
                tempPoint.y = y_center + yDist;
                tempPoint.z = z_center + zDist;
             

                
                endPointsCorr.push_back(tempPoint);
            }
            if (numBadPoints_end_corr > maxBadPoints) continue;
            
            PCAResults results_start_corr = DoPCA(startPointsCorr);
            PCAResults results_end_corr = DoPCA(endPointsCorr);
            
            Double_t dotProdCorr = results_start_corr.eVecs[0](0)*results_end_corr.eVecs[0](0) + results_start_corr.eVecs[0](1)*results_end_corr.eVecs[0](1) + results_start_corr.eVecs[0](2)*results_end_corr.eVecs[0](2);
            Double_t startMagCorr = sqrt(pow(results_start_corr.eVecs[0](0),2) + pow(results_start_corr.eVecs[0](1),2) + pow(results_start_corr.eVecs[0](2),2));
            Double_t endMagCorr = sqrt(pow(results_end_corr.eVecs[0](0),2) + pow(results_end_corr.eVecs[0](1),2) + pow(results_end_corr.eVecs[0](2),2));
            
            Double_t dThetaCorr = TMath::ACos(dotProdCorr/(startMagCorr*endMagCorr));
            
            //cout << "    ANGLE:  " << min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265));
            if(trackSeg == 0)
               corrAngHist->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByX[xBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByY[yBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByZ[zBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            
            corrAngMean[xBin][yBin][zBin] += std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265));
            corrAngNum[xBin][yBin][zBin]  += 1.0;
            if(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)) > 5.0 && trackNo < nTracks){
                //loop over track points
                for(int point = 0; point < *nElec_tracks; ++point){
                    //std::cout << elecX_tracks[point] << " " << elecY_tracks[point] << " " << elecZ_tracks[point] << std::endl;
                    //tracksOrig[trackNo]->Fill(elecZ_tracks[point], elecX_tracks[point], elecY_tracks[point]);
                    float xDistortion = -0.01*dist_dX->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                    float yDistortion = 0.01*dist_dY->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                    float zDistortion = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                    //tracksCorr[trackNo]->Fill((elecZ_tracks[point]+zDistortion), (elecX_tracks[point]+xDistortion), (elecY_tracks[point]+yDistortion));
                    
                }
                  
                  
                ++trackNo;
            }
              
              
            elecX_tracks.clear();
            elecY_tracks.clear();
            elecZ_tracks.clear();
           
        }//end of loop over track segments
       }//end of if number of track points
    } //end of loop over track segments
    double totalTracks = origAngHist->Integral();
    origAngHist->Scale(1.0/origAngHist->Integral());
    corrAngHist->Scale(1.0/corrAngHist->Integral());
    
    gStyle->SetTitleW(0.9);
    gStyle->SetOptStat(0);
    
    /*
    for (int x = 0; x <nCalibDivisions_x+1; x++ ) {
        for (int y = 0; y <nCalibDivisions_y+1; y++ ) {
            for (int z = 0; z <nCalibDivisions_z+1; z++ ) {
                double corrMean = -1.0;
                if(corrAngNum[x][y][z] >= 10){
                    corrMean = (double)corrAngMean[x][y][z]/corrAngNum[x][y][z];
                    meansCorr->SetBinContent(x+1, y+1, z+1, corrMean);
                }
                double origMean = -1.0;
                if(origAngNum[x][y][z] >= 10){
                    origMean = (double)origAngMean[x][y][z]/origAngNum[x][y][z];
                    meansOrig->SetBinContent(x+1, y+1, z+1, origMean);
                }
                double meanDiff = -1.0;
                if(corrMean > 0.0 && origMean > 0.0)
                   meanDiff = origMean - corrMean;
                meansDiff->SetBinContent(x+1, y+1, z+1, meanDiff);
                numTracks->SetBinContent(x+1, y+1, z+1, corrAngNum[x][y][z]);
                
                
            }
        }
    }*/
    /*
    for(int i=0; i < numDivisions_x; ++i){
        draw1DPlot(origAngHistByX[i], corrAngHistByX[i], i, "Angle", "AnglePlots/combinedAngHistInX", axisType::xAxis);
        drawPlanarPlot(PointsYZ[i], i, "Points", "Points", axisType::xAxis);
    }
    
    for(int i=0; i < numDivisions_y; ++i){
        draw1DPlot(origAngHistByY[i], corrAngHistByY[i], i, "Angle", "AnglePlots/combinedAngHistInY", axisType::yAxis);
    }
    
    for(int i=0; i < numDivisions_z; ++i){
        draw1DPlot(origAngHistByZ[i], corrAngHistByZ[i], i, "Angle", "AnglePlots/combinedAngHistInZ", axisType::zAxis);
    }*/
    
    for(int i=0; i<maxTrackSegments; ++i){
        TCanvas *can = new TCanvas(Form("can"),"",900,900);
        can->cd();
        h_xBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_xBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_xBins[i]->GetXaxis()->SetTitle("X Voxel");
        h_xBins[i]->GetYaxis()->SetTitle("Number of Tracks");
        h_xBins[i]->SetTitle("Start Point Voxel");
        h_xBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_xBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_xBins[i]->SetLineWidth(2.0);
        h_xBins[i]->SetLineColor(kBlack);
        h_xBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/StartXBin_%d.png", i));
        can->Update();
        can->Clear();
        
        h_yBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_yBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_xBins[i]->GetXaxis()->SetTitle("Y Voxel");
        h_xBins[i]->GetYaxis()->SetTitle("Number of Tracks");
        h_xBins[i]->SetTitle("Start Point Voxel");
        h_yBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_yBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_yBins[i]->SetLineWidth(2.0);
        h_yBins[i]->SetLineColor(kBlack);
        h_yBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/StartYBin_%d.png", i));
        can->Update();
        can->Clear();
        
        h_zBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_zBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_xBins[i]->GetXaxis()->SetTitle("Z Voxel");
        h_xBins[i]->GetYaxis()->SetTitle("Number of Tracks");
        h_xBins[i]->SetTitle("Start Point Voxel");
        h_zBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_zBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_zBins[i]->SetLineWidth(2.0);
        h_zBins[i]->SetLineColor(kBlack);
        h_zBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/StartZBin_%d.png", i));
        can->Update();
        can->Clear();
        
        can->cd();
        h_last_xBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_last_xBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_last_xBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_last_xBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_last_xBins[i]->SetLineWidth(2.0);
        h_last_xBins[i]->SetLineColor(kBlack);
        h_last_xBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/LastXBin_%d.png", i));
        can->Update();
        can->Clear();
        
        h_last_yBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_last_yBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_last_yBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_last_yBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_last_yBins[i]->SetLineWidth(2.0);
        h_last_yBins[i]->SetLineColor(kBlack);
        h_last_yBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/LastYBin_%d.png", i));
        can->Update();
        can->Clear();
        
        h_last_zBins[i]->GetXaxis()->SetTitleOffset(0.95);
        h_last_zBins[i]->GetXaxis()->SetTitleSize(0.045);
        h_last_zBins[i]->GetYaxis()->SetTitleOffset(0.95);
        h_last_zBins[i]->GetYaxis()->SetTitleSize(0.045);
        h_last_zBins[i]->SetLineWidth(2.0);
        h_last_zBins[i]->SetLineColor(kBlack);
        h_last_zBins[i]->Draw("hist");
        can->SaveAs(Form("AnglePlots/LastZBin_%d.png", i));
        can->Update();
        can->Clear();
        
    }
    
    TCanvas *can0 = new TCanvas(Form("can0"),"",900,900);
    can0->cd();
    double scale = 1.0;
    h_xBins[0]->Sumw2();
    h_last_xBins[0]->Sumw2();
    if(h_xBins[0]->Integral() > 0.0)
        scale = 1 / (h_xBins[0]->Integral());
    h_xBins[0]->Scale(scale);
    scale = 1.0;
    if(h_last_xBins[0]->Integral() > 0.0)
        scale = 1 / (h_last_xBins[0]->Integral());
    
    h_last_xBins[0]->Scale(scale);
    double maximum = 0.0;
    double headroom = 1.1;
    if(h_xBins[0]->GetMaximum() > maximum )
        maximum = h_xBins[0]->GetMaximum();
    if(h_last_xBins[0]->GetMaximum() > maximum )
        maximum = h_last_xBins[0]->GetMaximum();
    maximum = headroom*maximum;
    h_xBins[0]->SetMarkerStyle(kFullDotMedium);
    h_xBins[0]->SetMarkerColor(kRed+2);
    h_xBins[0]->SetLineColor(kRed+2);
    h_last_xBins[0]->SetMarkerStyle(kFullDotMedium);
    h_last_xBins[0]->SetMarkerColor(kBlue+2);
    h_last_xBins[0]->SetLineColor(kBlue+2);
    h_xBins[0]->SetMaximum(maximum);
    h_xBins[0]->Draw("E0");
    h_last_xBins[0]->Draw("E0 same");
    can0->SaveAs(Form("AnglePlots/XBinByDeflection.png"));
    can0->Update();
    can0->Clear();
    
    TCanvas *c_orig = new TCanvas(Form("can_orig"),"",900,900);
    c_orig->cd();
    origAngHist->GetXaxis()->SetTitle("#Delta #theta [degrees]");
    origAngHist->GetXaxis()->SetTitleOffset(0.95);
    origAngHist->GetXaxis()->SetTitleSize(0.045);
    origAngHist->GetYaxis()->SetTitle("Arb. Units");
    origAngHist->GetYaxis()->SetTitleOffset(0.95);
    origAngHist->GetYaxis()->SetTitleSize(0.05);
    origAngHist->GetYaxis()->SetNoExponent(kTRUE);
    origAngHist->SetStats(0);
    origAngHist->SetLineWidth(2.0);
    origAngHist->SetLineColor(kRed);
    origAngHist->Draw("HIST");
    origAngHist->Draw("AXISsame");
    origAngHist->SetMinimum(0.0001);
    c_orig->SaveAs("origAngHist.png");
    
    TCanvas *c_corr = new TCanvas();
    c_corr->cd();
    corrAngHist->GetXaxis()->SetTitle("#Delta#theta [degrees]");
    corrAngHist->GetXaxis()->SetTitleOffset(0.95);
    corrAngHist->GetXaxis()->SetTitleSize(0.045);
    corrAngHist->GetYaxis()->SetTitle("Arb. Units");
    corrAngHist->GetYaxis()->SetTitleOffset(0.95);
    corrAngHist->GetYaxis()->SetTitleSize(0.05);
    corrAngHist->GetYaxis()->SetNoExponent(kTRUE);
    corrAngHist->SetStats(0);
    corrAngHist->SetLineWidth(2.0);
    corrAngHist->SetLineColor(kGreen+2);
    corrAngHist->Draw("HIST");
    corrAngHist->Draw("AXISsame");
    corrAngHist->SetMinimum(0.0001);
    //c_corr->SaveAs("corrAngHist.png");
    
    TCanvas *c_combined = new TCanvas();
    c_combined->cd();
    TLegend *leg_combined = new TLegend(0.55,0.70,0.88,0.85);
    leg_combined->SetLineColor(kWhite);
    leg_combined->AddEntry(origAngHist,"Before Calibration","L");
    leg_combined->AddEntry(corrAngHist,"After Calibration","L");
    origAngHist->GetXaxis()->SetTitle("#Delta#theta [degrees]");
    origAngHist->Draw("HIST");
    corrAngHist->Draw("HISTsame");
    std::cout <<"Cal: " << corrAngHist->GetMean() << std::endl;
    std::cout <<"Uncal: " << origAngHist->GetMean() << std::endl;
    std::cout << "NTracks: " << totalTracks << std::endl;
    leg_combined->Draw("same");
    origAngHist->Draw("AXISsame");
    origAngHist->SetMaximum(1.1*std::max(origAngHist->GetMaximum(),corrAngHist->GetMaximum()));
    origAngHist->SetMinimum(0.0001);
    c_combined->SaveAs(plotName.c_str());
    //c_combined->SaveAs("/nashome/j/joelam/combinedAngHist.pdf");
    
    //fileMeans->Write();
    //fileMeans->Close();
    return_vec.push_back(corrAngHist->GetMean());
    return_vec.push_back(origAngHist->GetMean());
    return_vec.push_back(totalTracks);
    
    return return_vec;
}

void eFieldCalculator::plotSpaceChargeBoundary(std::string inputMapFileName, std::string meanFileName){
    bool containsLaser = false;
    bool isMC = false;
    
    if(inputMapFileName.find("MC") != std::string::npos)
        isMC = true;
        if(inputMapFileName.find("Laser") != std::string::npos)
            containsLaser = true;
    
            std::cout << "Is MC: " << isMC << " Is laser: " << containsLaser << std::endl;
    
            //Bin along track segments.
            //Compare to the far edge.
    
    
    
            //Look at larger track slices
            //const double bufferLength = 0.05*100;
            double driftV = 0.001;
            //double driftV = 0.01; //for optimization studies
            if(isMC)
               driftV   = 0.0;
    
            const int minTrackPoints = 75;
            /*const int numTrackSegPoints = 15;
            const int lastSegment       = 2;
            const int firstSegment      = 0;
            const int maxTrackSegments    = 1;
            const int stepsBack           = 1;
            const int maxBadPoints        = 5;*/
    
            TGaxis::SetMaxDigits(3);
    

    
           /* const double xLow  = 0.0;
            const double xHigh = Lx;
            const double yLow = 0.0;
            const double yHigh = Ly;
    
            const double zLow = 0.0;
            const double zHigh = Lz;*/
    
            TFile *distortionMapInput = new TFile(inputMapFileName.c_str());
    
            TH3F* dist_dX  = (TH3F*) distortionMapInput->Get("Reco_Displacement_X");
            TH3F* dist_dY  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Y");
            TH3F* dist_dZ  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Z");
    
            TH3F* dist_dX_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_X_Error");
            TH3F* dist_dY_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_Y_Error");
            TH3F* dist_dZ_err = (TH3F*) distortionMapInput->Get("Reco_Displacement_Z_Error");
    
            std::string TrackFileName = "/uboone/data/users/mrmooney/ForJoel/output_data_tracks.root";
            if(isMC)
               TrackFileName = "/uboone/data/users/mrmooney/ForJoel/output_MC_tracks.root";
    
            TFile* inputFile = new TFile(TrackFileName.c_str());

    
    
            TTreeReader readerTracks("SpaCEtree_tracks", inputFile);
            TTreeReaderValue<Int_t> nElec_tracks(readerTracks, "nElec_tracks");
            TTreeReaderArray<Double_t> elecX_tracks_orig(readerTracks, "elecX_tracks");
            TTreeReaderArray<Double_t> elecY_tracks_orig(readerTracks, "elecY_tracks");
            TTreeReaderArray<Double_t> elecZ_tracks_orig(readerTracks, "elecZ_tracks");
    
            //std::cout << dist_dX->GetNbinsX() << " " << dist_dX->GetNbinsY() << " " << dist_dX->GetNbinsZ() << std::endl;
            TFile* fileMeans = new TFile(meanFileName.c_str(), "RECREATE");
    
            TH3F* h_tpcVolume = new TH3F("tpcVolume","",1000,0,10.0,250,0,2.5,250,0,2.5);
            TH3F* h_spaceChargeVolume = new TH3F("spaceChargeVolume","",1000,0,10.0,250,0,2.5,250,0,2.5);
    
            TH2F *h_tpcVolume_XY = new TH2F("tpcVolume_XY", "", 250,0,2.5,250,0,2.5);
            TH2F *h_spaceChargeVolume_XY = new TH2F("spaceChargeVolume_XY", "", 250,0,2.5,250,0,2.5);
    
            TH2F *h_tpcVolume_YZ = new TH2F("tpcVolume_YZ", "", 1000,0,10.0,250,0,2.5);
            TH2F *h_spaceChargeVolume_YZ = new TH2F("spaceChargeVolume_YZ", "", 1000,0,10.0,250,0,2.5);
    
            TH2F *h_tpcVolume_XZ = new TH2F("tpcVolume_XZ", "", 1000,0,10.0,250,0,2.5);
            TH2F *h_spaceChargeVolume_XZ = new TH2F("spaceChargeVolume_XZ", "", 1000,0,10.0,250,0,2.5);
    
            h_tpcVolume_XY->SetMarkerColor(kRed);
            h_tpcVolume_XY->SetLineColor(kRed);
            h_tpcVolume_XY->SetLineWidth(2.0);
            h_tpcVolume_XY->GetXaxis()->SetTitle("Distance from Cathode (m)");
            h_tpcVolume_XY->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
            h_spaceChargeVolume_XY->SetLineColor(kGreen+2);
            h_spaceChargeVolume_XY->SetMarkerColor(kGreen+2);
            h_spaceChargeVolume_XY->SetLineWidth(2.0);
            h_spaceChargeVolume_XY->GetXaxis()->SetTitle("Distance from Cathode (m)");
            h_spaceChargeVolume_XY->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
    
            h_tpcVolume_YZ->SetMarkerColor(kRed);
            h_tpcVolume_YZ->GetXaxis()->SetTitle("Distance from TPC Endcap (m)");
            h_tpcVolume_YZ->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
            h_spaceChargeVolume_YZ->SetMarkerColor(kGreen+2);
            h_spaceChargeVolume_YZ->GetXaxis()->SetTitle("Distance from TPC Endcap (m)");
            h_spaceChargeVolume_YZ->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
    
            h_tpcVolume_XZ->SetMarkerColor(kRed);
            h_tpcVolume_XZ->GetXaxis()->SetTitle("Distance from TPC Endcap (m)");
            h_tpcVolume_XZ->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
            h_spaceChargeVolume_XZ->SetMarkerColor(kGreen+2);
            h_spaceChargeVolume_XZ->GetXaxis()->SetTitle("Distance from TPC Endcap (m)");
            h_spaceChargeVolume_XZ->GetYaxis()->SetTitle("Distance from TPC Bottom (m)");
    
            TLegend *legend = new TLegend(0.1, 0.3, 0.5, 0.55);
            legend->AddEntry(h_tpcVolume_XY, "No Space Charge Correction", "l");
            legend->AddEntry(h_spaceChargeVolume_XY, "With Space Charge Correction", "l");
    
    
    
    
    
            int counter = 0;
            while ((readerTracks.Next()) && (counter < 213000)){
              //std::cout << "Go go go" << std::endl;
              if(counter % 10000 == 0) {std::cout << counter << std::endl;}
                counter++;
                
            
                if(*nElec_tracks >= minTrackPoints) {
                    
                    std::vector<double> elecX_tracks;
                    std::vector<double> elecY_tracks;
                    std::vector<double> elecZ_tracks;
                    
                    for(int i = 0; i < *nElec_tracks; i++){
                        //std::vector<double> points_laser;
                        //std::vector<double> points_cosmic;
                        double driftVFactor = cosmicDriftVScale;
                       
                        const double pointOffset = driftV*cosmicDriftVScale*(Lx-elecX_tracks_orig[i]);
                        
                        elecX_tracks.push_back(elecX_tracks_orig[i]  + pointOffset);
                        elecY_tracks.push_back(elecY_tracks_orig[i]);
                        elecZ_tracks.push_back(elecZ_tracks_orig[i]);
                        
                        int point = i;
                        
                        //Only interested in top / bottom of cathode
                        if(elecX_tracks[point] > 2.75 || (elecY_tracks[point] > 0.5 && elecY_tracks[point] < 2.0 ) || (elecZ_tracks[point] > 1.5 && elecZ_tracks[point] < 8.5 )  )
                            continue;
                        
                        h_tpcVolume->Fill(elecZ_tracks[point], elecX_tracks[point], elecY_tracks[point]);
                        h_tpcVolume_XY->Fill(elecX_tracks[point], elecY_tracks[point]);
                        h_tpcVolume_XZ->Fill(elecZ_tracks[point], elecX_tracks[point]);
                        h_tpcVolume_YZ->Fill(elecY_tracks[point], elecZ_tracks[point]);
                        
                        float xDistortion = -0.01*dist_dX->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                        float yDistortion = 0.01*dist_dY->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                        float zDistortion = 0.01*dist_dZ->Interpolate(doInvCoordTransformX(elecX_tracks[point]), doInvCoordTransformY(elecY_tracks[point]), doInvCoordTransformZ(elecZ_tracks[point]));
                        h_spaceChargeVolume->Fill((elecZ_tracks[point]+zDistortion), (elecX_tracks[point]+xDistortion), (elecY_tracks[point]+yDistortion));
                        
                        h_spaceChargeVolume_XY->Fill(elecX_tracks[point]+xDistortion, elecY_tracks[point]+yDistortion);
                        h_spaceChargeVolume_XZ->Fill(elecZ_tracks[point]+zDistortion, elecX_tracks[point]+xDistortion);
                        h_spaceChargeVolume_YZ->Fill(elecY_tracks[point]+yDistortion, elecZ_tracks[point]+zDistortion);
                       
                   }//end of loop over track points
                }//end of if nElec_tracks
            }//end of while loop
    gStyle->SetOptStat(0);
    TCanvas can("can", "", 600, 600);
    h_spaceChargeVolume_XY->Draw();
    h_tpcVolume_XY->Draw("same");
    
    legend->Draw("same");
    can.SaveAs("SpaceChargeVolumeXY.png");
    
    
    fileMeans->Write();
    fileMeans->Close();
    return;
}

void eFieldCalculator::compareDataMC(std::string inputData, std::string inputMC, std::string outputTag){

    TFile *data   = new TFile(inputData.c_str());
    TFile *mc     = new TFile(inputMC.c_str());
    TH1F *h1_data = (TH1F*) data->Get("h1");
    TH1F *h1_mc   = (TH1F*) mc->Get("h1");
    TH1F *h2_data = (TH1F*) data->Get("h2");
    TH1F *h2_mc   = (TH1F*) mc->Get("h2");
    
    TCanvas *c1 = new TCanvas("c1","D Map at Central Z",600,400);
    TLegend *leg_cal = new TLegend(0.55,0.70,0.88,0.85);
    double scale = 1.0;
    double sum = h2_mc->Integral();
    if(sum > 0)
        scale = (double) (h2_data->Integral() / sum);
    h2_mc->Scale(scale);
    h2_mc->GetXaxis()->SetRangeUser(0.0, 8.0);
    h2_mc->SetTitle("Residual of Track Hits After Correction");
    h2_mc->Draw("hist");
    h2_data->SetMarkerStyle(kFullDotLarge);
    h2_data->SetMarkerColor(kBlack);
    h2_data->Draw("E0 same");
    leg_cal->AddEntry(h2_data, "Data", "p");
    leg_cal->AddEntry(h2_mc, "MC", "l");
    leg_cal->Draw("same");
    c1->SaveAs(Form("DataMCResidualPlots/ResidualDataMCAfterCalib%s.png", outputTag.c_str()) );
    c1->Update();
    c1->Clear();
    
    TLegend *leg_uncal = new TLegend(0.55,0.70,0.88,0.85);
    scale = 1.0;
    sum = h1_mc->Integral();
    if(sum > 0)
        scale = (double) (h1_data->Integral() / sum);
    h1_mc->Scale(scale);
    h1_mc->GetXaxis()->SetRangeUser(0.0, 12.0);
    h1_mc->SetTitle("Residual of Track Hits Before Correction");
    h1_mc->Draw("hist");
    h1_data->SetMarkerStyle(kFullDotLarge);
    h1_data->SetMarkerColor(kBlack);
    h1_data->Draw("E0 same");
    leg_uncal->AddEntry(h1_data, "Data", "p");
    leg_uncal->AddEntry(h1_mc, "MC", "l");
    leg_uncal->Draw("same");
    c1->SaveAs(Form("DataMCResidualPlots/ResidualDataMCBeforeCalib%s.png", outputTag.c_str()));
    c1->Update();
    c1->Clear();
    
    
    

    

}

std::vector<double> eFieldCalculator::Residual_afterTrackCorr(std::string inputMapFileName, std::string plotName){
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.15);
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetPadTopMargin(0.12);
    gStyle->SetPadLeftMargin(0.11);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetLegendTextSize(0.04);
    
    std::vector<double> return_vec;
    
    const double epsilon = 1e-10;
    const int    pointReduction = 10;
    const int TailDivisions = 5;
    const int TailBin = 1;
    const double TailValue[TailDivisions] = {1.0, 2.0, 4.0, 6.0, 8.0};
    int badPoints_up = 0;
    int badPoints_down = 0;
    int oov_Hits = 0;
    double driftV = 0.001;
    bool doTriLin = false;
    
    bool containsLaser = false;
    bool isMC = false;
    const bool laserRegionOnly = false;
    
    
    
    if(inputMapFileName.find("MC") != std::string::npos)
        isMC = true;
    if(inputMapFileName.find("Laser") != std::string::npos)
        containsLaser = true;
    if(isMC)
        driftV = 0.0;
    else
        driftV = 0.015;
    
    
    std::cout << "Is MC: " << isMC << " Is laser: " << containsLaser << std::endl;
    std::cout << "Drift V: " << driftV << std::endl;
    
    gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h");
    gInterpreter->GenerateDictionary("vector<vector<vector<TVector3>>>","TVector3.h");
    
    // Input
    TChain *RecoTrackTree = new TChain("tracks");
    TChain *LaserInfoTree = new TChain("lasers");
    std::string inputTracks;
    std::string outputName;
    
    if(isMC){
        inputTracks = "/uboone/data/users/joelam/SCEInputFiles/lasersim-FullSim-lcs1.root";
        outputName  = "TrackResidualsLaserOnlyMC.root";
    }
    
    else{
        inputTracks = "/uboone/data/users/joelam/SCEInputFiles/laserbeams-data-newsmooth-7268.root";
        outputName  = "TrackResidualsLaserOnly.root";
    }
    
    
    RecoTrackTree->Add(inputTracks.c_str());
    LaserInfoTree->Add(inputTracks.c_str());
    
    // Laser truth Info
    TVector3 EntryPoint;
    TVector3 *pEntryPoint = &EntryPoint;
    TVector3 ExitPoint;
    TVector3 *pExitPoint = &ExitPoint;
    
    LaserInfoTree->SetBranchAddress("entry", &pEntryPoint);
    LaserInfoTree->SetBranchAddress("exit", &pExitPoint);
    
    // Track Info
    int EventNumber;
    std::vector<TVector3> TrackSamples;
    std::vector<TVector3>* pTrackSamples = &TrackSamples;
    
    RecoTrackTree->SetBranchAddress("track", &pTrackSamples);
    RecoTrackTree->SetBranchAddress("event", &EventNumber);
    
    // Read Distortion
    TFile *InFile = new TFile(inputMapFileName.c_str());
    //TFile *InFile = new TFile("output_hists_data_200k_Aug3.root","READ");
    //TFile *InFile = new TFile("/uboone/data/users/joelam/SCEInputFiles/MergedMapsCosmicAndLaser.root", "READ");
    //TFile *InFile = new TFile("/uboone/data/users/joelam/SCEInputFiles/MergedMapsCosmicOnly.root", "READ");
    //TFile *InFile = new TFile("/uboone/data/users/joelam/SCEInputFiles/MergedMapsLaserOnly.root", "READ");
    
    //TFile *InFile = new TFile("RecoCorr-N3-S50-Data-2side-Anode.root", "READ");
    
    TH3F *Dx = (TH3F*) InFile->Get("Reco_Displacement_X");
    TH3F *Dy = (TH3F*) InFile->Get("Reco_Displacement_Y");
    TH3F *Dz = (TH3F*) InFile->Get("Reco_Displacement_Z");
    
    /*
    TH3F *Dx = (TH3F*) InFile->Get("combined_dX");
    TH3F *Dy = (TH3F*) InFile->Get("combined_dY");
    TH3F *Dz = (TH3F*) InFile->Get("combined_dZ");
    */
    /*
    TH3F *Dx_err = (TH3F*) InFile->Get("combined_dX_Error");
    TH3F *Dy_err = (TH3F*) InFile->Get("combined_dY_Error");
    TH3F *Dz_err = (TH3F*) InFile->Get("combined_dZ_Error");
    */
    double corr_Dx[Dx->GetNbinsX()][Dx->GetNbinsY()][Dx->GetNbinsZ()];
    double corr_Dy[Dx->GetNbinsX()][Dx->GetNbinsY()][Dx->GetNbinsZ()];
    double corr_Dz[Dx->GetNbinsX()][Dx->GetNbinsY()][Dx->GetNbinsZ()];
    
    for (int x = 0; x <Dx->GetNbinsX(); x++ ) {
        for (int y = 0; y <Dx->GetNbinsY(); y++ ) {
            for (int z = 0; z <Dx->GetNbinsZ(); z++ ) {
                corr_Dx[x][y][z] = 0.0;
                corr_Dy[x][y][z] = 0.0;
                corr_Dz[x][y][z] = 0.0;
                
            }
        }
    }
    
    
    // Define number of tracks
    Int_t nElecX_tracks;
    int nElec_tracks;
    std::vector<double> elecX;
    std::vector<double> elecY;
    std::vector<double> elecZ;
    
    // Define histograms for residuals
    TCanvas *c1 = new TCanvas("c1","D Map at Central Z",600,400);
    TFile *OutFile = new TFile(outputName.c_str(), "RECREATE");
    TH1F *h1 = new TH1F("h1","Residuals of track hits before correction",100,-1,19);
    TH1F *h2 = new TH1F("h2","Residuals of track hits after correction",100,-1,19);
    TH1F *h_IntX = new TH1F("IntX", "Difference in Interpolation #Delta X", 250, -5.0, 5.0);
    TH1F *h_IntY = new TH1F("IntX", "Difference in Interpolation #Delta Y", 250, -5.0, 5.0);
    TH1F *h_IntZ = new TH1F("IntX", "Difference in Interpolation #Delta Z", 250, -5.0, 5.0);
    
    TH1F *h_trkX = new TH1F("trkX", "X Track Points", nCalibDivisions_x+1,xMin,xMax);
    TH1F *h_trkY = new TH1F("trkY", "Y Track Points", nCalibDivisions_y+1,yMin,yMax);
    TH1F *h_trkZ = new TH1F("trkZ", "Z Track Points", nCalibDivisions_z+1,zMin,zMax);
    
    TH2F *h_tracksXZ = new TH2F("tracksXZ", "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_x+1,xMin,xMax);
    TH2F *h_tracksZY = new TH2F("tracksZY", "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_y+1,yMin,yMax);
    
    TH2F *h_tracksXZInTail = new TH2F("tracksXZInTail", "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_x+1,xMin,xMax);
    TH2F *h_tracksZYInTail = new TH2F("tracksZYInTail", "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_y+1,yMin,yMax);
    
    TH2F *h_tracksXZ_tail[TailDivisions];
    TH2F *h_tracksZY_tail[TailDivisions];
    for(int i = 0; i < TailDivisions; ++i){
        h_tracksXZ_tail[i] = new TH2F(Form("tracksXZ_tail_%d"), "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_x+1,xMin,xMax);
        h_tracksZY_tail[i] = new TH2F(Form("tracksZY_tail_%d"), "Tracks XZ", 100*nCalibDivisions_z+1,zMin,zMax,100*nCalibDivisions_y+1,yMin,yMax);
        
    }
    
    
   
    
    TH3F *h_tracks   = new TH3F("tracks", "Tracks", nCalibDivisions_x+1,xMin,xMax, nCalibDivisions_y+1,yMin,yMax, nCalibDivisions_z+1,zMin,zMax);
    TH3F *h_tracks_tail   = new TH3F("tracks_tail", "Tracks", nCalibDivisions_x+1,xMin,xMax, nCalibDivisions_y+1,yMin,yMax, nCalibDivisions_z+1,zMin,zMax);
    
    
    TTree *T_tracks = new TTree("SpaCEtree_tracks","SpaCEtree_tracks");
    T_tracks->Branch("elecX_tracks",&elecX);
    T_tracks->Branch("elecY_tracks",&elecY);
    T_tracks->Branch("elecZ_tracks",&elecZ);
    T_tracks->Branch("nElec_tracks",&nElec_tracks,"nElec_tracks/I");
   // T_tracks->SetDirectory(OutFile);
    
    TF1 *f1 = new TF1("chi2L", "ROOT::Math::chisquared_pdf(x, [0])", 0, 10.0);
    
    f1->SetParameters(1,2);
    f1->SetLineColor(kBlack);
    f1->SetLineWidth(2);
    f1->SetLineStyle(9);
    
    // Loop over all tree entries
    if (LaserInfoTree->GetEntries() == RecoTrackTree->GetEntries()) {
        for (Size_t id = 0; id < RecoTrackTree->GetEntries(); id++) {
            
            RecoTrackTree->GetEntry(id);
            LaserInfoTree->GetEntry(id);
            
            TVector3 TrackVec = ExitPoint - EntryPoint;
            double LTrack = TrackVec.Mag();
            
            nElecX_tracks = TrackSamples.size();
            elecX.clear();
            elecY.clear();
            elecZ.clear();
            nElec_tracks= TMath::Floor(nElecX_tracks/pointReduction);
            
            for(int i = 0; i < nElecX_tracks; i++){
                
                // Pre-Correction
                TVector3 HitVec = TrackSamples[i] - EntryPoint;
                TVector3 CrossProduct = TrackVec.Cross(HitVec);
                h1->Fill(CrossProduct.Mag() / LTrack);
                
                double xPoint = TrackSamples[i][0];
                double yPoint = TrackSamples[i][1];
                double zPoint = TrackSamples[i][2];
                
                
                if(i % pointReduction == 0){
                  elecX.push_back(doCoordTransformX(xPoint));
                  elecY.push_back(doCoordTransformY(yPoint));
                  elecZ.push_back(doCoordTransformZ(zPoint));
                  h_trkX->Fill(xPoint);
                  h_trkY->Fill(yPoint);
                  h_trkZ->Fill(zPoint);
                }
                
                double xLength = xMax - xMin;
                double yLength = yMax - yMin;
                double zLength = zMax - zMin;
                
                //std::cout << (Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength) << " " << (Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)<< " " << (Int_t)TMath::Floor(zPoint*Dx->GetNbinsZ()/zLength) << std::endl;
                
                //std::cout << "Dx: " <<  corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dx->GetNbinsZ()/zLength)] << std::endl;
                
                float dXddd = corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXddu = corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXdud = corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXduu = corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXudd = corr_Dx[(Int_t)TMath::Ceil(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXudu = corr_Dx[(Int_t)TMath::Ceil(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXuud = corr_Dx[(Int_t)TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dx->GetNbinsZ()/zLength)];
                float dXuuu = corr_Dx[(Int_t)TMath::Ceil(xPoint*Dx->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dx->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dx->GetNbinsZ()/zLength)];
                
              //  std::cout << "Dy: " <<  corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dy->GetNbinsZ()/zLength)] << std::endl;
                
                float dYddd = corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYddu = corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYdud = corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYduu = corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYudd = corr_Dy[(Int_t)TMath::Ceil(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYudu = corr_Dy[(Int_t)TMath::Ceil(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYuud = corr_Dy[(Int_t)TMath::Floor(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dy->GetNbinsZ()/zLength)];
                float dYuuu = corr_Dy[(Int_t)TMath::Ceil(xPoint*Dy->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dy->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dy->GetNbinsZ()/zLength)];
                
               // std::cout << "Dz: " <<  corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength)] << std::endl;
                
                float dZddd = corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZddu = corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZdud = corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZduu = corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZudd = corr_Dz[(Int_t)TMath::Ceil(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZudu = corr_Dz[(Int_t)TMath::Ceil(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZuud = corr_Dz[(Int_t)TMath::Floor(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Floor(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)];
                float dZuuu = corr_Dz[(Int_t)TMath::Ceil(xPoint*Dz->GetNbinsX()/xLength)][(Int_t)TMath::Ceil(yPoint*Dz->GetNbinsY()/yLength)][(Int_t)TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)];
                
                float x_up   = TMath::Floor(xPoint*Dx->GetNbinsX()/xLength)*(Dx->GetNbinsX()/xLength)*100.0;
                float x_down = TMath::Ceil(xPoint*Dx->GetNbinsX()/xLength)*(Dx->GetNbinsX()/xLength)*100.0;
                float y_up   = TMath::Floor(yPoint*Dy->GetNbinsY()/yLength)*(Dy->GetNbinsY()/yLength)*100.0;
                float y_down = TMath::Ceil(yPoint*Dy->GetNbinsY()/yLength)*(Dy->GetNbinsY()/yLength)*100.0;
                float z_up   = TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength)*(Dz->GetNbinsZ()/zLength)*100.0;
                float z_down = TMath::Ceil(zPoint*Dz->GetNbinsZ()/zLength)*(Dz->GetNbinsZ()/zLength)*100.0;
                
                //std::cout << yPoint << " " << TMath::Ceil((yPoint)*Dy->GetNbinsY()/yLength) << " " <<TMath::Floor((yPoint)*Dy->GetNbinsY()/yLength) << " " << (Dy->GetNbinsY()/yLength) << std::endl;
                //std::cout << TMath::Floor(zPoint*Dz->GetNbinsZ()/zLength) << " " << (Dz->GetNbinsZ()/zLength) << std::endl;
                //std::cout << xPoint << " " << yPoint << " " << zPoint << std::endl;
               // std::cout << x_up << " " << x_down << " " << y_up << " " << y_down << " " << z_up << " " << z_down << std::endl;
                
                // Correction
                double distX = 0.0;
                double distY = 0.0;
                double distZ = 0.0;
                
                double distXRoot = Dx->Interpolate(TrackSamples[i][0], TrackSamples[i][1], TrackSamples[i][2]);
                double distYRoot = Dy->Interpolate(TrackSamples[i][0], TrackSamples[i][1], TrackSamples[i][2]);
                double distZRoot = Dz->Interpolate(TrackSamples[i][0], TrackSamples[i][1], TrackSamples[i][2]);
                
                double distXTriLin = TrilinInterp(xPoint, yPoint, zPoint, dXddd, dXddu, dXdud, dXduu, dXudd, dXudu, dXuud, dXuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                double distYTriLin = TrilinInterp(xPoint, yPoint, zPoint, dYddd, dYddu, dYdud, dYduu, dYudd, dYudu, dYuud, dYuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                double distZTriLin = TrilinInterp(xPoint, yPoint, zPoint, dZddd, dZddu, dZdud, dZduu, dZudd, dZudu, dZuud, dZuuu, x_up, x_down, y_up, y_down, z_up, z_down);
                
                if(doTriLin){
                    distX = distXTriLin - driftV*xPoint;
                    distY = distYTriLin;
                    distZ = distZTriLin;
                }
                
                else{
                   distX = distXRoot - driftV*TrackSamples[i][0];
                   distY = distYRoot;
                   distZ = distZRoot;
                }
                h_IntX->Fill(distXRoot - distXTriLin);
                h_IntY->Fill(distYRoot - distYTriLin);
                h_IntZ->Fill(distZRoot - distZTriLin);
                
               // std::cout << distXRoot << " " << distXTriLin << " " << distYRoot << " " << distYTriLin << " " << distZRoot << " " << distZTriLin << std::endl;
                
                // Define corrected hit
                TVector3 newHit;
                
                // If the correction doesn't exist, then skip this hit
                if(distX>1E10 || distY>1E10 || distZ>1E10){
                    ++badPoints_up;
                    continue;
                }
                
                //Only Look at hits in the center
                else if(laserRegionOnly && !isInLaserRegion(xPoint, yPoint, zPoint)){
                    ++oov_Hits;
                    continue;
                }
                
                // If the correction doesn't exist, then skip this hit
                else if(fabs(distX) < epsilon || fabs(distY) < epsilon || fabs(distZ) < epsilon){
                    ++badPoints_down;
                    continue;
                }
                
                

                
                // Otherwise, apply the correction
                else{
                    TVector3 dist(distX, distY, distZ);
                    newHit = TrackSamples[i] + dist;
                    
                }
                
                // After-Correction
                TVector3 newHitVec = newHit - EntryPoint;
                TVector3 newCrossProduct = TrackVec.Cross(newHitVec);
                double residual = (newCrossProduct.Mag() / LTrack);
                h2->Fill(residual);
                h_tracksXZ->Fill(TrackSamples[i][2],TrackSamples[i][0]);
                h_tracksZY->Fill(TrackSamples[i][2],TrackSamples[i][1]);
                h_tracks->Fill(TrackSamples[i][0], TrackSamples[i][1], TrackSamples[i][2]);
                
                if(residual > TailValue[TailBin]){
                    h_tracksXZInTail->Fill(TrackSamples[i][2],TrackSamples[i][0]);
                    h_tracksZYInTail->Fill(TrackSamples[i][2],TrackSamples[i][1]);
                }
                
                for(int bin=0; bin < TailDivisions; ++bin){
                    //std::cout << TailValue[i] << std::endl;
                     if(residual > TailValue[bin]){
                        h_tracksXZ_tail[bin]->Fill(TrackSamples[i][2],TrackSamples[i][0]);
                        h_tracksZY_tail[bin]->Fill(TrackSamples[i][2],TrackSamples[i][1]);
                     //h_tracks_tail->Fill(TrackSamples[i][0], TrackSamples[i][1], TrackSamples[i][2]);
                    
                }
              }
                
            } //end of loop over track pointss
            
            T_tracks->Fill();
            
        }
    }
    else{
        std::cout<<"The laser info and the track info do not match!"<<std::endl;
    }
    
    h1->SetLineColor(2);
    h1->SetLineWidth(2);
    
    h2->SetLineColor(4);
    h2->SetLineWidth(2);
    
    std::string Title;
    if(doTriLin)
        Title = "Residual of Track Hits before and after Cosmic Correction (Tri Linear Interpolation)";
    else
        Title = "Residual of Track Hits before and after Correction";
    h2->SetTitle(Title.c_str());
    h2->GetXaxis()->SetTitle("Distance to the true track [cm]");
    h2->GetYaxis()->SetTitle("Entries / 0.2 cm");
    //h2->GetYaxis()->SetRangeUser(0.0, 1.0);
    //h1->GetYaxis()->SetRangeUser(0.0, 1.0);
    //double scale = 1.0;
    /*if(h1->Integral() > 0)
        scale = 1 / (h1->Integral());
    h1->Scale(scale);
    if(h2->Integral() > 0)
        scale = 1 / (h2->Integral());*/
   //h2->Scale(scale);
    
    h2->SetMaximum(400e3);
    h2->Draw("hist");
    h1->Draw("hist SAME");
    //f1->Draw("SAME");
   // h2->GetYaxis()->SetRangeUser(0.0, 0.5);
   // h1->GetYaxis()->SetRangeUser(0.0, 0.5);
    
    TLegend *legend = new TLegend(0.4,0.65,0.8,0.85);
    legend->AddEntry(h1,"Residual before correction");
    legend->AddEntry(h2,"Residual after correction");
    legend->Draw("SAME");
    
    
    c1->Update();
    c1->SaveAs(plotName.c_str());
    
    std::cout << "Bad points (too big): " << badPoints_up << std::endl;
    std::cout << "Bad points (too small): " << badPoints_down << std::endl;
    std::cout << "Out of Fiducial Hits: " << oov_Hits << std::endl;
    std::cout << "Mean before correction: " <<  h1->GetMean() << std::endl;
    std::cout << "Mean after correction: " <<  h2->GetMean() << std::endl;
    
    return_vec.push_back(h2->GetMean());
    return_vec.push_back(h1->GetMean());
    return_vec.push_back(h1->Integral());
    
    c1->Clear();
    h_trkX->GetXaxis()->SetTitle("Track Point (cm)");
    h_trkX->GetYaxis()->SetTitle("Number of Points");
    h_trkX->Draw("hist");
    c1->Update();
    c1->SaveAs("TrackX.pdf");
    
    c1->Clear();
    h_trkY->GetXaxis()->SetTitle("Track Point (cm)");
    h_trkY->GetYaxis()->SetTitle("Number of Points");
    h_trkY->Draw("hist");
    c1->Update();
    c1->SaveAs("TrackY.pdf");
    
    c1->Clear();
    h_trkZ->GetXaxis()->SetTitle("Track Point (cm)");
    h_trkZ->GetYaxis()->SetTitle("Number of Points");
    h_trkZ->Draw("hist");
    c1->Update();
    c1->SaveAs("TrackZ.pdf");
    
    drawPlanarPlot(h_tracksXZ, 0, "Laser Tracks", "LaserTracksZX", axisType::yAxis);
    drawPlanarPlot(h_tracksZY, 0, "Laser Tracks", "LaserTracksZY", axisType::xAxis);
    
    drawPlanarPlot(h_tracksXZInTail, 0, Form("Laser Tracks Reisdual > %.2f cm", TailValue[TailBin]), Form("LaserTracksZXInTail"), axisType::yAxis);
    drawPlanarPlot(h_tracksZYInTail, 0, Form("Laser Tracks Reisdual > %.2f cm", TailValue[TailBin]), Form("LaserTracksZYInTail"), axisType::xAxis);
    for(int i=0; i < TailDivisions; ++i){
      drawPlanarPlot(h_tracksXZ_tail[i], 0, Form("Laser Tracks Reisdual > %.2f cm", TailValue[i]), Form("LaserTracksZXInTail%d", i), axisType::yAxis);
      drawPlanarPlot(h_tracksZY_tail[i], 0, Form("Laser Tracks Reisdual > %.2f cm", TailValue[i]), Form("LaserTracksZYInTail%d", i), axisType::xAxis);
    }
    OutFile->Write();
    OutFile->Close();
    

    return return_vec;
    
    
}

void eFieldCalculator::makeCSVMap(std::string inputFileName, std::string outputFileName){
    TFile *distortionMapInput = new TFile(inputFileName.c_str());
    
    std::ofstream outFileX;
    std::ofstream outFileY;
    std::ofstream outFileZ;
    
    std::string xMapName = outputFileName + ".csv";

    
    outFileX.open(xMapName.c_str());

    TH3F* dist_dX  = (TH3F*) distortionMapInput->Get("Reco_Displacement_X");
    TH3F* dist_dY  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Y");
    TH3F* dist_dZ  = (TH3F*) distortionMapInput->Get("Reco_Displacement_Z");
    
    const int startX = 1;
    const int endX   = dist_dX->GetNbinsX();
    const int startY = 1;
    const int endY   = dist_dX->GetNbinsY();
    const int startZ = 1;
    const int endZ   = dist_dX->GetNbinsZ();
    
    /*
    const int startX = 6;
    const int endX   = 21;
    const int startY = 8;
    const int endY   = 19;
    const int startZ = 21;
    const int endZ   = 71;
    */
    
    
    double x_reco;
    double y_reco;
    double z_reco;
    
    double xBins[dist_dX->GetNbinsX()];
    dist_dX->GetXaxis()->GetLowEdge(xBins);
    double yBins[dist_dX->GetNbinsY()];
    dist_dX->GetYaxis()->GetLowEdge(yBins);
    double zBins[dist_dX->GetNbinsZ()];
    dist_dX->GetZaxis()->GetLowEdge(zBins);
    
    for(Int_t i = startX; i <= endX; i++){
        if(i!=dist_dX->GetNbinsX())
            x_reco = (xBins[i-1] + xBins[i])/2.0;
        else
            x_reco = 256.0;
        for(Int_t j = startY; j <= endY; j++){
            if(j != dist_dX->GetNbinsY())
                y_reco = (yBins[j-1] + yBins[j])/2.0;
            else
                y_reco = 116.5;
            for(Int_t k = startZ; k <= endZ; k++){
                if(k != dist_dX->GetNbinsZ())
                    z_reco = (zBins[k-1] + zBins[k])/2.0;
                else
                    z_reco = 1037.0;
                
                
                
                outFileX << x_reco << ", " << y_reco << ", " << z_reco << ", " << dist_dX->GetBinContent(i, j, k) << ", " << dist_dY->GetBinContent(i, j, k)  << ", " << dist_dZ->GetBinContent(i, j, k)  << std::endl;
                //outFileY << x_reco << ", " << y_reco << ", " << z_reco << ", " << dist_dY->GetBinContent(i, j, k) << std::endl;
                //outFileZ << x_reco << ", " << y_reco << ", " << z_reco << ", " << dist_dZ->GetBinContent(i, j, k) << std::endl;
                
            }
        }
    }
}

float eFieldCalculator::SmoothBoundary(TH3 *h_cosmic, TH3 *h_laser, std::vector<int> firstBin, std::vector<int> secondBin, axisType axis){
    const double epsilon = 0.01;
    double startingWeight = 1.0;
    const int maxWeight = 600;
    double returnVal = 999.9;
    //double returnVal = 1.0;
    double binWidth = 10.0;
    /*
    if(axis == axisType::xAxis){
        binWidth = h_cosmic->GetXaxis()->GetBinWidth(firstBin[0]);
       
    }
    if(axis == axisType::yAxis)
            binWidth = h_cosmic->GetYaxis()->GetBinWidth(firstBin[1]);
    if(axis == axisType::zAxis){
            binWidth = h_cosmic->GetZaxis()->GetBinWidth(firstBin[2]);
        
    }*/
    //First, scan from 1.0 to max weight and see if we get some convergence.
    for(int i = 0; i <= maxWeight; ++i){
        double weight = (double) startingWeight - 0.01*i;
        double value = (weight*h_cosmic->GetBinContent(firstBin[0], firstBin[1], firstBin[2]) - h_laser->GetBinContent(secondBin[0], secondBin[1], secondBin[2]) ) / binWidth;
        if(0)
           std::cout << firstBin[0] << " " << firstBin[1] << " " <<firstBin[2] << " " << secondBin[0] << " " << secondBin[1] << " " << secondBin[2] << " " << weight << " " << value << std::endl;
        if(fabs(value) < epsilon)
            return weight;
    }
    
    for(int i = 0; i <= maxWeight; ++i){
        double weight = (double) startingWeight + 0.01*i;
        double value = (weight*h_cosmic->GetBinContent(firstBin[0], firstBin[1], firstBin[2]) - h_laser->GetBinContent(secondBin[0], secondBin[1], secondBin[2]) ) / binWidth;
        if(0)
           std::cout << firstBin[0] << " " << firstBin[1] << " " <<firstBin[2] << " " << secondBin[0] << " " << secondBin[1] << " " << secondBin[2] << " " << weight << " " << value << std::endl;
        if(fabs(value) < epsilon)
            return weight;
    }

    return returnVal;
}

float eFieldCalculator::TrilinInterp(TH3 *histo, std::vector<int> centerBin){
    //Create a dummy TH3 so I can use ROOTs TH3.
    TH3F *dummy = new TH3F("dummy", "", 4, 0.0, 4.0, 4, 0.0, 4.0, 4, 0.0, 4.0);
    int middleBin = 2;
    float returnVal = -999.0;
    //dummy->SetBinContent(middleBin,middleBin,middleBin, histo->GetBinContent(centerBin[0], centerBin[1], centerBin[2]));
    for(int i = -1; i <= 1; ++i){
        for(int j = -1; j <= 1; ++j){
            for(int k = -1; k<= 1; ++k){
                /*if(i == 0 && j == 0 && k == 0 )
                    continue;*/
                float something = histo->GetBinContent(centerBin[0]+i, centerBin[1]+j, centerBin[2]+k);
                dummy->SetBinContent(middleBin+i, middleBin+j, middleBin+k, something);
            }
        }
        
    }
    returnVal = dummy->Interpolate(2.0, 2.0, 2.0);
    //std::cout << returnVal << " " << histo->GetBinContent(centerBin[0]-1, centerBin[1], centerBin[2]) << " " << histo->GetBinContent(centerBin[0], centerBin[1]-1, centerBin[2]) << " " << histo->GetBinContent(centerBin[0], centerBin[1], centerBin[2]-1) << " " << histo->GetBinContent(centerBin[0]-1, centerBin[1]-1, centerBin[2]) << " " << histo->GetBinContent(centerBin[0]-1, centerBin[1], centerBin[2]-1) << " " << histo->GetBinContent(centerBin[0], centerBin[1]-1, centerBin[2]-1) << " " << histo->GetBinContent(centerBin[0]-1, centerBin[1]-1, centerBin[2]-1) << " " << histo->GetBinContent(centerBin[0]+1, centerBin[1], centerBin[2]) << " " << histo->GetBinContent(centerBin[0], centerBin[1]+1, centerBin[2]) << " " << histo->GetBinContent(centerBin[0], centerBin[1], centerBin[2]+1) << " " << histo->GetBinContent(centerBin[0]+1, centerBin[1]+1, centerBin[2]) << " " << histo->GetBinContent(centerBin[0]+1, centerBin[1], centerBin[2]+1) << " " << histo->GetBinContent(centerBin[0], centerBin[1]+1, centerBin[2]+1) << " " << histo->GetBinContent(centerBin[0]+1, centerBin[1]+1, centerBin[2]+1) << " " << histo->GetBinContent(centerBin[0], centerBin[1], centerBin[2]) <<  std::endl;
    
    delete dummy; //begone!
    return returnVal;
    
    
    
    
}

int main(int argc, char *argv[]){
  
  const bool doVoxIter = false;
  const bool doVelocityIter = false;
  eFieldCalculator *calculator = new eFieldCalculator();
  
    calculator->setDriftVScale(0.0, 0.0);
   // calculator->combineMaps(6, 21, 8, 19, 21, 79, true);
    
 // calculator->compareCalib(true);
//  calculator->compareCalibZXPlane(true);
  //calculator->compareTruth(false);
  //  calculator->combineWeightedMaps();
  //  calculator->combineMaps(true, true, false);
 //   calculator->combineMaps(true, false, false);
  //  calculator->combineMaps(true, false, true);
   // calculator->combineMaps(6, 21, 8, 19, 21, 79, true);
  //  calculator->studyResults2("MergedMapsCosmicOnlySmoothed.root", "AnglePlots/MergedMapsCosmicOnlySmoothed.png");
  //  calculator->studyResults2("MergedMapsCosmicAndLaserNoDriftV.root", "AnglePlots/MergedMapsCosmicAndLaserNoDriftV.png");
   // calculator->studyResults2("MergedMapsSmoothCosmicAndLaserV1098NoAvg.root", "AnglePlots/MergedMapsSmoothCosmicAndLaserV1098NoAvg.png");
   // calculator->studyResults2("MergedMapsLaserOnly.root", "AnglePlots/MergedMapsLaserOnly.png");
  //  calculator->studyResults2("MergedMapsCosmicOnlySmoothed.root", "AnglePlots/MergedMapsCosmicOnlySmoothed.png");
  //  calculator->studyResults2("MergedMapsCosmicOnlySmoothed.root", "AnglePlots/Dokus.png");
  //  calculator->studyResults2("/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.root", "AnglePlots/Dokus.png");
    calculator->studyResults2("/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsSmoothCosmicAndLaserV1098VolumeSmoothedMC.root", "AnglePlots/BadMC.png");
  //  calculator->plotSpaceChargeBoundary("/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.root", "SCESurface.root");
   // calculator->Residual_afterTrackCorr("MergedMapsSmoothCosmicAndLaserNoDriftVNoAvg.root", "ResidualPlots/MergedMapsSmoothCosmicAndLaserNoDriftVNoAvg.png");
   // calculator->Residual_afterTrackCorr("/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.root", "ResidualPlots/Dokus.png");
  //  calculator->Residual_afterTrackCorr("MergedMapsCosmicAndLaserNoDriftV.root", "ResidualPlots/MergedMapsCosmicAndLaserNoDriftVRes.png");
   // calculator->Residual_afterTrackCorr("MergedMapsSmoothCosmicAndLaserV1098VolumeSmoothed.root", "ResidualPlots/MergedMapsSmoothCosmicAndLaserV1098VolumeSmoothed.png");
    
  //  calculator->MakeDistortionHistograms(true);
  //  calculator->MakeDistortionHistograms(false);
  //  calculator->makeFwdMapPlots();
  //  calculator->makeSmoothMapPlots();
  //  calculator->makeSmoothMap("MergedMapsSmoothCosmicAndLaserNoDriftV.root", "EdgesSmoothed.root", false, true);
    //calculator->makeSmoothMap("EdgesSmoothed.root", "MergedMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.root", true, false);
  //  calculator->makeCSVMap("/uboone/data/users/joelam/SCEDistortionMaps/MergedMapsCosmicOnlySmoothed.root","MergedMapsCosmicOnlySmoothed");
   
    
    
  if(doVoxIter){
    const int maxVoxels   = 26;
    const int lowerZVoxel = 21;
    const int upperZVoxel = 79;
    const int largestVox = 10;
    int iteration = 0;
    TFile *file = new TFile("OptimizationHistograms.root", "RECREATE");
    TH2F *h_voxelsVAngle = new TH2F("voxelsVAngle", "Difference in Mean Angle (degrees)", maxVoxels, 0.0, maxVoxels, maxVoxels, 0.0, maxVoxels);
    TH2F *h_voxelsVResidual = new TH2F("voxelsVResidual", "Difference in Mean Residual (cm)", maxVoxels, 0.0, maxVoxels, maxVoxels, 0.0, maxVoxels);
    std::vector<double> angleMeans;
    std::vector<double> residualMeans;
    for(int xVox = 5; xVox < maxVoxels; ++xVox){
        
	int lowXVox = xVox;
        int highXVox = maxVoxels - xVox;
        if(lowXVox > largestVox) break;
        for(int yVox = 5; yVox < maxVoxels; ++yVox){
          std::cout << "Iteration: " << iteration << std::endl;
	  int lowYVox = yVox;
          int highYVox = maxVoxels - yVox;
          if(lowYVox > largestVox) break;
          
          calculator->combineMaps(lowXVox, highXVox, lowYVox, highYVox, lowerZVoxel, upperZVoxel, true);
          angleMeans = calculator->studyResults2();
          //std::cout << angleMeans[0] << " " << angleMeans[1] << std::endl;
          residualMeans = calculator->Residual_afterTrackCorr();
          //std::cout << residualMeans[0] << " " << residualMeans[1] << std::endl;
          h_voxelsVAngle->Fill(lowXVox, lowYVox, (angleMeans[0] - angleMeans[1]));
          h_voxelsVResidual->Fill(lowXVox, lowYVox, (residualMeans[0] - residualMeans[1]) );
	  ++iteration;
        }
        
    }
    
    
    file->Write();
    file->Close();
      delete file;
    }
    if(doVelocityIter){
        const double driftV = 0.01;
        const int maxVelocity   = 4;
        int iteration = 0;
        TFile *file = new TFile("VelocityOptimizationHistograms.root", "RECREATE");
        TH2F *h_voxelsVAngle = new TH2F("voxelsVAngle", "Difference in Mean Angle (degrees)", 7,  -0.03, 0.04, 7, -0.03, 0.04);
        TH2F *h_voxelsVResidual = new TH2F("voxelsVResidual", "Difference in Mean Residual (cm)", 7, -0.03, 0.04, 7, -0.03, 0.04);
        h_voxelsVAngle->GetXaxis()->SetTitle("Laser Drift V");
        h_voxelsVAngle->GetYaxis()->SetTitle("Cosmic Drift V");
        h_voxelsVResidual->GetXaxis()->SetTitle("Laser Drift V");
        h_voxelsVResidual->GetYaxis()->SetTitle("Cosmic Drift V");
        std::vector<double> angleMeans;
        std::vector<double> residualMeans;
        for(int xVox = -3; xVox <maxVelocity; ++xVox){
            for(int yVox = -3; yVox < maxVelocity; ++yVox){
                std::cout << "Iteration: " << iteration << std::endl;
                calculator->setDriftVScale(xVox, yVox); //puts laser drift V on the X axis
                //calculator->combineMaps(6, 21, 8, 19, 21, 79, true);
                calculator->combineMaps(6, 21, 8, 19, 1, 100, true);
                //angleMeans = calculator->studyResults2();
               
                residualMeans = calculator->Residual_afterTrackCorr();
                //h_voxelsVAngle->SetBinContent(xVox+4, yVox+4, (angleMeans[0] - angleMeans[1]));
                h_voxelsVResidual->SetBinContent(xVox+4, yVox+4, (residualMeans[0] - residualMeans[1]) );
                ++iteration;
            }
            
        }
        
        
        file->Write();
        file->Close();
        delete file;
    }
    
  return 0;
} //end of main

