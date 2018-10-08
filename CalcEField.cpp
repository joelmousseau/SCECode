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


void eFieldCalculator::compareCalib(bool isData)
{
  
  const double xScale = TPC_X/Lx;
  const double yScale = TPC_Y/Ly;
  const double zScale = TPC_Z/Lz;
  const int   zRegions = 3; //upstream, middle, downstream
  int       zRegion = 0;
  const double zMaximum = 3.0;

  double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};
  int   zCuts[2] = {35, 65};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);
  
  std::string inputLaser;
  std::string inputCosmic;
  if(isData){
     inputLaser = "RecoCorr-N3-S50-Data-2side-Anode.root";
     inputCosmic = "output_hists_data_200k_Aug3.root";
     
  }
  
  else{
     inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "output_hists_MC_200k_Aug3.root";
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
  
  bool plotX = false;
  bool plotY = false;
  bool plotZ = false;
  double driftV = 0.0;
  if(isData)
    driftV = 0.01;

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
	
	
	double x_correction = i*driftV*diff_dX->GetXaxis()->GetBinWidth(i);
	/*if(j== 1 && k == 1)
	  std::cout << i*diff_dX->GetXaxis()->GetBinWidth(i) << std::endl;*/
	//x_correction = 0.0;
	
	cosmic_Ds.push_back(cosmic_dX->GetBinContent(i,j,k)+x_correction);
	cosmic_Ds.push_back(cosmic_dY->GetBinContent(i,j,k));
	cosmic_Ds.push_back(cosmic_dZ->GetBinContent(i,j,k));
	
	laser_Ds.push_back(laser_dX->GetBinContent(i,j,k));
	laser_Ds.push_back(laser_dY->GetBinContent(i,j,k));
	laser_Ds.push_back(laser_dZ->GetBinContent(i,j,k));
	
	double cosine = getAngle(cosmic_Ds, laser_Ds);
			
	if(plotX)
	{
          diff_dX->SetBinContent(i,j,k,(laser_dX->GetBinContent(i,j,k)-cosmic_dX->GetBinContent(i,j,k)+x_correction) );
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
        
	plotX = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dY->GetBinContent(i,j,k), cosmic_dY_err->GetBinContent(i,j,k)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && goodCosmic(cosmic_dZ->GetBinContent(i,j,k), cosmic_dZ_err->GetBinContent(i,j,k)) );
	
	
	if(plotX){	
	  laser_2D_dX.SetBinContent(i,j,laser_dX->GetBinContent(i,j,k));
	  cosmic_2D_dX.SetBinContent(i,j,cosmic_dX->GetBinContent(i,j,k));
	  
	  h_diffXByRegion[zRegion]->Fill(diff_dX->GetBinContent(i,j,k));
	  if(zRegion == 1 && fabs(diff_dX->GetBinContent(i,j,k)) <= maxDiff)
	    good_dX->SetBinContent(i,j,k,1.0);
	}
	
	if(plotY){
	  
          laser_2D_dY.SetBinContent(i,j,laser_dY->GetBinContent(i,j,k));
	  cosmic_2D_dY.SetBinContent(i,j,cosmic_dY->GetBinContent(i,j,k));
	  h_diffYByRegion[zRegion]->Fill(diff_dY->GetBinContent(i,j,k));
	  if(zRegion == 1 && fabs(diff_dY->GetBinContent(i,j,k)) <= maxDiff)
	    good_dY->SetBinContent(i,j,k,1.0);
	}
	
	if(plotZ){
	  laser_2D_dZ.SetBinContent(i,j,laser_dZ->GetBinContent(i,j,k));
	  cosmic_2D_dZ.SetBinContent(i,j,cosmic_dZ->GetBinContent(i,j,k));
	  h_diffZByRegion[zRegion]->Fill(diff_dZ->GetBinContent(i,j,k));
	  if(zRegion == 1 && fabs(diff_dZ->GetBinContent(i,j,k)) <= maxDiff)
	    good_dZ->SetBinContent(i,j,k,1.0);
	
	
	}
	
	if(plotX && plotY && plotZ){
	   cosine_2D.SetBinContent(i,j, cosine_angles->GetBinContent(i,j,k));
	}  
	
        diff_2D_dX.SetBinContent(i,j,diff_dX->GetBinContent(i,j,k));
        diff_2D_dY.SetBinContent(i,j,diff_dY->GetBinContent(i,j,k));
        diff_2D_dZ.SetBinContent(i,j,diff_dZ->GetBinContent(i,j,k));
	  
      }
    }

    
    
   drawPlanarPlot(laser_2D_dX, k, "Laser #DeltaX", "laser_2D_dX", axisType::zAxis);
   drawPlanarPlot(laser_2D_dY, k, "Laser #DeltaY", "laser_2D_dY", axisType::zAxis);
   drawPlanarPlot(laser_2D_dZ, k, "Laser #DeltaZ", "laser_2D_dZ", axisType::zAxis);
   drawPlanarPlot(cosmic_2D_dX, k, "Cosmic #DeltaX", "cosmic_2D_dX", axisType::zAxis);
   drawPlanarPlot(cosmic_2D_dY, k, "Cosmic #DeltaY", "cosmic_2D_dY", axisType::zAxis);
   drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic #DeltaZ", "cosmic_2D_dZ", axisType::zAxis);
   drawPlanarPlot(diff_2D_dX, k, "Laser-Cosmic #DeltaX", "diff_2D_dX", axisType::zAxis, zMaximum);
   drawPlanarPlot(diff_2D_dY, k, "Laser-Cosmic #DeltaY", "diff_2D_dY", axisType::zAxis, zMaximum);
   drawPlanarPlot(diff_2D_dZ, k, "Laser-Cosmic #DeltaZ", "diff_2D_dZ", axisType::zAxis, zMaximum);
   drawPlanarPlot(cosine_2D,  k, "Cosine of Laser - Cosmic", "cos_2D", axisType::zAxis, 1.0);
      
   

    

    
   

    
    
    
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
  c_1D_diff_dX.SaveAs("plotdump/diff_1D_dX.png");
  
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
  c_1D_diff_dY.SaveAs("plotdump/diff_1D_dY.png");
  
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
  c_1D_diff_dZ.SaveAs("plotdump/diff_1D_dZ.png");
  
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

  double stops[5] = {0.00,0.34,0.61,0.84,1.00};
  double red[5] = {0.00,0.00,0.87,1.00,0.51};
  double green[5] = {0.00,0.81,1.00,0.20,0.00};
  double blue[5] = {0.51,1.00,0.12,0.00,0.00};
  int   zCuts[2] = {9, 17};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);
  
  std::string inputLaser;
  std::string inputCosmic;
  if(isData){
     inputLaser = "RecoCorr-N3-S50-Data-2side-Anode.root";
     inputCosmic = "output_hists_data_200k_Aug3.root";
  }
  
  else{
     inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "output_hists_MC_200k_Aug3.root";
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
  
  bool plotX = false;
  bool plotY = false;
  bool plotZ = false;
  double driftV = 0.0;
  if(isData)
    driftV = 0.01;

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
          
	double x_correction = i*driftV*diff_dX->GetXaxis()->GetBinWidth(i);
	/*if(j== 1 && k == 1)
	  std::cout << i*diff_dX->GetXaxis()->GetBinWidth(i) << std::endl;*/
	//x_correction = 0.0;
    
     cosmic_Ds.push_back(cosmic_dX->GetBinContent(i,j,k)+x_correction);
     cosmic_Ds.push_back(cosmic_dY->GetBinContent(i,j,k));
     cosmic_Ds.push_back(cosmic_dZ->GetBinContent(i,j,k));
          
     laser_Ds.push_back(laser_dX->GetBinContent(i,j,k));
     laser_Ds.push_back(laser_dY->GetBinContent(i,j,k));
     laser_Ds.push_back(laser_dZ->GetBinContent(i,j,k));
          
          double cosine = getAngle(cosmic_Ds, laser_Ds);
			
	if(plotX)
	{
          diff_dX->SetBinContent(i,j,k,(laser_dX->GetBinContent(i,j,k)-cosmic_dX->GetBinContent(i,j,k)+x_correction) );
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
        
	plotX = (goodLaser(laser_dX->GetBinContent(i,k,j), laser_dX_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dX->GetBinContent(i,k,j), cosmic_dX_err->GetBinContent(i,k,j)) );
	plotY = (goodLaser(laser_dY->GetBinContent(i,k,j), laser_dY_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dY->GetBinContent(i,k,j), cosmic_dY_err->GetBinContent(i,k,j)) );
	plotZ = (goodLaser(laser_dZ->GetBinContent(i,k,j), laser_dZ_err->GetBinContent(i,k,j)) && goodCosmic(cosmic_dZ->GetBinContent(i,k,j), cosmic_dZ_err->GetBinContent(i,k,j)) );
	
	
	if(plotX){	
	  laser_2D_dX.SetBinContent(j,i,laser_dX->GetBinContent(i,k,j));
	  cosmic_2D_dX.SetBinContent(j,i,cosmic_dX->GetBinContent(i,k,j));
	  diff_2D_dX.SetBinContent(j,i,diff_dX->GetBinContent(i,k,j));
	  h_diffXByRegion[zRegion]->Fill(diff_dX->GetBinContent(i,k,j));
	}
	
	if(plotY){
	  
          laser_2D_dY.SetBinContent(j,i,laser_dY->GetBinContent(i,k,j));
	  cosmic_2D_dY.SetBinContent(j,i,cosmic_dY->GetBinContent(i,k,j));
	  diff_2D_dY.SetBinContent(j,i,diff_dY->GetBinContent(i,k,j));
	  h_diffYByRegion[zRegion]->Fill(diff_dY->GetBinContent(i,k,j));
	}
	
	if(plotZ){
	  laser_2D_dZ.SetBinContent(j,i,laser_dZ->GetBinContent(i,k,j));
	  cosmic_2D_dZ.SetBinContent(j,i,cosmic_dZ->GetBinContent(i,k,j));
	  diff_2D_dZ.SetBinContent(j,i,diff_dZ->GetBinContent(i,k,j));
	  h_diffZByRegion[zRegion]->Fill(diff_dZ->GetBinContent(i,k,j));
	
	
	}
          
    if(plotX && plotY && plotZ){
       cosine_2D.SetBinContent(j,i,cosine_angles->GetBinContent(i,k,j));
    }
		  
      }
    }

    drawPlanarPlot(laser_2D_dX, k, "Laser #DeltaX", "laser_2D_dX", axisType::yAxis);
    drawPlanarPlot(laser_2D_dY, k, "Laser #DeltaY", "laser_2D_dY", axisType::yAxis);
    drawPlanarPlot(laser_2D_dZ, k, "Laser #DeltaZ", "laser_2D_dZ", axisType::yAxis);
    drawPlanarPlot(cosmic_2D_dX, k, "Cosmic #DeltaX", "laser_2D_dX", axisType::yAxis);
    drawPlanarPlot(cosmic_2D_dY, k, "Cosmic #DeltaY", "laser_2D_dY", axisType::yAxis);
    drawPlanarPlot(cosmic_2D_dZ, k, "Cosmic #DeltaZ", "laser_2D_dZ", axisType::yAxis);
    drawPlanarPlot(diff_2D_dX, k, "Laser-Cosmic #DeltaX", "diff_2D_dX", axisType::yAxis, zMaximum);
    drawPlanarPlot(diff_2D_dY, k, "Laser-Cosmic #DeltaY", "diff_2D_dY", axisType::yAxis, zMaximum);
    drawPlanarPlot(diff_2D_dZ, k, "Laser-Cosmic #DeltaZ", "diff_2D_dZ", axisType::yAxis, zMaximum);
    drawPlanarPlot(cosine_2D,  k, "Cosine of Laser - Cosmic", "cos_2D", axisType::yAxis, 1.0);
    
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
  c_1D_diff_dX.SaveAs("plotdump/diff_1D_dX.png");
  
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
  c_1D_diff_dY.SaveAs("plotdump/diff_1D_dY.png");
  
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
  c_1D_diff_dZ.SaveAs("plotdump/diff_1D_dZ.png");   
  
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
 /* TFile* fileCosmic = new TFile("output_hists_MC_200k_Aug3.root");
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("Reco_Displacement_X");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("Reco_Displacement_Y");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("Reco_Displacement_Z");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("Reco_Displacement_X_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Y_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("Reco_Displacement_Z_Error");*/
    
  TFile* fileCosmic = new TFile("MergedMaps.root");
  TH3F* cosmic_dX = (TH3F*) fileCosmic->Get("combined_dX");
  TH3F* cosmic_dY = (TH3F*) fileCosmic->Get("combined_dY");
  TH3F* cosmic_dZ = (TH3F*) fileCosmic->Get("combined_dZ");
  TH3F* cosmic_dX_err = (TH3F*) fileCosmic->Get("combined_dX_Error");
  TH3F* cosmic_dY_err = (TH3F*) fileCosmic->Get("combined_dY_Error");
  TH3F* cosmic_dZ_err = (TH3F*) fileCosmic->Get("combined_dZ_Error");
    
    
  
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

    TCanvas c_laser_dX(Form("canv_2D_dX_%d",k),"",600,600);
    c_laser_dX.cd();
    laser_2D_dX.SetTitle(Form("Laser - Truth #DeltaX [cm]:  ^{}Z_{reco} = %.2f cm",(((double) k-1)/100.0)*1036.8));
    laser_2D_dX.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dX.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dX.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dX.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dX.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dX.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dX.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dX.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dX.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    laser_2D_dX.SetStats(0);
    laser_2D_dX.Draw("COLZ");
    c_laser_dX.SaveAs(Form("Truth_plots/laser_2D_dX_%d.png",k));

    TCanvas c_laser_dY(Form("canv_2D_dY_%d",k),"",600,600);
    c_laser_dY.cd();
    laser_2D_dY.SetTitle(Form("Laser - Truth #DeltaY [cm]:  ^{}Z_{reco} = %.2f cm",(((double) k-1)/100.0)*1036.8));
    laser_2D_dY.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dY.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dY.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dY.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dY.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dY.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dY.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dY.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dY.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    laser_2D_dY.SetStats(0);
    laser_2D_dY.Draw("COLZ");
    c_laser_dY.SaveAs(Form("Truth_plots/laser_2D_dY_%d.png",k));
    
    TCanvas c_laser_dZ(Form("canv_2D_dZ_%d",k),"",600,600);
    c_laser_dZ.cd();
    laser_2D_dZ.SetTitle(Form("Laser - Truth #DeltaZ [cm]:  ^{}Z_{reco} = %.2f cm",(((double) k-1)/100.0)*1036.8));
    laser_2D_dZ.GetXaxis()->SetTitle("X_{reco} [cm]");
    laser_2D_dZ.GetXaxis()->SetTitleOffset(1.0);
    laser_2D_dZ.GetXaxis()->SetTitleSize(0.04);
    laser_2D_dZ.GetYaxis()->SetTitle("Y_{reco} [cm]");
    laser_2D_dZ.GetYaxis()->SetTitleOffset(1.2);
    laser_2D_dZ.GetYaxis()->SetTitleSize(0.04);
    laser_2D_dZ.GetZaxis()->SetNoExponent(kTRUE);
    laser_2D_dZ.GetZaxis()->SetLabelSize(0.025);
    laser_2D_dZ.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    laser_2D_dZ.SetStats(0);
    laser_2D_dZ.Draw("COLZ");
    c_laser_dZ.SaveAs(Form("Truth_plots/laser_2D_dZ_%d.png",k));

    TCanvas c_cosmic_dX(Form("canv_2D_dX_%d",k),"",600,600);
    c_cosmic_dX.cd();
    cosmic_2D_dX.SetTitle(Form("Cosmic - Truth #DeltaX [cm]:  ^{}Z_{reco} = %.2f cm",(((Double_t) k-1)/100.0)*1036.8));
    cosmic_2D_dX.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dX.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dX.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dX.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dX.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dX.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dX.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dX.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dX.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    cosmic_2D_dX.SetStats(0);
    cosmic_2D_dX.Draw("COLZ");
    c_cosmic_dX.SaveAs(Form("Truth_plots/cosmic_2D_dX_%d.png",k));

    TCanvas c_cosmic_dY(Form("canv_2D_dY_%d",k),"",600,600);
    c_cosmic_dY.cd();
    cosmic_2D_dY.SetTitle(Form("Cosmic - Truth #DeltaY [cm]:  ^{}Z_{reco} = %.2f cm",(((double) k-1)/100.0)*1036.8));
    cosmic_2D_dY.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dY.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dY.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dY.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dY.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dY.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dY.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dY.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dY.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    cosmic_2D_dY.SetStats(0);
    cosmic_2D_dY.Draw("COLZ");
    c_cosmic_dY.SaveAs(Form("Truth_plots/cosmic_2D_dY_%d.png",k));
    
    TCanvas c_cosmic_dZ(Form("canv_2D_dZ_%d",k),"",600,600);
    c_cosmic_dZ.cd();
    cosmic_2D_dZ.SetTitle(Form("Cosmic - Truth #DeltaZ [cm]:  ^{}Z_{reco} = %.2f cm",(((double) k-1)/100.0)*1036.8));
    cosmic_2D_dZ.GetXaxis()->SetTitle("X_{reco} [cm]");
    cosmic_2D_dZ.GetXaxis()->SetTitleOffset(1.0);
    cosmic_2D_dZ.GetXaxis()->SetTitleSize(0.04);
    cosmic_2D_dZ.GetYaxis()->SetTitle("Y_{reco} [cm]");
    cosmic_2D_dZ.GetYaxis()->SetTitleOffset(1.2);
    cosmic_2D_dZ.GetYaxis()->SetTitleSize(0.04);
    cosmic_2D_dZ.GetZaxis()->SetNoExponent(kTRUE);
    cosmic_2D_dZ.GetZaxis()->SetLabelSize(0.025);
    cosmic_2D_dZ.GetZaxis()->SetRangeUser(-zMaximum, zMaximum);
    cosmic_2D_dZ.SetStats(0);
    cosmic_2D_dZ.Draw("COLZ");
    c_cosmic_dZ.SaveAs(Form("Truth_plots/cosmic_2D_dZ_%d.png",k));

    
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
      sprintf(plotTitle, "%s [cm]:  ^{}Y_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_Y) - Ly/2.0) );
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
        sprintf(plotTitle, "%s [cm]:  ^{}Z_{reco} = %.2f cm", label,(((double) planeNum-1)/100.0)*TPC_Z);
    }
    
    //ZX Plane
    if(axis == axisType::yAxis){
        canWidth = 900;
        xAxisLabel = "Z_{reco} [cm]";
        yAxisLabel = "X_{reco} [cm]";
        dir        = "2DPlots_ZX";
        sprintf(plotTitle, "%s [cm]:  ^{}Y_{reco} = %.2f cm", label,(((((double) planeNum-1)/25.0)*TPC_Y) - Ly/2.0) );
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
    
    histOne->Scale(1.0/histOne->Integral());
    std::cout << histOne->GetMaximum() << std::endl;
    histTwo->Scale(1.0/histTwo->Integral());
    
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
            
            //bool  plotLaser  =  goodLaser(laser_dY->GetBinContent(i,25,j), laser_dY_err->GetBinContent(i,25,j));
            //bool  plotCosmic =  goodCosmic(cosmic_dY->GetBinContent(i,25,j), cosmic_dY_err->GetBinContent(i,25,j));
            bool plotLaser  = true;
            bool plotCosmic = true;
            /*
            if(plotLaser && plotCosmic)
              std::cout << laserVal << " " << cosmicVal << " " << fabs(face_dY->GetBinContent(i,j)) << std::endl;
            */
            if(plotLaser)
              laser_2D_dY.SetBinContent(j,i,  (laserVal  - fabs(face_dY->GetBinContent(i,j) ) ) );
            if(plotCosmic)
              cosmic_2D_dY.SetBinContent(j,i, (cosmicVal - fabs(face_dY->GetBinContent(i,j) ) ) );
            
            
            
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
            double laserVal  =   laser_dX->Interpolate(xVal, yVal, zVal);
            double cosmicVal =  cosmic_dX->Interpolate(xVal, yVal, zVal);
            
            //bool  plotLaser  =  goodLaser(laser_dY->GetBinContent(i,25,j), laser_dY_err->GetBinContent(i,25,j));
            //bool  plotCosmic =  goodCosmic(cosmic_dY->GetBinContent(i,25,j), cosmic_dY_err->GetBinContent(i,25,j));
            bool plotLaser  = true;
            bool plotCosmic = true;
            
            if(plotLaser && plotCosmic && j == 24)
                 std::cout << zVal << ", " << yVal << ": " <<  laserVal << " " << cosmicVal << " " <<  -1.0*face_dX->GetBinContent(i,j) << std::endl;
             
            if(plotLaser)
               laser_2D_dX.SetBinContent(j,i,  (laserVal  - face_dX->GetBinContent(i,j) ) ) ;
            if(plotCosmic)
               cosmic_2D_dX.SetBinContent(j,i, (cosmicVal - face_dX->GetBinContent(i,j) ) ) ;
            
            
            
        }
    }
    
    drawPlanarPlot(laser_2D_dY, 26, "Laser - Face #DeltaY (top)", "laser_2D_top_dY", axisType::yAxis, 6.0);
    drawPlanarPlot(cosmic_2D_dY, 26, "Cosmic - Face #DeltaY (top)", "cosmic_2D_top_dY", axisType::yAxis, 6.0);
    
    drawPlanarPlot(laser_2D_dX, 26,  "Laser - Face #DeltaX (cathode)", "laser_2D_cathode_dX", axisType::xAxis, 15.0);
    drawPlanarPlot(cosmic_2D_dX, 26, "Cosmic - Face #DeltaX (cathode)", "cosmic_2D_cathode_dX", axisType::xAxis, 5.0);
}

void eFieldCalculator::combineMaps(bool isData, bool skipLaser){
  const int lowZ    = 21;
  const int highZ   = 87;
  const int lowX    = 0;
  const int highX   = 40;
  const int lowY    = 4;
  const int highY   = 21;
  double driftV = 0.0;
  //TF1 *f1 = new TF1("chi2", calcChi2, 0, 100.0, 2);
  //Minimizer *min = new Minimizer("Minuit", "Migrad");
    
    if(isData)
      driftV = 0.01;
  
  
  std::string inputLaser;
  std::string inputCosmic;
  std::string inputGoodVoxels;
  std::string inputTruth;
  std::string outputName;
  if(skipLaser)
      outputName = "MergedMapsCosmicOnly.root";
  else
      outputName = "MergedMapsCosmicAndLaser.root";
  if(isData){
     inputLaser = "RecoCorr-N3-S50-Data-2side-Anode.root";
     inputCosmic = "output_hists_data_200k_Aug3.root";
     inputGoodVoxels = "GoodVoxels.root";
     inputTruth      = "output_truth_hists.root";
  }
  
  else{
     inputLaser = "RecoCorr-N3-S50-LaserMC-2side-Anode.root";
     inputCosmic = "output_hists_MC_200k_Aug3.root";
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
    
  TFile* fileGoodVoxels = new TFile(inputGoodVoxels.c_str());
  TH3F*  good_dX = (TH3F*) fileGoodVoxels->Get("good_dX");
  TH3F*  good_dY = (TH3F*) fileGoodVoxels->Get("good_dY");
  TH3F*  good_dZ = (TH3F*) fileGoodVoxels->Get("good_dZ");
  
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
	 
	 useCosmic = goodCosmic(cosmic_dX->GetBinContent(i,j,k), cosmic_dX_err->GetBinContent(i,j,k));
	 useLaser  = (goodLaser(laser_dX->GetBinContent(i,j,k), laser_dX_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
     goodAgreement = true;
     
     if(!goodAgreement){
         numerator = 0.0;
         denominator = 1.0;
         error = 0.0;
     }
	 else if(useCosmic && useLaser){
       numerator = ( (cosmic_dX->GetBinContent(i,j,k)+x_correction)/cosmic_dX_err->GetBinContent(i,j,k) + laser_dX->GetBinContent(i,j,k)/laser_dX_err->GetBinContent(i,j,k) );
	   denominator = 1/cosmic_dX_err->GetBinContent(i,j,k) + 1/laser_dX_err->GetBinContent(i,j,k);
	   error = sqrt(pow(cosmic_dX_err->GetBinContent(i,j,k), 2) + pow(laser_dX_err->GetBinContent(i,j,k), 2));
         double weights[2] = {0.5, 0.5};
         double vals[4]    = {cosmic_dX->GetBinContent(i,j,k)+x_correction, laser_dX->GetBinContent(i,j,k), truth_dX->GetBinContent(i,j,k), error};
         float chi2 = calcChi2(weights, vals);
         if(chi2 != 0) chi2_plot->Fill(chi2);
         else std::cout << "Chi2 equals 0!" << std::endl;
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
	 useLaser  = (goodLaser(laser_dY->GetBinContent(i,j,k), laser_dY_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
     goodAgreement = true;
     
    if(!goodAgreement){
      numerator = 0.0;
      denominator = 1.0;
      error = 0.0;
    }
	 
	else if(useCosmic && useLaser){
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
	 useLaser  = (goodLaser(laser_dZ->GetBinContent(i,j,k), laser_dZ_err->GetBinContent(i,j,k)) && k >= lowZ && k < highZ && j >= lowY && j < highY && !skipLaser);
     //useLaser = false;
     goodAgreement = true;
     
     if(!goodAgreement){
              numerator = 0.0;
              denominator = 1.0;
              error = 0.0;
     }
	 
	 else if(useCosmic && useLaser){
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

float eFieldCalculator::LinInterp(float x, float x1, float x2, float q00, float q01) {
    return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
}

float eFieldCalculator::calcChi2(double *scale, double *vals){
    float numerator = pow(((scale[0]*vals[0]+scale[1]*vals[1]) - vals[2]),2);
    float denominator = pow(vals[3],2);
    
    if(denominator > 0.0)
        return numerator/denominator;
    else
        return 0.0;
    
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

void eFieldCalculator::studyResults2(bool skipLaser)
{
    //Bin along track segments.
    //Compare to the far edge.
    
    //Look at larger track slices
    //const double bufferLength = 0.05*100;
    double driftV = 0.01;
    //double driftV   = 0.0;
    const double bufferLength = 0.0;
    
    const int minTrackPoints = 50;
    const int numTrackSegPoints = 15;
    
    TGaxis::SetMaxDigits(3);
    
    double numDivisions_x = nCalibDivisions_x;
    double numDivisions_y = nCalibDivisions_y;
    double numDivisions_z = nCalibDivisions_z;
    
    //const double zMax = Lz;
    //const double zMin = 0.0;
    
    const double zLow = 2.1;
    const double zHigh = 9.0;
    
    
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
    
    //std::string inputMapFileName = "MergedMapsCosmicAndLaser.root";
    std::string inputMapFileName;
    std::string plotName;
    std::string meanFileName;
    
    if(skipLaser){
       inputMapFileName = "MergedMapsCosmicOnly.root";
       plotName = "AnglePlots/CombinedAngCosmicOnly.png";
       meanFileName = "MeansCosmicOnly.root";
    }
    
    else{
        std::cout << "Using Laser..." << std::endl;
        inputMapFileName = "MergedMapsCosmicAndLaser.root";
        plotName = "AnglePlots/CombinedAngCosmicAndLaser.png";
        meanFileName = "MeansCosmicAndLaser.root";
    }
    //TFile* inputFileInterp = new TFile("/uboone/data/users/mrmooney/ForJoel/output_MC_30k.root");

    
    //Change this to TH3F input
    /*TTreeReader readerCalib("SpaCEtree_calib", inputFileInterp);
    TTreeReaderValue<Double_t> Dx_calib(readerCalib, "Dx.data_calib");
    TTreeReaderValue<Double_t> Dy_calib(readerCalib, "Dy.data_calib");
    TTreeReaderValue<Double_t> Dz_calib(readerCalib, "Dz.data_calib");
    TTreeReaderValue<Double_t> xReco_calib(readerCalib, "x_reco.data_calib");
    TTreeReaderValue<Double_t> yReco_calib(readerCalib, "y_reco.data_calib");
    TTreeReaderValue<Double_t> zReco_calib(readerCalib, "z_reco.data_calib");
    TTreeReaderValue<Double_t> xTrue_calib(readerCalib, "x_true.data_calib");
    TTreeReaderValue<Double_t> yTrue_calib(readerCalib, "y_true.data_calib");
    TTreeReaderValue<Double_t> zTrue_calib(readerCalib, "z_true.data_calib");
    TTreeReaderValue<Int_t> elecFate_calib(readerCalib, "elecFate.data_calib");*/
    
    TFile *distortionMapInput = new TFile(inputMapFileName.c_str());
    TH3F* dist_dX = (TH3F*) distortionMapInput->Get("combined_dX");
    TH3F* dist_dY = (TH3F*) distortionMapInput->Get("combined_dY");
    TH3F* dist_dZ = (TH3F*) distortionMapInput->Get("combined_dZ");
    
    TH3F* dist_dX_err = (TH3F*) distortionMapInput->Get("combined_dX_Error");
    TH3F* dist_dY_err = (TH3F*) distortionMapInput->Get("combined_dY_Error");
    TH3F* dist_dZ_err = (TH3F*) distortionMapInput->Get("combined_dZ_Error");
    
    double corr_Dx[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    double corr_Dy[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    double corr_Dz[dist_dX->GetNbinsX()][dist_dX->GetNbinsY()][dist_dX->GetNbinsZ()];
    
    double origAngMean[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    double origAngNum[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    
    double corrAngMean[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    double corrAngNum[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
    
    
    
    //TFile* inputFile = new TFile("output.root");
    TFile* inputFile = new TFile("/uboone/data/users/mrmooney/ForJoel/output_data_tracks.root");
    ////TFile* inputFile = new TFile("output_MC_NoSCE_tracks.root");
    //TFile* inputFile = new TFile("output_data_tracks.root");
    ////TFile* inputFile = new TFile("output_data_tracks_new.root");
    ////TFile* inputFile = new TFile("output_data_tracks_noMCScut.root");
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
                if(!goodLaser(dist_dX->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1), dist_dX_err->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1)))
                    corr_Dx[i][j][k] = 0.0;
                corr_Dy[i][j][k] = 0.01*dist_dY->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1);
                if(!goodLaser(dist_dY->GetBinContent(dist_dY->GetNbinsX()-i,j+1,k+1), dist_dY_err->GetBinContent(dist_dY->GetNbinsX()-i,j+1,k+1)))
                    corr_Dy[i][j][k] = 0.0;
                corr_Dz[i][j][k] = 0.01*dist_dZ->GetBinContent(dist_dX->GetNbinsX()-i,j+1,k+1);
                if(!goodLaser(dist_dZ->GetBinContent(dist_dZ->GetNbinsX()-i,j+1,k+1), dist_dZ_err->GetBinContent(dist_dZ->GetNbinsX()-i,j+1,k+1)))
                    corr_Dz[i][j][k] = 0.0;
                //std::cout << corr_Dz[i][j][k] << std::endl;
        /*corr_Dx[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dx_calib;
        corr_Dy[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dy_calib;
        corr_Dz[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dz_calib;*/
            }
        }
    }
    
    int counter = 0;
    //while (readerTracks.Next())
    //while ((readerTracks.Next()) && (counter < 100000))
    while ((readerTracks.Next()) && (counter < 213000))
    {
        if(counter % 1000 == 0) {std::cout << counter << std::endl;}
        counter++;
        
        if (*nElec_tracks >= minTrackPoints) {
            
            PointCloud startPoints;
            PointCloud endPoints;
            
            std::vector<double> elecX_tracks;
            std::vector<double> elecY_tracks;
            std::vector<double> elecZ_tracks;
            
            int xBin = TMath::Nint(elecX_tracks_orig[0]*numDivisions_x/Lx);
            int yBin = TMath::Nint(elecY_tracks_orig[0]*numDivisions_y/Ly);
            int zBin = TMath::Nint(elecZ_tracks_orig[0]*numDivisions_z/Lz);
            
            if(xBin > nCalibDivisions_x)
                xBin = nCalibDivisions_x;
            if(yBin > nCalibDivisions_y)
                yBin = nCalibDivisions_y;
            if(zBin > nCalibDivisions_z)
                zBin = nCalibDivisions_z;
            
            for(int i = 0; i < *nElec_tracks; i++){
                std::vector<double> points_laser;
                //std::vector<double> points_cosmic;
                
                const double pointOffset = driftV*(Lx-elecX_tracks_orig[i]);
                
                points_laser.push_back(elecX_tracks_orig[i]);
                points_laser.push_back(elecY_tracks_orig[i]);
                points_laser.push_back(elecZ_tracks_orig[i]);
                
                //points_cosmic = cosmicToLaser(points_laser);
                //std::cout << points_laser[0] << " " << pointOffset << std::endl;
                elecX_tracks.push_back(points_laser[0] + pointOffset);
                elecY_tracks.push_back(points_laser[1] );
                elecZ_tracks.push_back(points_laser[2] );
                
                
            }
            
            
            Int_t numBadPoints_start = 0;
            for (int i = 0; i < numTrackSegPoints; i++) {
                if ((elecX_tracks[i] <= 0.0-bufferLength) || (elecX_tracks[i] >= Lx+bufferLength) || (elecY_tracks[i] <= (0 - bufferLength)) || (elecY_tracks[i] >= Ly+bufferLength) || (elecZ_tracks[i] <= 0-bufferLength) || (elecZ_tracks[i] >= Lz+bufferLength)) {
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
            if (numBadPoints_start > 5) continue;
            
            Int_t numBadPoints_end = 0;
            for (int i = *nElec_tracks-1; i > *nElec_tracks-1-numTrackSegPoints; i--) {
                if ((elecX_tracks[i] <= 0.0-bufferLength) || (elecX_tracks[i] >= Lx+bufferLength) || (elecY_tracks[i] <=(0 - bufferLength)) || (elecY_tracks[i] >= Ly+bufferLength) || (elecZ_tracks[i] <= 0-bufferLength) || (elecZ_tracks[i] >= Lz+bufferLength)) {
                    numBadPoints_end++;
                    continue;
                }
                
                Point tempPoint;
                tempPoint.x = elecX_tracks[i];
                tempPoint.y = elecY_tracks[i];
                tempPoint.z = elecZ_tracks[i];
                
                endPoints.push_back(tempPoint);
            }
            if (numBadPoints_end > 5) continue;
            
            PCAResults results_start = DoPCA(startPoints);
            PCAResults results_end = DoPCA(endPoints);
            
            Double_t dotProd = results_start.eVecs[0](0)*results_end.eVecs[0](0) + results_start.eVecs[0](1)*results_end.eVecs[0](1) + results_start.eVecs[0](2)*results_end.eVecs[0](2);
            Double_t startMag = sqrt(pow(results_start.eVecs[0](0),2) + pow(results_start.eVecs[0](1),2) + pow(results_start.eVecs[0](2),2));
            Double_t endMag = sqrt(pow(results_end.eVecs[0](0),2) + pow(results_end.eVecs[0](1),2) + pow(results_end.eVecs[0](2),2));
            
            Double_t dTheta = TMath::ACos(dotProd/(startMag*endMag));
            
            origAngHist->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByX[xBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByY[yBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngHistByZ[zBin]->Fill(std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));
            origAngMean[xBin][yBin][zBin] += std::min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265));
            origAngNum[xBin][yBin][zBin]  += 1.0;
            
            PointCloud startPointsCorr;
            PointCloud endPointsCorr;
            
            Int_t numBadPoints_start_corr = 0;
            for (int i = 0; i < numTrackSegPoints; i++) {
                if ((elecX_tracks[i] <= 0.0) || (elecX_tracks[i] >= Lx) || (elecY_tracks[i] <= 0) || (elecY_tracks[i] >= Ly ) || (elecZ_tracks[i] <= 0) || (elecZ_tracks[i] >= Lz)) {
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
                tempPoint.x = elecX_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                tempPoint.y = elecY_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                tempPoint.z = elecZ_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                
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
            if (numBadPoints_start_corr > 5) continue;
            
            Int_t numBadPoints_end_corr = 0;
            for (int i = *nElec_tracks-1; i > *nElec_tracks-1-numTrackSegPoints; i--) {
                if ((elecX_tracks[i] <= 0.0) || (elecX_tracks[i] >= Lx) || (elecY_tracks[i] <= 0) || (elecY_tracks[i] >= Ly) || (elecZ_tracks[i] <= 0) || (elecZ_tracks[i] >= Lz) ) {
                    numBadPoints_end_corr++;
                    continue;
                }
                
                /*
                 cout << counter << " " << i << "   ";
                 cout << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
                 << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << endl;
                 */
                
                Point tempPoint;
                tempPoint.x = elecX_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                tempPoint.y = elecY_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                tempPoint.z = elecZ_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
                
                endPointsCorr.push_back(tempPoint);
            }
            if (numBadPoints_end_corr > 5) continue;
            
            PCAResults results_start_corr = DoPCA(startPointsCorr);
            PCAResults results_end_corr = DoPCA(endPointsCorr);
            
            Double_t dotProdCorr = results_start_corr.eVecs[0](0)*results_end_corr.eVecs[0](0) + results_start_corr.eVecs[0](1)*results_end_corr.eVecs[0](1) + results_start_corr.eVecs[0](2)*results_end_corr.eVecs[0](2);
            Double_t startMagCorr = sqrt(pow(results_start_corr.eVecs[0](0),2) + pow(results_start_corr.eVecs[0](1),2) + pow(results_start_corr.eVecs[0](2),2));
            Double_t endMagCorr = sqrt(pow(results_end_corr.eVecs[0](0),2) + pow(results_end_corr.eVecs[0](1),2) + pow(results_end_corr.eVecs[0](2),2));
            
            Double_t dThetaCorr = TMath::ACos(dotProdCorr/(startMagCorr*endMagCorr));
            
            //cout << "    ANGLE:  " << min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265));
            
            corrAngHist->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByX[xBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByY[yBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            corrAngHistByZ[zBin]->Fill(std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
            
            corrAngMean[xBin][yBin][zBin] += std::min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265));
            corrAngNum[xBin][yBin][zBin]  += 1.0;
           
        }
    }
    origAngHist->Scale(1.0/origAngHist->Integral());
    corrAngHist->Scale(1.0/corrAngHist->Integral());
    
    gStyle->SetTitleW(0.9);
    gStyle->SetOptStat(0);
    
    TFile* fileMeans = new TFile(meanFileName.c_str(), "RECREATE");
    
    TH3F* meansOrig = new TH3F("meansOrig","",dist_dX->GetNbinsX(),0,25.0,dist_dX->GetNbinsY(),0,25.0,dist_dX->GetNbinsZ(),0,100.0);
    TH3F* meansCorr = new TH3F("meansCorr","",dist_dX->GetNbinsX(),0,25.0,dist_dX->GetNbinsY(),0,25.0,dist_dX->GetNbinsZ(),0,100.0);
    TH3F* meansDiff = new TH3F("meansDiff","",dist_dX->GetNbinsX(),0,25.0,dist_dX->GetNbinsY(),0,25.0,dist_dX->GetNbinsZ(),0,100.0);
    
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
                double meanDiff = origMean - corrMean;
                meansDiff->SetBinContent(x+1, y+1, z+1, meanDiff);
                
                
            }
        }
    }
    
    for(int i=0; i < numDivisions_x; ++i){
        draw1DPlot(origAngHistByX[i], corrAngHistByX[i], i, "Angle", "AnglePlots/combinedAngHistInX", axisType::xAxis);
        drawPlanarPlot(PointsYZ[i], i, "Points", "Points", axisType::xAxis);
    }
    
    for(int i=0; i < numDivisions_y; ++i){
        draw1DPlot(origAngHistByY[i], corrAngHistByY[i], i, "Angle", "AnglePlots/combinedAngHistInY", axisType::yAxis);
    }
    
    for(int i=0; i < numDivisions_z; ++i){
        draw1DPlot(origAngHistByZ[i], corrAngHistByZ[i], i, "Angle", "AnglePlots/combinedAngHistInZ", axisType::zAxis);
    }
    
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
    std::cout << corrAngHist->GetMean() << std::endl;
    leg_combined->Draw("same");
    origAngHist->Draw("AXISsame");
    origAngHist->SetMaximum(1.1*std::max(origAngHist->GetMaximum(),corrAngHist->GetMaximum()));
    origAngHist->SetMinimum(0.0001);
    c_combined->SaveAs(plotName.c_str());
    //c_combined->SaveAs("/nashome/j/joelam/combinedAngHist.pdf");
    
    fileMeans->Write();
    fileMeans->Close();
}

void eFieldCalculator::compareMeans(){
    TFile *meansCosmicOnly = new TFile("MeansCosmicOnly.root");
    TH3F* diff_cosmic = (TH3F*) meansCosmicOnly->Get("meansDiff");
    
    TFile *meansCosmicAndLaser = new TFile("MeansCosmicAndLaser.root");
    TH3F* diff_cosmicAndLaser = (TH3F*) meansCosmicAndLaser->Get("meansDiff");
    
    TFile *file3D = new TFile("goodStartPoints.root", "RECREATE");
    TH3F *hist_3D = new TH3F("hist_3D","",diff_cosmic->GetNbinsX(),0,25.0,diff_cosmic->GetNbinsY(),0,25.0,diff_cosmic->GetNbinsZ(),0,100.0);
    
    for(int x = 1; x < diff_cosmic->GetNbinsX()+1; ++x){
        for(int y = 1; y < diff_cosmic->GetNbinsY()+1; ++y){
            for(int z = 1; z < diff_cosmic->GetNbinsZ()+1; ++z){
               double diffCosmic         = diff_cosmic->GetBinContent(x,y,z);
               double diffCosmicAndLaser = diff_cosmicAndLaser->GetBinContent(x,y,z);
               if(diffCosmicAndLaser > diffCosmic)
                   hist_3D->SetBinContent(x,y,z,1.0);
            }
       }
    
    
   }
    file3D->Write();
    file3D->Close();

}

int main(int argc, char *argv[]){

  eFieldCalculator *calculator = new eFieldCalculator();
  //calculator->compareCalib(true);
 // calculator->compareCalibZXPlane(true);
  //calculator->compareTruth(false);
  //calculator->combineMaps(true, true);
  //  calculator->compareFaces(true);
    calculator->studyResults2(true);
  //  calculator->compareMeans();

  return 0;
} //end of main

