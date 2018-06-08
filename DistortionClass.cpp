#include "DistortionClass.h"


distortionVoxel::distortionVoxel(){}

distortionVoxel::distortionVoxel(int x, int y, int z){
  xVoxel = x;
  yVoxel = y;
  zVoxel = z;

}

distortionVoxel::~distortionVoxel(){}

std::vector<int> distortionVoxel::getVoxelNumber() const{
   std::vector<int> returnVec;
   returnVec.push_back(xVoxel);
   returnVec.push_back(yVoxel);
   returnVec.push_back(zVoxel);
   
   return returnVec;


}

int distortionVoxel::setVoxelNumber(std::vector<int> index){
   if(index.size() != 3){
     std::cerr << "ERROR: Must specify a size 3 vector for setVoxelNumber!" <<std::endl;
     return 1;}
   
   xVoxel = index[0];
   yVoxel = index[1];
   zVoxel = index[2];
    
   //all is well
   return 0;

} 

int distortionVoxel::addDistortion(float dist){
   distortions.push_back(dist);
   return 0;

}

int distortionVoxel::clearDistortions(){
   distortions.clear();
   return 0;

}

int distortionVoxel::addWeight(float dist){
   weights.push_back(dist);
   return 0;

}

int distortionVoxel::clearWeights(){
   weights.clear();
   return 0;

}

std::vector<float> distortionVoxel::getDistortions() const{
   return distortions;

}

std::vector<double> distortionVoxel::getWeights() const{
   return weights;

}

unsigned int distortionVoxel::getNPairs() const{
   return distortions.size();
}

float distortionVoxel::calculateMean() const{
    float mean = 0.0;
    if(distortions.size() == 0 || weights.size() == 0)
       return mean;

    float wSum = calculateWeightSum();    

    if( wSum < 1e-6 && wSum > -1e6)
       return mean; 
           
    mean = TMath::Mean(distortions.begin(), distortions.end(), weights.begin() );
    return mean;

}

float distortionVoxel::calculateRMS() const{
    float rms = 0.0;
    
    if(distortions.size() == 0 || weights.size() == 0)
       return rms;
    
    float wSum = calculateWeightSum();    

    if( wSum < 1e-6 && wSum > -1e6)
       return rms;
    
    rms = TMath::RMS(distortions.begin(), distortions.end(), weights.begin() );
    return rms;

}

float distortionVoxel::calculateMedian() const{
    float median = 0.0;
    if(distortions.size() == 0 || weights.size() == 0)
       return median;
    float wSum = calculateWeightSum();    

    if( wSum < 1e-6 && wSum > -1e6)
       return median;   
    int pairs = getNPairs();
    
    median = TMath::Median(pairs, &distortions[0], &weights[0] );
    return median;
}

float distortionVoxel::calculateWeightSum() const{
    float wSum = 0.0;    
    for(unsigned int i = 0; i != weights.size(); ++i){
       wSum += weights[i];
    
    }
    return wSum;
}

distortionMap::distortionMap(){
   //intialize this beast
   for(int x=0; x < nXVoxels; ++x){
      std::vector< std::vector<distortionVoxel> > rowY;
      for(int y=0; y <nYVoxels; ++y){
         std::vector<distortionVoxel> rowZ;
	 for(int z=0; z < nZVoxels; ++z){  
            rowZ.push_back( distortionVoxel(x, y, z) );
         }
	 rowY.push_back(rowZ);
      }
      theMap.push_back(rowY);
   }

}

distortionMap::~distortionMap(){}

std::vector<  std::vector < std::vector<distortionVoxel> > > distortionMap::getMap(){
  return theMap;

}

distortionVoxel * distortionMap::getVoxel(int x, int y, int z){
  distortionVoxel *voxelPointer = &(theMap[x][y][z]);
  return voxelPointer;

}



std::vector< std::vector < std::vector<float> > > distortionMap::calculateMap(){
    float cv = 0.0;
    std::vector< std::vector < std::vector<float> > > return_vec;

    
    for(int x=0; x < nXVoxels; ++x){
       std::vector< std::vector<float> > rowY;
       for(int y=0; y < nYVoxels; ++y){
          std::vector<float> rowZ;
          for(int z=0; z < nZVoxels; ++z){
	     if(statMode == "mean")
	       cv = theMap[x][y][z].calculateMean();
	     else if(statMode == "median")
	       cv = theMap[x][y][z].calculateMedian();
	     else
	       cv = 0.0;  
	     rowZ.push_back( cv );	     	  
	  }
	  rowY.push_back(rowZ);
       }
       return_vec.push_back(rowY);
    }
    return return_vec;   	  
}

int distortionMap::FaceToVoxel(int face) const{
   //TOP
   if(face == 0)
     return nYVoxels;
   //BOTTOM
   else if(face == 1)
     return 0; 
   //UPSTREAM
   else if(face == 2)
     return 0;   
   //DOWNSTREAM
   else if(face == 3)
    return nZVoxels;
   //CATHODE IS AT 0
   else if(face == 4)  
    return 0;
   else{
     std::cout << "ERROR! Unkown TPC facing: " << face << std::endl;
     return -1;
   }  
      
}

std::vector < std::vector<float> > distortionMap::calculateMap(int face){
    float cv = 0.0;
    std::vector < std::vector<float> > return_vec;
    
    int fixedVox = FaceToVoxel(face);
    
    if(face == 0 || face == 1){         
       for(int x=0; x < nXVoxels; ++x){
         std::vector<float> rowZ;
         for(int z=0; z < nZVoxels; ++z){
	     if(statMode == "mean")
	       cv = theMap[x][fixedVox][z].calculateMean();
	     else if(statMode == "median")
	       cv = theMap[x][fixedVox][z].calculateMedian();
	     else
	       cv = 0.0;  
	     rowZ.push_back( cv );	     	  
	  }
	  return_vec.push_back(rowZ);
       }
       return return_vec;
    }
    
    else if(face == 2 || face == 3){         
       for(int x=0; x < nXVoxels; ++x){
         std::vector<float> rowY;
         for(int y=0; y < nYVoxels; ++y){
	     if(statMode == "mean")
	       cv = theMap[x][y][fixedVox].calculateMean();
	     else if(statMode == "median")
	       cv = theMap[x][y][fixedVox].calculateMedian();
	     else
	       cv = 0.0;  
	     rowY.push_back( cv );	     	  
	  }
	  return_vec.push_back(rowY);
       }
       return return_vec;
    }
    
    else if(face == 4){         
       for(int y=0; y < nYVoxels; ++y){
        std::vector<float> rowY;
         for(int z=0; z < nZVoxels; ++z){
	     if(statMode == "mean")
	       cv = theMap[fixedVox][y][z].calculateMean();
	     else if(statMode == "median")
	       cv = theMap[fixedVox][y][z].calculateMedian();
	     else
	       cv = 0.0;  
	     rowY.push_back( cv );	     	  
	  }
	  return_vec.push_back(rowY);
       }
       return return_vec;
    }
    
    else{
       std::cout << "ERROR! Unkown TPC facing: " << face << std::endl;
       return return_vec;    
    }
              	  
}

std::vector< std::vector < std::vector<float> > > distortionMap::calculateMapErrors(){
    float cv  = 0.0;
    float rms = 0.0;
    std::vector< std::vector < std::vector<float> > > return_vec; 
    
    for(int x=0; x < nXVoxels; ++x){
       std::vector< std::vector<float> > rowY;
       for(int y=0; y < nYVoxels; ++y){
          std::vector<float> rowZ;
          for(int z=0; z < nZVoxels; ++z){	     
	     cv = sqrt(theMap[x][y][z].getNPairs());
	     rms = theMap[x][y][z].calculateRMS();
	     if(cv != 0.0)
	        rms = (float)(rms / cv);
	     else
	        rms = 0.0;	
	       
	     rowZ.push_back( rms );	     	  
	  }
	  rowY.push_back(rowZ);
       }
       return_vec.push_back(rowY);
    }
    return return_vec;   	  
}

std::vector < std::vector<float> > distortionMap::calculateMapErrors(int face){
    float cv  = 0.0;
    float rms = 0.0;
    std::vector < std::vector<float> > return_vec;
    
    int fixedVox = FaceToVoxel(face);
    
    if(face == 0 || face == 1){         
       for(int x=0; x < nXVoxels; ++x){
         std::vector<float> rowZ;
         for(int z=0; z < nZVoxels; ++z){
	     cv  = sqrt(theMap[x][fixedVox][z].getNPairs() );
	     rms = theMap[x][fixedVox][z].calculateRMS();
	     if(cv != 0.0)
	       rms = (float)(rms / cv );
	     else
	       rms = 0.0;    
	     rowZ.push_back( rms );	     	  
	  }
	  return_vec.push_back(rowZ);
       }
       return return_vec;
    }
    
    else if(face == 2 || face == 3){         
       for(int x=0; x < nXVoxels; ++x){
         std::vector<float> rowY;
         for(int y=0; y < nYVoxels; ++y){
	     cv  = sqrt(theMap[x][y][fixedVox].getNPairs() );
	     rms = theMap[x][y][fixedVox].calculateRMS();
	     if(cv != 0.0)
	       rms = (float)(rms / cv );
	     else
	       rms = 0.0;    
	     rowY.push_back( rms );	     	  
	  }
	  return_vec.push_back(rowY);
       }
       return return_vec;
    }
    
    else if(face == 4){         
       for(int y=0; y < nYVoxels; ++y){
        std::vector<float> rowY;
         for(int z=0; z < nZVoxels; ++z){
	     cv  = sqrt(theMap[fixedVox][y][z].getNPairs() );
	     rms = theMap[fixedVox][y][z].calculateRMS();
	     if(cv != 0.0)
	       rms = (float)(rms / cv );
	     else
	       rms = 0.0;    
	     rowY.push_back( rms );	     	  
	  }
	  return_vec.push_back(rowY);
       }
       return return_vec;
    }
    
    else{
       std::cout << "ERROR! Unkown TPC facing: " << face << std::endl;
       return return_vec;    
    }
              	  
}

void PCAResult::doPCA(const std::vector<elecInfo> &points) {
  std::cout << "Do PCA..." << std::endl;
  
  TVector3 outputCentroid;
  std::cout << outputCentroid[0] << " " << outputCentroid[1] << " " << outputCentroid[2] << std::endl;
  std::pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  std::vector<TVector3> outputEigenVecs;
  
  double meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;

  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }

  if (nThreeDHits == 0) {
    
    return; // FAIL FROM NO INPUT POINTS
  }

  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid.SetX(meanPosition[0]);
  outputCentroid.SetY(meanPosition[1]);
  outputCentroid.SetZ(meanPosition[2]);

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
  /*
  std::cout << "Point: " << outputCentroid[0] << " " << outputCentroid[1] << " " << outputCentroid[2] << std::endl;
  centroid = outputCentroid;
  endPoints = outputEndPoints;
  length = outputLength;
  eVals = outputEigenValues;
  std::cout << "Size: " << outputEigenVecs.size() << std::endl;
  eVecs = outputEigenVecs;
  */
  std::cout << startX << std::endl;
  if (outputEndPoints.first(0) > outputEndPoints.second(0)) {
    startX =  endPoint1(0);
    startY =  endPoint1(1);
    startZ =  endPoint1(2);
  }
  
  else {
    std::cout << endPoint2(0) << std::endl;
    startX =  endPoint2(0);
    startY =  endPoint2(1);
    startZ =  endPoint2(2);
  }

  unitX = outputEigenVecs[0](0);
  unitY = outputEigenVecs[0](1);
  unitZ = outputEigenVecs[0](2);
  
  if (unitX > 0.0) {
      unitX *= -1.0;
      unitY *= -1.0;
      unitZ *= -1.0;
  }
  
  std::cout << "Done with PCA..." << std::endl;

}

std::vector<double> PCAResult::getStartPoints(){
  std::vector<double> return_vec;
  return_vec.push_back(startX);
  return_vec.push_back(startY);
  return_vec.push_back(startZ);
  return return_vec;   

}

std::vector<double> PCAResult::getUnits(){
  std::vector<double> return_vec;
  return_vec.push_back(unitX);
  return_vec.push_back(unitY);
  return_vec.push_back(unitZ);
  return return_vec;   

}

void SCECalib::loadTruthMap(bool isFwd)
{
  TFile* fileTruth = new TFile("data/dispOutput_MicroBooNE_E273.root");
  
  std::string treeVar = "bkwd";
  std::string recoVar = "reco";
  
  if(isFwd){
    treeVar = "fwd";
    recoVar = "true";
  }
    
  TTreeReader reader(Form("SpaCEtree_%sDisp", treeVar.c_str() ), fileTruth);
  TTreeReaderValue<Double_t> reco_x(reader, Form("x_%s.data_%sDisp", recoVar.c_str(), treeVar.c_str() ) );
  TTreeReaderValue<Double_t> reco_y(reader, Form("y_%s.data_%sDisp", recoVar.c_str(), treeVar.c_str() ) );
  TTreeReaderValue<Double_t> reco_z(reader, Form("z_%s.data_%sDisp", recoVar.c_str(), treeVar.c_str() ) );
  TTreeReaderValue<Double_t> Dx(reader, Form("Dx.data_%sDisp", treeVar.c_str() ) );
  TTreeReaderValue<Double_t> Dy(reader, Form("Dy.data_%sDisp", treeVar.c_str() ) );
  TTreeReaderValue<Double_t> Dz(reader, Form("Dz.data_%sDisp", treeVar.c_str() ) );
  TTreeReaderValue<Int_t> elecFate(reader, Form("elecFate.data_%sDisp", treeVar.c_str() ) );
       
  double trueDeltaX[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
  double trueDeltaY[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
  double trueDeltaZ[nCalibDivisions_x+1][nCalibDivisions_y+1][nCalibDivisions_z+1];
  

  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
      {
        trueDeltaX[x][y][z] = 0.0;
        trueDeltaY[x][y][z] = 0.0;
        trueDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  while (reader.Next())
  {
    int xIdx = (int)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx));
    int yIdx = (int)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly));
    int zIdx = (int)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz));
    
    //std::cout << xIdx << " " << yIdx << " " << zIdx << " " << *Dx << " " << *Dy << " " << *Dz << std::endl;
    
    if (*elecFate == 1) {
      trueDeltaX[xIdx][yIdx][zIdx] = *Dx;
      trueDeltaY[xIdx][yIdx][zIdx] = *Dy;
      trueDeltaZ[xIdx][yIdx][zIdx] = *Dz;
                  
    }
    
    else {
      trueDeltaX[xIdx][yIdx][zIdx] = -999;
      trueDeltaY[xIdx][yIdx][zIdx] = -999;
      trueDeltaZ[xIdx][yIdx][zIdx] = -999;
    }
  }

  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
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

  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
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

  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
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
  
  
  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
      {
	 //if(isFwd) std::cout << x << " " << y << " " << z << " " << trueDeltaX[x][y][z]<< " " << trueDeltaY[x][y][z] << " " << trueDeltaZ[x][y][z] << std::endl;
	 distortionVoxel *vox;
	
	 if(isFwd) 
	   vox = xFwdTruthMap.getVoxel(x, y, z);
	 else
	   vox = xTruthMap.getVoxel(x, y, z);

	 vox->addDistortion(trueDeltaX[x][y][z]);
	 
	 if(isFwd) 
	   vox = yFwdTruthMap.getVoxel(x, y, z);
	 else
	   vox = yTruthMap.getVoxel(x, y, z);
	 
	 vox->addDistortion(trueDeltaY[x][y][z]);
	 
	 if(isFwd) 
	   vox = zFwdTruthMap.getVoxel(x, y, z);
	 else
	   vox = zTruthMap.getVoxel(x, y, z);
	   
	 vox->addDistortion(trueDeltaZ[x][y][z]);
      }
    }
  }    
  
  return;
}

void SCECalib::loadTruthMap(std::vector<distortionMap>& inputMaps){
  for(int x = 0; x <= nCalibDivisions_x; x++)
  {
    for(int y = 0; y <= nCalibDivisions_y; y++)
    {
      for(int z = 0; z <= nCalibDivisions_z; z++)
      {
	 distortionVoxel *vox;		 
	 
	 vox = inputMaps[0].getVoxel(x, y, z);
	 float cv = vox->calculateMedian();
	 vox = xTruthMap.getVoxel(x, y, z);
	 vox->clearDistortions();
	 vox->clearWeights();
	 vox->addDistortion(cv);
	 
	 vox = inputMaps[1].getVoxel(x, y, z);
	 cv = vox->calculateMedian();
	 vox = yTruthMap.getVoxel(x, y, z);
	 vox->clearDistortions();
	 vox->clearWeights();
	 vox->addDistortion(cv);
	 
	 vox = inputMaps[2].getVoxel(x, y, z);
	 cv = vox->calculateMedian();
	 vox = zTruthMap.getVoxel(x, y, z);
	 vox->clearDistortions();
	 vox->clearWeights();
	 vox->addDistortion(cv);



      }
    }
  }

}

std::vector<trackInfo> SCECalib::getLaserTrackSet(){
  std::vector<trackInfo> tracks;

  TFile* mapFile;

  mapFile = new TFile(inputFileLaser.c_str(),"READ");

  TTreeReader readerLasers("lasers", mapFile);
  TTreeReader readerTracks("tracks", mapFile);

  TTreeReaderValue<TVector3> track_start(readerLasers, "entry");
  TTreeReaderValue<TVector3> track_end(readerLasers, "exit");
  TTreeReaderArray<TVector3> track_points(readerTracks, "track");

  int nTracks = 0;
  while (readerLasers.Next())
  {
    readerTracks.Next();
  
    nTracks++;
  
    std::vector<double> s0, e0;
  
    if ((*track_start).y() > (*track_end).y()) {
      s0.push_back( doCoordTransformX((*track_start).x()) );
      s0.push_back( doCoordTransformY((*track_start).y()) );
      s0.push_back( doCoordTransformZ((*track_start).z()) );
  
      e0.push_back( doCoordTransformX((*track_end).x()) );
      e0.push_back( doCoordTransformY((*track_end).y()) );
      e0.push_back( doCoordTransformZ((*track_end).z()) );
    }
    
    else {
      s0.push_back( doCoordTransformX((*track_end).x()) );
      s0.push_back( doCoordTransformY((*track_end).y()) );
      s0.push_back( doCoordTransformZ((*track_end).z()) );
  
      e0.push_back( doCoordTransformX((*track_start).x()) );
      e0.push_back( doCoordTransformY((*track_start).y()) );
      e0.push_back( doCoordTransformZ((*track_start).z()) );
    }
    
    double trackLength = sqrt(pow(s0[0]-e0[0],2.0)+pow(s0[1]-e0[1],2.0)+pow(s0[2]-e0[2],2.0));
    
    double theta = acos((e0[1]-s0[1])/trackLength);
    double phi = acos((e0[2]-s0[2])/(trackLength*sin(theta)));
    if (e0[0] > s0[0]) {
      phi = -1.0*fabs(phi);
    }
    else {
      phi = fabs(phi);
    }
  
    trackInfo track;
    elecInfo electron;
    track.updateEndpoints(s0, e0);
    track.updateAngles(theta, phi);    
  
    for(int j = 0; j < track_points.GetSize(); j++)
    {
      if (((j % 10) != 0) && (j != track_points.GetSize()-1)) continue;
      std::vector<double> xFormedPoints;
      xFormedPoints.push_back(std::min(std::max(0.0,doCoordTransformX(track_points[j].x())),Lx) );
      xFormedPoints.push_back(std::min(std::max(0.0,doCoordTransformY(track_points[j].y())),Ly) );
      xFormedPoints.push_back(std::min(std::max(0.0,doCoordTransformZ(track_points[j].z())),Lz) );
      
      electron.updateSMod( xFormedPoints  );
      track.electrons.push_back(electron);
    }
    
    tracks.push_back(track);
  }

  return tracks;
  }//end of getLArSoftTrackSet()
  
  
std::vector<trackInfo> SCECalib::getCosmicTrackSet(bool isCalibrated = false){
  
    std::vector<trackInfo> tracks;
    
    //Check to see if the 2D face maps are empty
    //If they are, and isCalibrated is set to 0, throw an error
    
    if(isCalibrated && (xTopMap.size() == 0 || yTopMap.size() == 0 || zTopMap.size() == 0 || xBottomMap.size() == 0 || yBottomMap.size() == 0 || zBottomMap.size() == 0 || xUpstreamMap.size() == 0 || yUpstreamMap.size() == 0 || zUpstreamMap.size() == 0 || xDownstreamMap.size() == 0 || yDownstreamMap.size() == 0 || zDownstreamMap.size() == 0 || xCathodeMap.size() == 0 || yCathodeMap.size() == 0 || zCathodeMap.size() == 0) ){
    
       std::cout << "ERROR! Cannot get calibrated offsets with empty calibration maps" << std::endl;
       return tracks;
    }
    
    std::vector<double> distFromFaces;
    distFromFaces.push_back(maxXdist);
    distFromFaces.push_back(maxYdist);
    distFromFaces.push_back(maxZdist);
    
    std::vector<double> tpcEdges;
    tpcEdges.push_back(Lx);
    tpcEdges.push_back(Ly);
    tpcEdges.push_back(Lz);
    
    TFile* mapFile;
    
    mapFile = new TFile(inputFileCosmic.c_str(),"READ");

    TTreeReader reader("SCEtree", mapFile);
    TTreeReaderValue<Double_t> track_startX(reader, "track_startX");
    TTreeReaderValue<Double_t> track_startY(reader, "track_startY");
    TTreeReaderValue<Double_t> track_startZ(reader, "track_startZ");
    TTreeReaderValue<Double_t> track_endX(reader, "track_endX");
    TTreeReaderValue<Double_t> track_endY(reader, "track_endY");
    TTreeReaderValue<Double_t> track_endZ(reader, "track_endZ");
    TTreeReaderValue<Int_t> nPoints(reader, "track_nPoints");
    TTreeReaderArray<Double_t> pointX(reader, "track_pointX");
    TTreeReaderArray<Double_t> pointY(reader, "track_pointY");
    TTreeReaderArray<Double_t> pointZ(reader, "track_pointZ");
    TTreeReaderValue<Double_t> track_MCS(reader, "track_MCS_measurement");
    TTreeReaderValue<Double_t> track_t0(reader, "track_t0");
    
    trackInfo track;
    elecInfo electron;
    
    int containedTracks = 0;
    
    TF1 *weightFunc;
    if (isMC == true) {
      weightFunc = new TF1("weightFunc","(1.0/75.0)*(x-5.0)^2 + (2.0/3.0)",0.0,10.0);
    }
    else {
      weightFunc = new TF1("weightFunc","(2.0/75.0)*(x-5.0)^2 + (1.0/3.0)",0.0,10.0);
    }
    
    TRandom3 *rand = new TRandom3(randSeed);
      
    int inputTrackNum = -1;
    int nTracks = 0;
    while (reader.Next())
    {
      inputTrackNum++;
      if ((inputTrackNum < minInputTrackNum) || (inputTrackNum > maxInputTrackNum)) continue;
      
      if ((maxCosmicTracks != -1) && (nTracks >= maxCosmicTracks)) continue;
      
      if (*nPoints < 3) continue;
      
      std::vector<double> s0, e0;
      std::vector<double> s1, e1;      
    
      //Double_t x_offset = 0.002564*(*track_t0); // TEMP OFFSET, CHRIS WILL FIX BUG SOON
      double x_offset = 0.0;
      
      if (*track_startY > *track_endY) {        
	
	s0.push_back( doCoordTransformX(*track_startX + x_offset) );
	s0.push_back( doCoordTransformY(*track_startY) );
	s0.push_back( doCoordTransformZ(*track_startZ) );
	
	e0.push_back( doCoordTransformX(*track_endX + x_offset) );
	e0.push_back( doCoordTransformY(*track_endY) );
	e0.push_back( doCoordTransformZ(*track_endZ) );
	
	
      }
      
      else {
        s0.push_back( doCoordTransformX(*track_endX + x_offset) );
	s0.push_back( doCoordTransformY(*track_endY) );
	s0.push_back( doCoordTransformZ(*track_endZ) );
    
        e0.push_back( doCoordTransformX(*track_startX + x_offset) );
	e0.push_back( doCoordTransformY(*track_startY) );
	e0.push_back( doCoordTransformZ(*track_startZ) );
	
	//std::cout << s0[0] << " " << s0[1] << " " << s0[2] << " " << e0[0] << " " << e0[1] << " " << e0[2] << std::endl;
      }      

      //Zhist_1s->Fill(s0[2]);
      //Zhist_1e->Fill(e0[2]);
      
      track.updateEndpoints(s0, e0); // update the endpoints of the track. This enables us to use track methods 
                        
      
      
      
      
      if(track.isNotContained(distFromFaces, tpcEdges) ) continue;
      
      //Zhist_2s->Fill(s0[2]);
      //Zhist_2e->Fill(e0[2]);
      
      
      
      if (track.isSomething(distFromFaces, tpcEdges) ) continue;
     
                   
      //Zhist_3s->Fill(s0[2]);
      //Zhist_3e->Fill(e0[2]);
    
      if (track.isCathodeExiting(distFromFaces, tpcEdges) && (*track_MCS < 1000.0*minTrackMCS_anode)) continue;
      
       

      //Zhist_4s->Fill(s0[2]);
      //Zhist_4e->Fill(e0[2]);
      
      if (track.isAnodeExiting(distFromFaces, tpcEdges) && (*track_MCS < 1000.0*minTrackMCS_cathode)) continue;

      //Zhist_5s->Fill(s0[2]);
      //Zhist_5e->Fill(e0[2]);

      if (track.isCrossing(distFromFaces, tpcEdges) && (*track_MCS < 1000.0*minTrackMCS_crossing)) continue;
      
      //Zhist_6s->Fill(s0[2]);
      //Zhist_6e->Fill(e0[2]);

      double randNum = rand->Uniform(1.0);
      if (randNum > std::max(weightFunc->Eval(s0[2]),weightFunc->Eval(e0[2]))) continue;
      
      
      nTracks++;
      ++containedTracks;
      

      //Zhist_7s->Fill(s0[2]);
      //Zhist_7e->Fill(e0[2]);
      
      // Correct track end point for cathode-piercing track (end furthest from cathode)
      //if (((s0[0] < (Lx - maxXdist)) && (e0[0] < maxXdist)) || ((e0[0] < (Lx - maxXdist)) && (s0[0] < maxXdist))) {
      if(track.isAnodeExiting(distFromFaces, tpcEdges) ){
        
	
	
	if (s0[0] < e0[0]) {
          std::vector<double> endPoints;
	  endPoints.push_back(0.0);
	  endPoints.push_back(s0[1]+ ( (isCalibrated ? getCalibOffset(s0,axisType::yAxis, 4) : getTruthOffset(s0,axisType::yAxis) ) ) );
	  endPoints.push_back(s0[2]+ ( (isCalibrated ? getCalibOffset(s0,axisType::zAxis, 4) : getTruthOffset(s0,axisType::zAxis) ) ) );
	  e0[0] += SCEFactor*getTruthOffset(endPoints, axisType::xAxis, true);
	  //std::cout << getTruthOffset(endPoints, axisType::xAxis, true) << std::endl;
	}
	
	else {
          std::vector<double> endPoints;
	  endPoints.push_back(0.0);
	  endPoints.push_back(e0[1]+ ( (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 4) : getTruthOffset(e0,axisType::yAxis) ) ) );
	  endPoints.push_back(e0[2]+ ( (isCalibrated ? getCalibOffset(e0,axisType::zAxis, 4) : getTruthOffset(e0,axisType::zAxis) ) ) );	  
	  s0[0] += SCEFactor*getTruthOffset(endPoints, axisType::xAxis, true);
	  
	}
      }

      
      //Fixes track endpoints, based on 1D face calib.
      //Change this s.t. you use face calibrations, at each x,y,z point once they are calibrated.  
      if(s0[0] < maxXdist) {
        //std::cout << getTruthOffset(s0,axisType::yAxis) << std::endl;
	s1.push_back( 0.0 );
        s1.push_back( s0[1] + ( isCalibrated ? getCalibOffset(s0,axisType::yAxis, 4) : getTruthOffset(s0,axisType::yAxis) ) );
        s1.push_back( s0[2] + ( isCalibrated ? getCalibOffset(s0,axisType::zAxis, 4) : getTruthOffset(s0,axisType::zAxis) ) );
      }
      else if (s0[0] > (Lx - maxXdist)) {
        s1.push_back( Lx );
        s1.push_back( s0[1] + (isCalibrated ? getCalibOffset(s0,axisType::yAxis, 4) : getTruthOffset(s0,axisType::yAxis) )  );
        s1.push_back( s0[2] + (isCalibrated ? getCalibOffset(s0,axisType::zAxis, 4) : getTruthOffset(s0,axisType::zAxis) )  );
      }
      
      else if (std::min( fabs(s0[1]),fabs(Ly-s0[1]) )  < std::min(fabs(s0[2]),fabs(Lz-s0[2])) ) {
        if (fabs(s0[1]) < fabs(Ly-s0[1])) {
          s1.push_back( s0[0] + (isCalibrated ? getCalibOffset(s0,axisType::xAxis, 1) : getTruthOffset(s0,axisType::xAxis) ) );
    	  s1.push_back( 0.0 );
    	  s1.push_back( s0[2] + (isCalibrated ? getCalibOffset(s0,axisType::zAxis, 1) : getTruthOffset(s0,axisType::zAxis) ) );
        }
        else {
          
	  
	  s1.push_back( s0[0] + (isCalibrated ? getCalibOffset(s0,axisType::xAxis, 0) : getTruthOffset(s0,axisType::xAxis) ) );
    	  s1.push_back( Ly );
    	  s1.push_back( s0[2] + (isCalibrated ? getCalibOffset(s0,axisType::zAxis, 0) : getTruthOffset(s0,axisType::zAxis) ) );
        }
      }
      
      else {
        if (fabs(s0[2]) < fabs(Lz-s0[2])) {
          s1.push_back( s0[0] + (isCalibrated ? getCalibOffset(s0,axisType::xAxis, 2) : getTruthOffset(s0,axisType::xAxis) ) );
    	  s1.push_back( s0[1] + (isCalibrated ? getCalibOffset(s0,axisType::yAxis, 2) : getTruthOffset(s0,axisType::yAxis) ) );
    	  s1.push_back( 0.0 );
        }
        else {
          s1.push_back( s0[0] + (isCalibrated ? getCalibOffset(s0,axisType::xAxis, 3) : getTruthOffset(s0,axisType::xAxis) ) );
    	  s1.push_back( s0[1] + (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 3) : getTruthOffset(s0,axisType::yAxis) ) );
    	  s1.push_back( Lz );
        }
      }
    
      if(e0[0] < maxXdist) {
        e1.push_back( 0.0 );
        e1.push_back( e0[1] + (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 4) : getTruthOffset(e0,axisType::yAxis) ) );
        e1.push_back( e0[2] + (isCalibrated ? getCalibOffset(e0,axisType::zAxis, 4) : getTruthOffset(e0,axisType::zAxis) ) );
      }
      
      else if (e0[0] > (Lx - maxXdist)) {
        e1.push_back( Lx );
        e1.push_back( e0[1] + (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 4) : getTruthOffset(e0,axisType::yAxis) ) );
        e1.push_back( e0[2] + (isCalibrated ? getCalibOffset(e0,axisType::zAxis, 4) : getTruthOffset(e0,axisType::zAxis) ) );
      }
      
      else if (std::min(fabs(e0[1]),fabs(Ly-e0[1])) < std::min(fabs(e0[2]),fabs(Lz-e0[2]))) {
        if (fabs(e0[1]) < fabs(Ly-e0[1])) {
          e1.push_back( e0[0] + (isCalibrated ? getCalibOffset(e0,axisType::xAxis, 1) : getTruthOffset(e0,axisType::xAxis) ) );
    	  e1.push_back( 0.0 );
    	  e1.push_back( e0[2] + (isCalibrated ? getCalibOffset(e0,axisType::zAxis, 1) : getTruthOffset(e0,axisType::zAxis) ) );
        }
        else {
          //std::cout << getTruthOffset(s0,axisType::xAxis) << std::endl;
	  e1.push_back( e0[0] + (isCalibrated ? getCalibOffset(e0,axisType::xAxis, 0) : getTruthOffset(e0,axisType::xAxis) ) );
    	  e1.push_back( Ly );
    	  e1.push_back( e0[2] + (isCalibrated ? getCalibOffset(e0,axisType::zAxis, 0) : getTruthOffset(e0,axisType::zAxis) ) );
        }
      }
      
      else {
        if (fabs(e0[2]) < fabs(Lz-e0[2])) {
          e1.push_back( e0[0] + (isCalibrated ? getCalibOffset(e0,axisType::xAxis, 2) : getTruthOffset(e0,axisType::xAxis) ) );
    	  e1.push_back( e0[1] + (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 2) : getTruthOffset(e0,axisType::yAxis) ) );
    	  e1.push_back( 0.0 );
        }
        else {
          e1.push_back( e0[0] + (isCalibrated ? getCalibOffset(e0,axisType::xAxis, 3) : getTruthOffset(e0,axisType::xAxis) ) );
    	  e1.push_back( e0[1] + (isCalibrated ? getCalibOffset(e0,axisType::yAxis, 3) : getTruthOffset(e0,axisType::yAxis) ) );
    	  e1.push_back( Lz );
        }
      }
      
      double trackLength = sqrt(pow(s1[0]-e1[0],2.0)+pow(s1[1]-e1[1],2.0)+pow(s1[2]-e1[2],2.0));
      //std::cout << s1[0] << " " << s1[1] << " " << s1[2] << " : " << e1[0] << " " << e1[1] << " " << e1[2] << std::endl;
      //std::cout << "Track L: " << trackLength << std::endl;
      track.electrons.clear();
      
      for(int j = 0; j < *nPoints; j++)
      {
        
	if (((j % 10) != 0) && (j != *nPoints-1)) continue;
        std::vector<double> xFormedPoints;
	if (track.isAnodeExiting(distFromFaces, tpcEdges)) {
	  if (s1[0] < e1[0]) {
	    xFormedPoints.push_back( doCoordTransformX(pointX[j]+x_offset) + (isCalibrated ? getCalibOffset(s1, axisType::xAxis, 4) : getTruthOffset(s1, axisType::xAxis, true) ) );
	  }
	  
	  else {
	    xFormedPoints.push_back( doCoordTransformX(pointX[j]+x_offset) + (isCalibrated ? getCalibOffset(e1, axisType::xAxis, 4) : getTruthOffset(e1 ,axisType::xAxis, true) ) );
	  }
        
	}
        
	else {
          xFormedPoints.push_back( doCoordTransformX(pointX[j]+x_offset) );
        }
        
       
       xFormedPoints.push_back(doCoordTransformY(pointY[j]) );
       xFormedPoints.push_back(doCoordTransformZ(pointZ[j]) );       
        
       electron.updateSMod( xFormedPoints );
       track.electrons.push_back(electron);
	
      }
      
      double theta = acos( (e1[1]-s1[1])/trackLength );
      double phi = acos( (e1[2]-s1[2])/(trackLength*sin(theta)) );
      if (e1[0] > s1[0]) {
        phi = -1.0*fabs(phi);
      }
      else {
        phi = fabs(phi);
      }
      
      track.energy = *track_MCS/1000.0;
      track.pdgID = 13;
      track.updateEndpoints(s1, e1);//update the start and endpoints after "gluing" to the faces.
      track.updateAngles(theta, phi);
      
      
      tracks.push_back(track);
    }

    outputFile->cd();
    /*
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
*/
    std::cout << "N Tracks: " << nTracks << std::endl;
    std::cout << "Contained tracks: " << containedTracks << std::endl;
    std::cout << "Input Track Num: " << inputTrackNum << std::endl;
    return tracks;
}  

double SCECalib::doCoordTransformX(double inputX) const
{
  double outputX;

  //outputX = Lx - (Lx/2.56)*inputX/100.0;
  if(isMC) {
    outputX = Lx - (Lx/2.5524)*inputX/100.0;
  }
  else {
    outputX = Lx - (Lx/2.58)*inputX/100.0;
  }

  return outputX;
}

double SCECalib::doCoordTransformY(double inputY) const 
{
  double outputY;

  if(isMC) {
    outputY = (Ly/(1.173+1.154))*(inputY+115.4)/100.0;
  }
  else {
    outputY = (Ly/(1.170+1.151))*(inputY+115.1)/100.0;
  }

  return outputY;
}

double SCECalib::doCoordTransformZ(double inputZ) const
{
  double outputZ;

  if(isMC) {
    outputZ = (Lz/(10.368-0.003))*(inputZ-0.3)/100.0;
  }
  else {
    outputZ = (Lz/(10.365+0.007))*(inputZ+0.7)/100.0;
  }
  
  return outputZ;
}

double SCECalib::getTruthOffset(std::vector<double> sVec, axisType comp, bool isFwd)
{
  double offset = 0.0;
  
  if (sVec[0] < 0.0) {
    sVec[0] = 0.0;
  }
  if (sVec[0] > Lx) {
    sVec[0] = Lx;
  }

  if (sVec[1] < 0.0) {
    sVec[1] = 0.0;
  }
  if (sVec[1] > Ly) {
    sVec[1] = Ly;
  }

  if (sVec[2] < 0.0) {
    sVec[2] = 0.0;
  }
  if (sVec[2] > Lz) {
    sVec[2] = Lz;
  }
  
  int x = (int)TMath::Nint(nCalibDivisions_x*(sVec[0]/Lx));
  int y = (int)TMath::Nint(nCalibDivisions_y*(sVec[1]/Ly));
  int z = (int)TMath::Nint(nCalibDivisions_z*(sVec[2]/Lz));
  
  distortionVoxel *vox;
  if(isFwd){
    if(comp == 0){
      vox = xFwdTruthMap.getVoxel(x, y, z);
    }
  
    else if(comp == 1){
      vox = yFwdTruthMap.getVoxel(x, y, z);
    }
  
    else if(comp == 2){
      vox = zFwdTruthMap.getVoxel(x, y, z);
    }
  
    else
      return 0.0;
  }
  
  else{
    if(comp == 0){
      vox = xTruthMap.getVoxel(x, y, z);
    }
  
    else if(comp == 1){
      vox = yTruthMap.getVoxel(x, y, z);
    }
  
    else if(comp == 2){
      vox = zTruthMap.getVoxel(x, y, z);
    }
  
    else
      return 0.0;
  }
       
  vox->getDistortions();
  
  offset = vox->getDistortions()[0];
     
  if (!isFwd && (comp == 0) && (offset < 0.0)) {
    offset = 0.0;
  }
  
  if (isFwd && (comp == 1) && (offset > 0.0)) {
    offset = 0.0;
  }
  
  return offset;
}

double SCECalib::getCalibOffset(std::vector<double> sVec, axisType comp, int face)
{
  double offset = 0.0;
  
  if (sVec[0] < 0.0) {
    sVec[0] = 0.0;
  }
  if (sVec[0] > Lx) {
    sVec[0] = Lx;
  }

  if (sVec[1] < 0.0) {
    sVec[1] = 0.0;
  }
  if (sVec[1] > Ly) {
    sVec[1] = Ly;
  }

  if (sVec[2] < 0.0) {
    sVec[2] = 0.0;
  }
  if (sVec[2] > Lz) {
    sVec[2] = Lz;
  }
  
  int x = (int)TMath::Nint(nCalibDivisions_x*(sVec[0]/Lx));
  int y = (int)TMath::Nint(nCalibDivisions_y*(sVec[1]/Ly));
  int z = (int)TMath::Nint(nCalibDivisions_z*(sVec[2]/Lz));
  if(face == 0){
    if(comp == 0){
      offset = xTopMap[x][z];
    }
  
    else if(comp == 1){
      offset = yTopMap[x][z];
    }
  
    else if(comp == 2){
      offset = zTopMap[x][z];
    }
  
    else
      return 0.0;
  
  }
  
  else if(face == 1){
    if(comp == 0){
      offset = xBottomMap[x][z];
    }
  
    else if(comp == 1){
      offset = yBottomMap[x][z];
    }
  
    else if(comp == 2){
      offset = zBottomMap[x][z];
    }
  
    else
      return 0.0;
  }

  else if(face == 2){
    if(comp == 0){
      offset = xUpstreamMap[x][y];
    }
  
    else if(comp == 1){
      offset = yUpstreamMap[x][y];
    }
  
    else if(comp == 2){
      offset = zUpstreamMap[x][y];
    }
  
    else
      return 0.0;
  }
  
  else if(face == 3){
    if(comp == 0){
      offset = xDownstreamMap[y][y];
    }
  
    else if(comp == 1){
      offset = yDownstreamMap[y][y];
    }
  
    else if(comp == 2){
      offset = zDownstreamMap[y][y];
    }
  
    else
      return 0.0;
  }
  
  else if(face == 4){
    if(comp == 0){
      offset = xCathodeMap[y][z];
    }
  
    else if(comp == 1){
      offset = yCathodeMap[y][z];
    }
  
    else if(comp == 2){
      offset = zCathodeMap[y][z];
    }
  
    else
      return 0.0;
  }
  
  else
    return 0.0;    
  
     
  if ((comp == 0) && (offset < 0.0)) {
    offset = 0.0;
  }
  
  return offset;
}

void SCECalib::Intialize(){
  inputFileLaser = "data/laserDataSCE_NEW.root";

  inputFileCosmic = "/uboone/data/users/joelam/MCExtForSCE.root";

  outputFile = new TFile("output.root","RECREATE");

  Lx = 2.5;
  Ly = 2.5;
  Lz = 10.0;
  
  cosmicTruthMode = 1;

  isMC = true;
  doBulk = false;


  relAngleCut = 20.0;
  maxXdist = 0.05;
  maxYdist = 0.20;
  maxZdist = 0.20;
  maxDistFactor = 3.0;
  distScale     = 0.01;

  minInputTrackNum = 0;
  maxInputTrackNum = 1150000;
  
  SCEFactor = 1.0;

  maxCosmicTracks = 100;

  minTrackMCS_anode = 3.3;
  minTrackMCS_cathode = 1.7;
  minTrackMCS_crossing = 1.4;

  nCalibDivisions = 25; // MICROBOONE
//Int_t nCalibDivisions = 18; // PROTODUNE-SP

  piVal = TMath::Pi();
  
  nCalibDivisions_x = nCalibDivisions;
  nCalibDivisions_y = TMath::Nint((Ly/Lx)*((double)nCalibDivisions));
  nCalibDivisions_z = TMath::Nint((Lz/Lx)*((double)nCalibDivisions));

  randSeed = 1; //set this to 0 to get a random seed

}
/*
int SCECalib::doCalibration(const std::vector<trackInfo> &laserTracks, const std::vector<trackInfo> &cosmicTracks){
  
  
  std::vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);

  
  std::cout << track.x0 << " " << track.y0 << " " << track.z0 << std::endl;
  std::cout << calibTrack.x0 << " " << calibTrack.y0 << " " << calibTrack.z0 << std::endl;
  
  std::vector<elecInfo> electronsNonCalib = track.electrons;
  std::vector<elecInfo> electronsCalib    = calibTrack.electrons;
  
  std::cout << electronsNonCalib.size() << ", " << electronsCalib.size() << std::endl;
  
  for(unsigned int i = 0; i != electronsCalib.size(); ++i){
     std::cout << electronsNonCalib.at(i).x << ", " << electronsNonCalib.at(i).y << ", " << electronsNonCalib.at(i).z << std::endl;
     std::cout << electronsCalib.at(i).x << ", " << electronsCalib.at(i).y << ", " << electronsCalib.at(i).z << std::endl;
  }
  
  double x_true, y_true, z_true;
  double x_reco, y_reco, z_reco;
  double Dx, Dy, Dz;
  int elecFate;
  
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
  
  std::cout << "Number of Tracks: " << cosmicCalibTracks.size() << std::endl;
  
  switch (stepsToRun)  {
   case(bulkOnly) : {
                                      
    std::vector<distortionMap> calibDistortions = doCosmicCosmicCalib(cosmicCalibTracks);
    
  
    std::vector< std::vector < std::vector<float> > > xMap     =  calibDistortions[0].calculateMap();
    std::vector< std::vector < std::vector<float> > > xMapErrs =  calibDistortions[0].calculateMapErrors();
  
    std::vector< std::vector < std::vector<float> > > yMap     =  calibDistortions[1].calculateMap();
    std::vector< std::vector < std::vector<float> > > yMapErrs =  calibDistortions[1].calculateMapErrors();
  
    std::vector< std::vector < std::vector<float> > > zMap     =  calibDistortions[2].calculateMap();
    std::vector< std::vector < std::vector<float> > > zMapErrs =  calibDistortions[2].calculateMapErrors();
  
    return 0;
   }
   
   case(faceOnly) : {
     std::vector<distortionMap> faceCalibDistortions = doCalibFaces(cosmicCalibTracks, 50, 15);
     calculate2DMaps(faceCalibDistortions);
     return 0;
   
   }
   
   case(bulkAndFace) : {
     std::vector<distortionMap> calibDistortions = doCosmicCosmicCalib(cosmicCalibTracks);
     loadTruthMap(calibDistortions);     
     std::vector<distortionMap> faceCalibDistortions = doCalibFaces(cosmicCalibTracks, 50, 15);
     calculate2DMaps(faceCalibDistortions);
     faceCalibrated = true;
     return 0;
     
     
   
   }
   
   case(fullCalib) : {
     std::vector<distortionMap> calibDistortions = doCosmicCosmicCalib(cosmicCalibTracks);
     loadTruthMap(calibDistortions);     
     std::vector<distortionMap> faceCalibDistortions = doCalibFaces(cosmicCalibTracks, 50, 15);
     faceCalibrated = true;
     calculate2DMaps(faceCalibDistortions);              
     std::vector<distortionMap> finalDistortions = doCosmicCosmicCalib(cosmicCalibTracks);
     return 0;
     
   
   }
   
   default : return 1;
   
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
        
        if(calibWeight[x][y][z] > 0.0)
          elecFate = 1;
	else
          elecFate = 0;
    
        T_calib->Fill();
      }
    }
  }
  
  //probably shouldn't get here...
  return -1;
}*/

std::vector<distortionMap> SCECalib::doCosmicCosmicCalib(const std::vector<calibTrackInfo> &cosmicCalibTracks)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;
  
  distortionMap deltaXMap;
  distortionMap deltaYMap;
  distortionMap deltaZMap;

  unsigned int numCosmicTracks = cosmicCalibTracks.size();

  std::vector<double> POAparams;
  std::vector<double> POAparamsDistorted;
  double distVal;
  double distValDistorted;
  double distWeight;
  double xVal;
  double yVal;
  double zVal;
  double xValDistorted;
  double yValDistorted;
  double zValDistorted;
  int xCalibLowIndex;
  int xCalibHighIndex;
  int yCalibLowIndex;
  int yCalibHighIndex;
  int zCalibLowIndex;
  int zCalibHighIndex;
  double xCalibFrac;
  double yCalibFrac;
  double zCalibFrac;
  double tempFactor;

  double crossDistX;
  double crossDistY;
  double crossDistZ;
  double crossDistX_mod;
  double crossDistY_mod;
  double crossDistZ_mod;
  
  int trackNum1;
  int trackNum2;
  int crossType = 3;
  
  double rawDeltaX = 0.0;
  double rawDeltaY = 0.0;
  double rawDeltaZ = 0.0;
  
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
    

  for(int i = 0; i < numCosmicTracks; i++)
  {
    std::cout << "COSMIC-COSMIC " << i << std::endl;

    calibTrackA = cosmicCalibTracks.at(i);
    if(calibTrackA.electrons.size() < 3) continue;
    if((cosmicTruthMode == 0) && (calibTrackA.calibFlag == false)) continue;

    for(int j = i+1; j < numCosmicTracks; j++)
    {
      calibTrackB = cosmicCalibTracks.at(j);
      if(calibTrackB.electrons.size() < 3) continue;
      if((cosmicTruthMode == 0) && (calibTrackB.calibFlag == false)) continue;      
      
      double dTheta = fabs(acos(((calibTrackA.x1-calibTrackA.x0)*(calibTrackB.x1-calibTrackB.x0) + (calibTrackA.y1-calibTrackA.y0)*(calibTrackB.y1-calibTrackB.y0) + (calibTrackA.z1-calibTrackA.z0)*(calibTrackB.z1-calibTrackB.z0))/(sqrt(pow((calibTrackA.x1-calibTrackA.x0),2.0)+pow((calibTrackA.y1-calibTrackA.y0),2.0)+pow((calibTrackA.z1-calibTrackA.z0),2.0))*sqrt(pow((calibTrackB.x1-calibTrackB.x0),2.0)+pow((calibTrackB.y1-calibTrackB.y0),2.0)+pow((calibTrackB.z1-calibTrackB.z0),2.0)))));
      
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;
      
      //std::cout << i << " " << calibTrackA.x0 << " " << calibTrackA.y0 << " " << calibTrackA.z0 << " " << calibTrackA.theta << " " << calibTrackA.phi << std::endl;
      
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

      xCalibLowIndex  = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex  = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex  = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);
      
      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((double) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((double) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((double) zCalibLowIndex);
      

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
      
      //std::cout << xValDistorted << " " << Lx << " " nCalibDivisions_x << " " << xCalibLowIndex << std::endl;
      std::cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << std::endl;

      distortionVoxel *vox; 
      
      if(tempFactor < 0.0)
        std::cout << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << std::endl;
      
      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);            
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibLowIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
     
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);


      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibHighIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);      

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibLowIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
      

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibHighIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
      

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;

      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibLowIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
	


      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibHighIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
      	

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibLowIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
      	

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibHighIndex);              
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal - xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal - yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal - zValDistorted);
    }
  }
  std::vector<distortionMap> returnVec;
  returnVec.push_back(deltaXMap);
  returnVec.push_back(deltaYMap);
  returnVec.push_back(deltaZMap);
  
  return returnVec;
}


calibTrackInfo::calibTrackInfo(trackInfo track){
  pdgID = track.pdgID;
  energy = track.energy;
  x0 = track.x0;
  y0 = track.y0;
  z0 = track.z0;
  x1 = track.x1;
  y1 = track.y1;
  z1 = track.z1;
  theta = track.theta;
  phi = track.phi;
  electrons = track.electrons;
  /*
  std::cout << track.x0 << " " << x0 << std::endl;
  std::cout << track.x1 << " " << x1 << std::endl;
  */
  for(unsigned int j = 0; j < electrons.size(); j++){
      DxVec.push_back(0.0);
      DyVec.push_back(0.0);
      DzVec.push_back(0.0);
  }

    calibFlag = false;    
}

std::vector<calibTrackInfo> SCECalib::makeCalibTracks(const std::vector<trackInfo> &tracks)
{
  std::vector<calibTrackInfo> calibTracks;
  calibTrackInfo calibTrack;

  for(unsigned int i = 0; i < tracks.size(); i++)
  {
   trackInfo track = tracks.at(i);

    calibTrack = calibTrackInfo(track);
        
    calibTracks.push_back(calibTrack);
  }

  return calibTracks;
}

std::vector<double> SCECalib::findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB) const
{ 
  //   (xA,yA,zA)+t(xA_step,yA_step,zA_step)
  
  /*
  double xA = calibTrackA.x0_calib, yA = calibTrackA.y0_calib, zA = calibTrackA.z0_calib;
  double xB = calibTrackB.x0_calib, yB = calibTrackB.y0_calib, zB = calibTrackB.z0_calib;

  double xA_step = -1.*sin(calibTrackA.theta_calib)*sin(calibTrackA.phi_calib);
  double yA_step = cos(calibTrackA.theta_calib);
  double zA_step = sin(calibTrackA.theta_calib)*cos(calibTrackA.phi_calib);

  double xB_step = -1.*sin(calibTrackB.theta_calib)*sin(calibTrackB.phi_calib);
  double yB_step = cos(calibTrackB.theta_calib);
  double zB_step = sin(calibTrackB.theta_calib)*cos(calibTrackB.phi_calib);
*/
  double xA = calibTrackA.x0, yA = calibTrackA.y0, zA = calibTrackA.z0;
  
  double xB = calibTrackB.x0, yB = calibTrackB.y0, zB = calibTrackB.z0;
  
  
  
  double xA_step = -1.*sin(calibTrackA.theta)*sin(calibTrackA.phi);
  double yA_step = cos(calibTrackA.theta);
  double zA_step = sin(calibTrackA.theta)*cos(calibTrackA.phi);
  
  //std::cout << xA_step << " " << yA_step << " " << zA_step << " " << calibTrackA.theta << " " << calibTrackA.phi <<  std::endl;

  double xB_step = -1.*sin(calibTrackB.theta)*sin(calibTrackB.phi);
  double yB_step = cos(calibTrackB.theta);
  double zB_step = sin(calibTrackB.theta)*cos(calibTrackB.phi);
  //perpendicular line between the two tracks
  double x_prep = (yA_step*zB_step)-(yB_step*zA_step);
  double y_prep = (xB_step*zA_step)-(xA_step*zB_step);
  double z_prep = (xA_step*yB_step)-(xB_step*yA_step);
 
  // if cross product is zero then the lines are parallel so return distance = -2
  if (x_prep == 0 && y_prep == 0 && z_prep == 0) {
    std::vector<double> return_vector;
    return_vector.push_back(-2.);
    return_vector.push_back((xA+xB)/2.);
    return_vector.push_back((yA+yB)/2.);
    return_vector.push_back((zA+zB)/2.);   
    return return_vector;
  }
  //normalize the perpendicular line
  double mag_prep = sqrt(pow(x_prep,2)+pow(y_prep,2)+pow(z_prep,2));
  double x_prep_norm = x_prep / mag_prep;
  double y_prep_norm = y_prep / mag_prep;
  double z_prep_norm = z_prep / mag_prep;
 
  //defined to make the math simplier
  double a = y_prep_norm*(xA-xB);
  double b = x_prep_norm*(yA-yB);
  double c = z_prep_norm*(xA-xB);
  double d = x_prep_norm*(zA-zB);
 
  double g = y_prep_norm*xA_step;
  double h = y_prep_norm*xB_step;
  double i = x_prep_norm*yA_step;
  double j = x_prep_norm*yB_step;
  double k = z_prep_norm*xA_step;
  double l = z_prep_norm*xB_step;
  double m = x_prep_norm*zA_step;
  double n = x_prep_norm*zB_step;
 
  double chi = (l-n)/(h-j);
 
  //alpha: "t" for the first line //beta:  "t" for the second line
  double alpha = (chi*(a-b)-c+d)/(k-m-(chi*(g-i)));
  double beta = (c-d+alpha*(k-m))/(l-n);

  double cpa_xA = xA+alpha*xA_step;
  double cpa_yA = yA+alpha*yA_step;
  double cpa_zA = zA+alpha*zA_step;

  double cpa_xB = xB+beta*xB_step;
  double cpa_yB = yB+beta*yB_step;
  double cpa_zB = zB+beta*zB_step;
 
  //distance between the closest points on the lines\tracks
  //double distance = sqrt(pow((xA+alpha*xA_step)-(xB+beta*xB_step),2)+pow((yA+alpha*yA_step)-(yB+beta*yB_step),2)+pow((zA+alpha*zA_step)-(zB+beta*zB_step),2));
  double distance = sqrt(pow(cpa_xA-cpa_xB,2)+pow(cpa_yA-cpa_yB,2)+pow(cpa_zA-cpa_zB,2));

  //midpoint between points of closest approach on both lines
  double x_mid = (cpa_xA+cpa_xB)/2.;
  double y_mid = (cpa_yA+cpa_yB)/2.;
  double z_mid = (cpa_zA+cpa_zB)/2.;
 
  //check to see if this point is outside the detector
  //  TODO: < or <=, Is 0 or Lx "outside"?
  //if (x_mid > Lx || x_mid < 0 || y_mid > Ly || y_mid < 0 || z_mid > Lz || z_mid < 0) distance = -1;
  if (cpa_xA > Lx || cpa_xA < 0 || cpa_yA > Ly || cpa_yA < 0 || cpa_zA > Lz || cpa_zA < 0 ||
      cpa_xB > Lx || cpa_xB < 0 || cpa_yB > Ly || cpa_yB < 0 || cpa_zB > Lz || cpa_zB < 0) distance = -1;
 
  std::vector<double> return_vector;

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

std::vector<double> SCECalib::findDistortedClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB) const
{
  std::vector<double> return_vector;
 
  if (calibTrackA.electrons.size() < 3) {
    std::cout << "Less than three track points are inside the detector for calibTrackA." << std::endl;
    return return_vector;
  }
  if (calibTrackB.electrons.size() < 3) {
    std::cout << "Less than three track points are inside the detector for calibTrackB." << std::endl;
    return return_vector;
  }

  double min_distance = -1.;//-1 means no value has yet been set; for the first iteration
  double min_a, min_b;
  /*
  for(int iter_a = 0; iter_a < calibTrackA.track.electrons.size(); iter_a++){
    for(int iter_b = 0; iter_b < calibTrackB.track.electrons.size(); iter_b++) {
      double distance = sqrt(pow(calibTrackA.track.electrons.at(iter_a).x_mod-calibTrackB.track.electrons.at(iter_b).x_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).y_mod-calibTrackB.track.electrons.at(iter_b).y_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).z_mod-calibTrackB.track.electrons.at(iter_b).z_mod,2));
      if (min_distance == -1 || distance < min_distance) {
	min_distance = distance;
	min_a = iter_a; min_b = iter_b;
      }
    }
  }*/
  
  
  for(int iter_a = 0; iter_a < calibTrackA.electrons.size(); iter_a++){
    for(int iter_b = 0; iter_b < calibTrackB.electrons.size(); iter_b++) {
      double distance = sqrt(pow(calibTrackA.electrons.at(iter_a).x_mod-calibTrackB.electrons.at(iter_b).x_mod,2)+pow(calibTrackA.electrons.at(iter_a).y_mod-calibTrackB.electrons.at(iter_b).y_mod,2)+pow(calibTrackA.electrons.at(iter_a).z_mod-calibTrackB.electrons.at(iter_b).z_mod,2));
      if (min_distance == -1 || distance < min_distance) {
	min_distance = distance;
	min_a = iter_a; min_b = iter_b;
      }
    }
  }
    
  //TODO if min_distance is zero then no need to fit a parabola

  std::vector<elecInfo> parabola_points_calibTrackA, parabola_points_calibTrackB;
 
  if (min_a != 0) parabola_points_calibTrackA.push_back(calibTrackA.electrons.at(min_a-1));
  else parabola_points_calibTrackA.push_back(calibTrackA.electrons.at(min_a+2));
  parabola_points_calibTrackA.push_back(calibTrackA.electrons.at(min_a));
  if (min_a != calibTrackA.electrons.size() - 1) parabola_points_calibTrackA.push_back(calibTrackA.electrons.at(min_a+1));
  else parabola_points_calibTrackA.push_back(calibTrackA.electrons.at(min_a-2));
 
  if (min_b != 0) parabola_points_calibTrackB.push_back(calibTrackB.electrons.at(min_b-1));
  else parabola_points_calibTrackB.push_back(calibTrackB.electrons.at(min_b+2));
  parabola_points_calibTrackB.push_back(calibTrackB.electrons.at(min_b));
  if (min_b != calibTrackB.electrons.size()-1) parabola_points_calibTrackB.push_back(calibTrackB.electrons.at(min_b+1));
  else parabola_points_calibTrackB.push_back(calibTrackB.electrons.at(min_b-2));
 
  std::vector<double> parabolaParameters_calibTrackA = getParabolaParameters(parabola_points_calibTrackA);
  std::vector<double> parabolaParameters_calibTrackB = getParabolaParameters(parabola_points_calibTrackB);
 
  double aA = parabolaParameters_calibTrackA.at(0);
  double bA = parabolaParameters_calibTrackA.at(1);
  double dA = parabolaParameters_calibTrackA.at(2);
  double eA = parabolaParameters_calibTrackA.at(3);
  double fA = parabolaParameters_calibTrackA.at(4);
  double phiA = parabolaParameters_calibTrackA.at(5);
  double x_midA = parabolaParameters_calibTrackA.at(6);
  double y_midA = parabolaParameters_calibTrackA.at(7);
 
  double aB = parabolaParameters_calibTrackB.at(0);
  double bB = parabolaParameters_calibTrackB.at(1);
  double dB = parabolaParameters_calibTrackB.at(2);
  double eB = parabolaParameters_calibTrackB.at(3);
  double fB = parabolaParameters_calibTrackB.at(4);
  double phiB = parabolaParameters_calibTrackB.at(5);
  double x_midB = parabolaParameters_calibTrackB.at(6);
  double y_midB = parabolaParameters_calibTrackB.at(7);
 
  double stepSize_x2 = 2.5*pow(10,-4); //0.25mm
 
  double min_distance_parabola = -1.;
 
  double min_Ax = 0., min_Ay = 0., min_Az = 0., min_Bx = 0., min_By = 0., min_Bz = 0.;
  if (min_a != 0 && min_a != calibTrackA.electrons.size()-1 && min_b != 0 && min_b != calibTrackB.electrons.size()-1) {
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    /*
    double A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    double A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    double B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    double B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
    */
    
    double A_x2_0 = (calibTrackA.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    double A_x2_2 = (calibTrackA.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    double B_x2_0 = (calibTrackB.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    double B_x2_2 = (calibTrackB.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);   
    //starting and ending points for the loop //values determined below
    double A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
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
   
    for (double A_x2 = A_x2_start; A_x2 < A_x2_end;  A_x2 += stepSize_x2) {
      double A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      //get x,y,z
      double A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      double A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      double A_z = dA*A_x+eA*A_y+fA;
      for (double B_x2 = B_x2_start; B_x2 < B_x2_end; B_x2 += stepSize_x2) {
	double B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	double B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	double B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	double B_z = dB*B_x+eB*B_y+fB;
   
	double distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;
     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }
  else {
    double stepSize_x2_A = stepSize_x2;
    double stepSize_x2_B = stepSize_x2;
   
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    
    /*
    double A_x2_0 = min_a == 0 ? -1. : (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    double A_x2_2 = min_a == calibTrackA.track.electrons.size()-1 ? -1. : (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    double B_x2_0 = min_b == 0 ? -1. : (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    double B_x2_2 = min_b == calibTrackB.track.electrons.size()-1 ? -1. : (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
    */
    
        double A_x2_0 = min_a == 0 ? -1. : (calibTrackA.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    double A_x2_2 = min_a == calibTrackA.electrons.size()-1 ? -1. : (calibTrackA.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    double B_x2_0 = min_b == 0 ? -1. : (calibTrackB.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    double B_x2_2 = min_b == calibTrackB.electrons.size()-1 ? -1. : (calibTrackB.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
    int nearBoundary_A = 1, nearBoundary_B = 1;
   
    double A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
   
    //If the first point is set to be -1 then start the loop form the final point and move towards the first point until the detector boundary is reached.
    if (A_x2_0 == -1.) A_x2_start = A_x2_2;
    else if (A_x2_2 == -1.) A_x2_start = A_x2_0;
    else {
      nearBoundary_A = 0;
     
      //A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
      //A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
      
      A_x2_0 = (calibTrackA.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
      A_x2_2 = (calibTrackA.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
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
     
     // B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
     // B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
      
      B_x2_0 = (calibTrackB.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
      B_x2_2 = (calibTrackB.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
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

    double A_x = 0., A_y = 0., A_z = 0.;
    int firstIter_A = 1; //1 means it's the first iteration
    for (double A_x2 = A_x2_start; nearBoundary_A*((A_x > 0 && A_x < Lx && A_y > 0 && A_y < Ly && A_z > 0 && A_z < Lz)+firstIter_A)+!nearBoundary_A*(A_x2 < A_x2_end); A_x2 += stepSize_x2_A) {
      firstIter_A = 0;
      
      double A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      A_z = dA*A_x+eA*A_y+fA;
     
      double B_y = 0., B_x = 0., B_z = 0.;
      int firstIter_B = 1;
      for (double B_x2 = B_x2_start; nearBoundary_B*((B_x > 0 && B_x < Lx && B_y > 0 && B_y < Ly && B_z > 0 && B_z < Lz)+firstIter_B)+!nearBoundary_B*(B_x2 < B_x2_end); B_x2 += stepSize_x2_B) {
	firstIter_B = 0;
   
	double B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	B_z = dB*B_x+eB*B_y+fB;

	double distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }

  //midpoint between points of closest approach on both lines
  double x_mid = (min_Ax+min_Bx)/2.;
  double y_mid = (min_Ay+min_By)/2.;
  double z_mid = (min_Az+min_Bz)/2.;

  if(std::isnan(min_distance_parabola)) {
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

std::vector<double> SCECalib::getParabolaParameters(const std::vector<elecInfo> &parabola_points_track) const
{
  if (parabola_points_track.size() < 3) std::cout << "Less than three points provided for the parameters of parabola." << std::endl;
 
  double x_middle = parabola_points_track.at(1).x_mod;
  double y_middle = parabola_points_track.at(1).y_mod;
  double z_middle = parabola_points_track.at(1).z_mod;
 
  //first/linear transformation
  double x_0 = parabola_points_track.at(0).x_mod-x_middle;
  double x_1 = parabola_points_track.at(1).x_mod-x_middle;
  double x_2 = parabola_points_track.at(2).x_mod-x_middle;
 
  double y_0 = parabola_points_track.at(0).y_mod-y_middle;
  double y_1 = parabola_points_track.at(1).y_mod-y_middle;
  double y_2 = parabola_points_track.at(2).y_mod-y_middle;
 
  double z_0 = parabola_points_track.at(0).z_mod;
  double z_1 = parabola_points_track.at(1).z_mod;
  double z_2 = parabola_points_track.at(2).z_mod;
 
  //angle of rotation where y_0_2 == y_2_2
  double phi = atan((y_2-y_0)/(x_0-x_2));
 
  //preform the second transformation
  double x_0_2 = x_0*cos(phi)-y_0*sin(phi);
  double x_1_2 = x_1*cos(phi)-y_1*sin(phi);
  double x_2_2 = x_2*cos(phi)-y_2*sin(phi);
 
  double y_0_2 = x_0*sin(phi)+y_0*cos(phi);
  double y_1_2 = x_1*sin(phi)+y_1*cos(phi);
  double y_2_2 = x_2*sin(phi)+y_2*cos(phi);
 
  //since x_1_2 = 0 and y_1_2 = 0, because that's the middle point, c = 0
  double a = ((y_0_2*x_2_2)-(y_2_2*x_0_2))/(x_0_2*x_2_2*(x_0_2-x_2_2));
  double b = (y_0_2-(a*pow(x_0_2,2)))/x_0_2;
 
  //find the plane defined by the three points without any transformations: z = d*x + e*y + f
  //any one of the three points dotted with the normal vector will give the plane: P_1*n = 0
 
  //comps. for a vector connecting point 2 to point 1
  double p1p2_x = parabola_points_track.at(1).x_mod-parabola_points_track.at(0).x_mod;
  double p1p2_y = parabola_points_track.at(1).y_mod-parabola_points_track.at(0).y_mod;
  double p1p2_z = parabola_points_track.at(1).z_mod-parabola_points_track.at(0).z_mod;
 
  //comps. for a vector connecting point 3 to point 1
  double p1p3_x = parabola_points_track.at(2).x_mod-parabola_points_track.at(0).x_mod;
  double p1p3_y = parabola_points_track.at(2).y_mod-parabola_points_track.at(0).y_mod;
  double p1p3_z = parabola_points_track.at(2).z_mod-parabola_points_track.at(0).z_mod;
 
  // normal vector, n, = P1P2xP1P3
  double norm_x = (p1p2_y*p1p3_z)-(p1p3_y*p1p2_z);
  double norm_y = (p1p3_x*p1p2_z)-(p1p2_x*p1p3_z);
  double norm_z = (p1p2_x*p1p3_y)-(p1p3_x*p1p2_y);
 
  // z = d*x + e*y + f
  double d = -norm_x/norm_z;
  double e = -norm_y/norm_z;
  double f = ((parabola_points_track.at(0).x_mod*norm_x)+(parabola_points_track.at(0).y_mod*norm_y)+(parabola_points_track.at(0).z_mod*norm_z))/norm_z;
 
  std::vector<double> return_vector;
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

std::vector<distortionMap> SCECalib::doCalibFaces(const std::vector<calibTrackInfo> &cosmicCalibTracks, int minTrackPoints, int numTrackSegPoints)
{


  calibTrackInfo calibTrack;
  int numCosmicTracks = cosmicCalibTracks.size();
  
  distortionMap deltaXMap;
  distortionMap deltaYMap;
  distortionMap deltaZMap;

  double xVal;
  double yVal;
  double zVal;
  double xValDistorted;
  double yValDistorted;
  double zValDistorted;
  int xCalibLowIndex;
  int xCalibHighIndex;
  int yCalibLowIndex;
  int yCalibHighIndex;
  int zCalibLowIndex;
  int zCalibHighIndex;
  double xCalibFrac;
  double yCalibFrac;
  double zCalibFrac;
  double tempFactor;

  int trackNum;
  int numElec;
  
  for(int i = 0; i < numCosmicTracks; i++)
  {
    std::cout << "COSMIC-FACE " << i << std::endl;

    calibTrack = cosmicCalibTracks.at(i);
    if(calibTrack.electrons.size() < minTrackPoints) continue;

    trackNum = i;
    numElec = calibTrack.electrons.size();
    std::vector<elecInfo> faceCalibPoints = calibTrack.electrons;
    

    int numNearAnode = 0;
    if (faceCalibPoints[0].x_mod > (Lx - maxXdist)) {
      numNearAnode++;
    }
    if (faceCalibPoints[numElec-1].x_mod > (Lx - maxXdist)) {
      numNearAnode++;
    }
    if((numNearAnode == 0) || (numNearAnode == 2)) continue;
	
    std::vector<elecInfo> calibPoints;
    
    double avgX1 = 0.0;
    double avgX2 = 0.0;

    for(int j = 0; j < 5; j++)
    {
      double elecX = faceCalibPoints[j].x_mod;
      avgX1 += elecX;
    }
    avgX1 /= 5.0;

    for(int j = numElec-1; j > numElec-6; j--)
    {
      double elecX = faceCalibPoints[j].x_mod;
      avgX2 += elecX;
    }
    avgX2 /= 5.0;

    int numBadPoints = 0;
    if(avgX1 < avgX2) {
      xValDistorted = faceCalibPoints[0].x_mod;
      yValDistorted = faceCalibPoints[0].y_mod;
      zValDistorted = faceCalibPoints[0].z_mod;

      for (int j = numElec-1; j > numElec-1-numTrackSegPoints; j--) {
        if ((faceCalibPoints[j].x_mod <= 0.0) || (faceCalibPoints[j].x_mod >= Lx) || (faceCalibPoints[j].y_mod <= 0.0) || (faceCalibPoints[j].y_mod >= Ly) || (faceCalibPoints[j].z_mod <= 0.0) || (faceCalibPoints[j].z_mod >= Lz)) {
          numBadPoints++;
	  continue;
	}
        std::vector<double> pointXYZ;
	pointXYZ.push_back( faceCalibPoints[j].x_mod );
	pointXYZ.push_back( faceCalibPoints[j].y_mod );
	pointXYZ.push_back( faceCalibPoints[j].z_mod );
	
	elecInfo tempPoint;
        tempPoint.x = faceCalibPoints[j].x_mod + getTruthOffset(pointXYZ,axisType::xAxis);
        tempPoint.y = faceCalibPoints[j].y_mod + getTruthOffset(pointXYZ,axisType::yAxis);
        tempPoint.z = faceCalibPoints[j].z_mod + getTruthOffset(pointXYZ,axisType::zAxis);

        calibPoints.push_back(tempPoint);
      }
    }
    
    else {
      xValDistorted = faceCalibPoints[numElec-1].x_mod;
      yValDistorted = faceCalibPoints[numElec-1].y_mod;
      zValDistorted = faceCalibPoints[numElec-1].z_mod;

      for (int j = 0; j < numTrackSegPoints; j++) {
	if ((faceCalibPoints[j].x_mod <= 0.0) || (faceCalibPoints[j].x_mod >= Lx) || (faceCalibPoints[j].y_mod <= 0.0) || (faceCalibPoints[j].y_mod >= Ly) || (faceCalibPoints[j].z_mod <= 0.0) || (faceCalibPoints[j].z_mod >= Lz)) {
          numBadPoints++;
	  continue;
	}
	
	elecInfo tempPoint;
	std::vector<double> pointXYZ;
	pointXYZ.push_back( faceCalibPoints[j].x_mod );
	pointXYZ.push_back( faceCalibPoints[j].y_mod );
	pointXYZ.push_back( faceCalibPoints[j].z_mod );
	
        tempPoint.x = faceCalibPoints[j].x_mod + getTruthOffset(pointXYZ,axisType::xAxis);
        tempPoint.y = faceCalibPoints[j].y_mod + getTruthOffset(pointXYZ,axisType::yAxis);
        tempPoint.z = faceCalibPoints[j].z_mod + getTruthOffset(pointXYZ,axisType::zAxis);

        calibPoints.push_back(tempPoint);
      }
    }
    if (numBadPoints > 5) continue;
    
    PCAResult results;
    results.doPCA(calibPoints);
    std::cout << "PCA Done." << std::endl;

    double startX = (results.getStartPoints() ).at(0);
    double startY = (results.getStartPoints() ).at(1);
    double startZ = (results.getStartPoints() ).at(2);
    

    double unitX = (results.getUnits() ).at(0);
    double unitY = (results.getUnits() ).at(1);
    double unitZ = (results.getUnits() ).at(2);
    
    int whichFace = -1;
    
    double t_Top = (Ly-startY)/unitY;
    double x_Top = unitX*t_Top + startX;
    double y_Top = Ly;
    double z_Top = unitZ*t_Top + startZ;
    if((t_Top > 0.0) && (x_Top > 0.0) && (x_Top < Lx) && (z_Top > 0.0) && (z_Top < Lz)) {
      xVal = x_Top;
      yVal = y_Top;
      zVal = z_Top;
      whichFace = 0;
    }

    double t_Bottom = (0.0-startY)/unitY;
    double x_Bottom = unitX*t_Bottom + startX;
    double y_Bottom = 0.0;
    double z_Bottom = unitZ*t_Bottom + startZ;
    if((t_Bottom > 0.0) && (x_Bottom > 0.0) && (x_Bottom < Lx) && (z_Bottom > 0.0) && (z_Bottom < Lz)) {
      xVal = x_Bottom;
      yVal = y_Bottom;
      zVal = z_Bottom;
      whichFace = 1;
    }

    double t_Upstream = (0.0-startZ)/unitZ;
    double x_Upstream = unitX*t_Upstream + startX;
    double y_Upstream = unitY*t_Upstream + startY;
    double z_Upstream = 0.0;
    if((t_Upstream > 0.0) && (x_Upstream > 0.0) && (x_Upstream < Lx) && (y_Upstream > 0.0) && (y_Upstream < Ly)) {
      xVal = x_Upstream;
      yVal = y_Upstream;
      zVal = z_Upstream;
      whichFace = 2;
    }

    double t_Downstream = (Lz-startZ)/unitZ;
    double x_Downstream = unitX*t_Downstream + startX;
    double y_Downstream = unitY*t_Downstream + startY;
    double z_Downstream = Lz;
    if((t_Downstream > 0.0) && (x_Downstream > 0.0) && (x_Downstream < Lx) && (y_Downstream > 0.0) && (y_Downstream < Ly)) {
      xVal = x_Downstream;
      yVal = y_Downstream;
      zVal = z_Downstream;
      whichFace = 3;
    }

    double t_Cathode = (0.0-startX)/unitX;
    double x_Cathode = 0.0;
    double y_Cathode = unitY*t_Cathode + startY;
    double z_Cathode = unitZ*t_Cathode + startZ;
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
    
    xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((double) xCalibLowIndex);
    yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((double) yCalibLowIndex);
    zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((double) zCalibLowIndex);
    
    distortionVoxel *vox;
    
    //TOP
    if (whichFace == 0) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
      

      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);      
      

      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      
      vox = deltaXMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, nCalibDivisions_y, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
    }
    
    //BOTTOM
    else if (whichFace == 1) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      
      vox = deltaXMap.getVoxel(xCalibHighIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, 0, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      
      vox = deltaXMap.getVoxel(xCalibHighIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, 0, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
    }
    
    //UPSTREAM
    else if (whichFace == 2) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibLowIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibHighIndex, 0);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
    }
    
    //DOWNSTREAM
    else if (whichFace == 3) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      
      vox = deltaXMap.getVoxel(xCalibLowIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibLowIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibLowIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
	
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibLowIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      
      vox = deltaXMap.getVoxel(xCalibHighIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(xCalibHighIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(xCalibHighIndex, yCalibHighIndex, nCalibDivisions_z);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
    }
    
    //CATHODE
    else if (whichFace == 4) {
      tempFactor = (1.0-yCalibFrac)*(1.0-zCalibFrac);
      
      vox = deltaXMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = (1.0-yCalibFrac)*zCalibFrac;
      
      vox = deltaXMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(nCalibDivisions_x, yCalibLowIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = yCalibFrac*(1.0-zCalibFrac);
      
      vox = deltaXMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibLowIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);

      tempFactor = yCalibFrac*zCalibFrac;
      
      vox = deltaXMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(xVal-xValDistorted);
      
      vox = deltaYMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(yVal-yValDistorted);
      
      vox = deltaZMap.getVoxel(nCalibDivisions_x, yCalibHighIndex, zCalibHighIndex);
      vox->addWeight(tempFactor);
      vox->addDistortion(zVal-zValDistorted);
    }

    if (whichFace >= 0) {
      std::cout << "  " << whichFace << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << std::endl;
    }
  }    
  
  std::vector<distortionMap> returnVec;
  returnVec.push_back(deltaXMap);
  returnVec.push_back(deltaYMap);
  returnVec.push_back(deltaZMap);
  
  return returnVec;
}

int SCECalib::calculate2DMaps(std::vector<distortionMap> faceCalibDistortions){
    xTopMap = faceCalibDistortions[0].calculateMap(0);
    yTopMap = faceCalibDistortions[1].calculateMap(0);
    zTopMap = faceCalibDistortions[2].calculateMap(0);
    
    xBottomMap = faceCalibDistortions[0].calculateMap(1);
    yBottomMap = faceCalibDistortions[1].calculateMap(1);
    zBottomMap = faceCalibDistortions[2].calculateMap(1);
    
    xUpstreamMap = faceCalibDistortions[0].calculateMap(2);
    yUpstreamMap = faceCalibDistortions[1].calculateMap(2);
    zUpstreamMap = faceCalibDistortions[2].calculateMap(2);
    
    xDownstreamMap = faceCalibDistortions[0].calculateMap(3);
    yDownstreamMap = faceCalibDistortions[1].calculateMap(3);
    zDownstreamMap = faceCalibDistortions[2].calculateMap(3);
    
    xCathodeMap = faceCalibDistortions[0].calculateMap(4);
    yCathodeMap = faceCalibDistortions[1].calculateMap(4);
    zCathodeMap = faceCalibDistortions[2].calculateMap(4);
    
    return 0;

}

int main(int argc, char *argv[]){
   
   calibSteps runType; 
   if(argc == 1)   
    runType = calibSteps::fullCalib;
   
   else if (argc == 2){
    int run = atoi(argv[1]);
    if(run < 0 || run > 3){
      std::cerr << "FATAL: Calibration modes can only be between 0 (face only) and 3 (fullCalib). You entered " << run << std::endl;
      return -1;  
    
    }
    runType = static_cast<calibSteps>(run);
   }
   
   else{
     std::cerr << "FATAL: Unrecongnized arguments. Aborting." << std::endl;
      return -1;
   
   } 
   bool faceCalibrated = false;
   
   SCECalib *calibrator = new SCECalib();   
         
   std::vector<trackInfo> laserTracks  = calibrator->getLaserTrackSet();
   std::vector<trackInfo> cosmicTracks = calibrator->getCosmicTrackSet(faceCalibrated);
   
   std::vector<calibTrackInfo> cosmicCalibTracks = calibrator->makeCalibTracks(cosmicTracks);
   
   std::cout << "Number of Tracks: " << cosmicCalibTracks.size() << std::endl;
  
  switch (runType)  {
   case(bulkOnly) : {
    std::cout << "Running Bulk Only Calibration." << std::endl;                                  
    
    std::vector<distortionMap> calibDistortions = calibrator->doCosmicCosmicCalib(cosmicCalibTracks);
    
  
    std::vector< std::vector < std::vector<float> > > xMap     =  calibDistortions[0].calculateMap();
    std::vector< std::vector < std::vector<float> > > xMapErrs =  calibDistortions[0].calculateMapErrors();
  
    std::vector< std::vector < std::vector<float> > > yMap     =  calibDistortions[1].calculateMap();
    std::vector< std::vector < std::vector<float> > > yMapErrs =  calibDistortions[1].calculateMapErrors();
  
    std::vector< std::vector < std::vector<float> > > zMap     =  calibDistortions[2].calculateMap();
    std::vector< std::vector < std::vector<float> > > zMapErrs =  calibDistortions[2].calculateMapErrors();
  
    return 0;
   }
   
   case(faceOnly) : {
     std::cout << "Running Face Only Calibration." << std::endl;
     std::vector<distortionMap> faceCalibDistortions = calibrator->doCalibFaces(cosmicCalibTracks, 50, 15);
     //std::cout << "Calculating 2D Maps." << std::endl;
     //calibrator->calculate2DMaps(faceCalibDistortions);
     faceCalibrated = true;
     return 0;
   
   }
   
   case(bulkAndFace) : {
     std::vector<distortionMap> calibDistortions = calibrator->doCosmicCosmicCalib(cosmicCalibTracks);
     calibrator->loadTruthMap(calibDistortions);     
     std::vector<distortionMap> faceCalibDistortions = calibrator->doCalibFaces(cosmicCalibTracks, 50, 15);
     calibrator->calculate2DMaps(faceCalibDistortions);
     faceCalibrated = true;
     
     return 0;
     
     
   
   }
   
   case(fullCalib) : {
     std::vector<distortionMap> calibDistortions = calibrator->doCosmicCosmicCalib(cosmicCalibTracks);
     calibrator->loadTruthMap(calibDistortions);     
     std::vector<distortionMap> faceCalibDistortions = calibrator->doCalibFaces(cosmicCalibTracks, 50, 15);
     calibrator->calculate2DMaps(faceCalibDistortions);
     faceCalibrated = true;
     calibrator->getCosmicTrackSet(faceCalibrated);              
     std::vector<distortionMap> finalDistortions = calibrator->doCosmicCosmicCalib(cosmicCalibTracks);
     return 0;
     
   
   }
   
   default : return 1;
  }
  
  //probably shouldn't get here 
  return -1;
   
   

}



