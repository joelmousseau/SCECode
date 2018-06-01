#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include "/uboone/app/users/joelam/SpaceChargeStudy/CalibSCE/Eigen/Dense.h"

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

enum axisType{
        xAxis,
	yAxis,
	zAxis
     
     };
     
enum calibSteps{
       bulkOnly,
       faceOnly,
       bulkAndFace,
       fullCalib          
     };

///Points along a track
/** Class for containing individual points along a track. Includes modified and unmodified points, as well as a method to set modified points */
class elecInfo{
  
  public:
    //!Default Constructor
    elecInfo() : x(-1.0), y(-1.0), z(-1.0), t(-1.0), x_mod(0.0), y_mod(0.0), z_mod(0.0), t_mod(0.0), fate(0.0)  {}; 
    
    //! Updates track points
    /*! 
         \param s vector of updated (x,y,z) points
	 \return nothing
    
    */
    void updateSMod(std::vector<double> s )  {x_mod = s[0]; y_mod = s[1]; z_mod = s[2];}; 
    
    double x; //!< X position before truth modification
    double y; //!< Y position before truth modification
    double z; //!< Z position before truth modification
    double t; //!< T position before truth modification
    double x_mod; //!< X position after truth modification
    double y_mod; //!< Y position after truth modification
    double z_mod; //!< Z position after truth modification
    double t_mod; //!< T position after truth modification
    int    fate; //!< 1 If the point is used in calibration. Currently not implemented

};

///Track without calibration info
/** Defines a track, and carries the track points (electrons) */
class trackInfo {
  
  public:
    //! Default Constructor
    trackInfo() : energy(1000000.0), pdgID(22), x0(0.0), y0(0.0), z0(0.0), x1(0.0), y1(0.0), z1(0.0), theta(0.0), phi(0.0)  {};
    //! Update Start and End Points
    /*! \param s vector of (x,y,z) start points
         \param e vector of (x,y,z) end points
	 \return nothing 
    */ 	 
    void updateEndpoints(std::vector<double> s, std::vector<double> e){x0 = s[0]; y0 = s[1]; z0 = s[2]; x1 = e[0]; y1 = e[1]; z1 = e[2]; }
    
    //! Update Phi and theta of a track
    /*! \param t track theta
          \param p track phi
	  \return nothing
      */	  
    void updateAngles(double t, double p){theta = t; phi = p;}
    
    bool isNotContained(std::vector<double> max, std::vector<double> L) {return (((x0 < max[0]) && (x1 < max[0])) || ((x0 > (L[0] - max[0])) && (x1 > (L[0] - max[0]))) || ((x0 < max[0]) && (y0 < max[1])) || ((x0 < max[0]) && (y0 > (L[1] - max[1]))) || ((x0 < max[0]) && (z0 < max[2])) || ((x0 < max[0]) && (z0 > (L[2] -max[2]))) || ((x1 < max[0]) && (y1 < max[1])) || ((x1 < max[0]) && (y1 > (L[1] - max[1]))) || ((x1 < max[0]) && (z1 < max[2])) || ((x1 < max[0]) && (z1 > (L[2] -max[2]))) || ((y0 < max[1]) && (z0 < max[2])) || ((y0 > (L[1] - max[1])) && (z0 < max[2])) || ((y0 < max[1]) && (z0 > (L[2] - max[2]))) || ((y0 > (L[1] - max[1])) && (z0 > (L[2] - max[2]))) || ((y1 < max[1]) && (z1 < max[2])) || ((y1 > (L[1] - max[1])) && (z1 < max[2])) || ((y1 < max[1]) && (z1 > (L[2] - max[2]))) || ((y1 > (L[1] - max[1])) && (z1 > (L[2] - max[2])))); }
    
    bool isSomething(std::vector<double> max, std::vector<double> L){return (((x0 > max[0]) && (x0 < (L[0] - max[0])) && (y0 > max[1]) && (y0 < (L[1] - max[1])) && (z0 > max[2]) && (z0 < (L[2] - max[2]))) || ((x1 > max[0]) && (x1 < (L[0] - max[0])) && (y1 > max[1]) && (y1 < (L[1] - max[1])) && (z1 > max[2]) && (z1 < (L[2] - max[2]))));}
    
    bool isCathodeExiting(std::vector<double> max, std::vector<double> L){return (((x0 > (L[0] - max[0])) && (x1 > max[0])) || ((x1 > (L[0] - max[0])) && (x0 > max[0]))); }
    
    bool isAnodeExiting(std::vector<double> max, std::vector<double> L){return (((x0 < (L[0] - max[0])) && (x1 < max[0])) || ((x1 < (L[0] - max[0])) && (x0 < max[0])));}
    
    bool isCrossing(std::vector<double> max, std::vector<double> L){return (((x0 > (L[0] - max[0])) && (x1 < max[0])) || ((x1 > (L[0] - max[0])) && (x0 < max[0])));}
    
    
    
    
  
    int pdgID; //!< PDG code of track
    double energy; //!< Energy of Track (MeV) 
    double x0; //!< Track X start point
    double y0; //!< Track Y start point
    double z0; //!< Track Z start point
    double x1; //!< Track X end point
    double y1; //!< Track Y end point
    double z1; //!< Track Z end point
    double theta; //!< Track theta
    double phi; //!< Track phi
    std::vector<elecInfo> electrons; //!< Container of points (electrons) that make up the track

};

///Track with calibration info
/** Defines a calibrated track. Contains the calibration factors as D{x,y,z}Vec*/ 
class calibTrackInfo : public trackInfo {
  
  public:
    //!Default constructor
    calibTrackInfo() {};
    
    //!Constructor using an uncalibrated track as input 
    /*! \param track Uncalibrated track. Calibrated track aquires uncalibrated tracks points
      */
    calibTrackInfo(trackInfo track);
                 
    std::vector<double> DxVec; //!< Vector of X deflections
    std::vector<double> DyVec; //!< Vector of Y deflections
    std::vector<double> DzVec; //!< Vector of Z deflections
    bool calibFlag; //!< True if track has been calibrated
};

///Encapsulation of distortions and weights within a voxel
/** Three dimensional voxel with its' distortion and weights. Voxels are numbered as a three-tuple of x,y,z. Numbering is not guaranteed to be unique, it's your job to make sure they are uniquely numbered*/ 
class  distortionVoxel{
   
   public:
     //! Default Constructor
     distortionVoxel();
     
     //! Constructor taking x, y, z voxel number
     /*! \param x X voxel number
           \param y Y voxel number
	   \param z Z voxel number
	*/   
     distortionVoxel(int x, int y, int z);
     
     ~distortionVoxel();
     //! Const gettor
     /*! \return (x,y,z) vector of voxel number
     */
     std::vector<int> getVoxelNumber() const;
     
     //! Set voxel number
     /*! \param index (x,y,z) vector of voxel number
          \return 0 if sucess
	*/  
     int setVoxelNumber(std::vector<int> index);
     
     //! Add a distortion to the vector of distortions
     /*! \param dist distortion to add
          \return 0 if success
	*/  
     int addDistortion(float dist);
     
     //! Clear the vector of distortions
     int clearDistortions();
     
     //! Add a weight to the vector of weights
     /*! \param dist weight to add
          \return 0 if success
       */	  
     int addWeight(float dist);
     
     //! Clear the vector of weights
     int clearWeights();
     
     //! Const gettor
     std::vector<float> getDistortions() const;
     
     //! Const gettor
     std::vector<double> getWeights() const;     
     
     //! Const gettor
     unsigned int getNPairs() const;
     
     /* \param vox distortionVoxel to calculate
          \return median of the distortions
	*/ 
     float calculateMedian() const;
     
     /*float calculateMode()*/;
     
     //! Calculate the weighted mean of all distortions in a voxel
     /* \param vox distortionVoxel to calculate
          \return mean of the distortions
	*/
     float calculateMean() const;
     
     //! Calculate the RMS of each voxel based on rms / number of entries. 
     /* \param vox distortionVoxel to caclulate
           \return uncertainty of distortions
       */
     float calculateRMS() const;
     
     //! Calculate the sum of weights of a particular voxel
     /*
         Some of ROOTs math functions will throw an error if the sum of weights == 0
	 This enables us to work around these errors by checking the weight sum ahead of time
	 \param vox distortionVoxel to calculate
	 \retrun sum of the weights in that voxel 
     */
     float calculateWeightSum() const;
   
   private:
     std::vector<float> distortions; //!< Vector of distortions. 
     std::vector<double> weights;    //!< Vector of weights. Double because ROOT is horrible.
     int xVoxel;                    //!< X index of voxel
     int yVoxel;                    //!< Y index of voxel
     int zVoxel;                    //!< Z index of voxel
               
};

///Container Class of voxels
/** Intializes a 3D array of voxels of fixed size. Contains methods to get pointers to voxels, and compute the calibration and error of a specific voxel based on mean or median*/ 
class distortionMap
{
   public:
     //!Default Constructor
     distortionMap();
     
     ~distortionMap();
     
     //!Get the 3D array of voxels. Can be slow.
     /*! \return 3D array of voxels
       */     
     std::vector<  std::vector < std::vector<distortionVoxel> > > getMap();
     
     //! Get a pointer to a specific voxel
     /* \param x X index
          \param y Y index
	  \param z Z index
	  \return pointer to a distortionVoxel
      */	  
     distortionVoxel * getVoxel(int x, int y, int z);  
     
     //! Calculate a calibration map based on the distortions within each voxels
     std::vector< std::vector < std::vector<float> > > calculateMap();
     
     //! Calculate the errors on a map based on the number of track crossings per voxel 
     std::vector< std::vector < std::vector<float> > > calculateMapErrors();
     
     //! Calculate a face calibration map based on the distortions within each voxels
     std::vector < std::vector<float> >  calculateMap(int face);
     
     //! Calculate the errors on a face map based on the number of track crossings per voxel 
     std::vector < std::vector<float> >  calculateMapErrors(int face);         
     
     private:
       int FaceToVoxel(int face) const; //!< Map face index to a voxel number
       
       std::vector<  std::vector < std::vector<distortionVoxel> > > theMap; //!< 3D array of voxels. Do not access directly! Use getVoxel method!
       
       
       int nXVoxels = 26;  //!< Number of X voxels, geometry dependent!
       int nYVoxels = 26;  //!< Number of Y voxels, geometry dependent!
       int nZVoxels = 101; //!< Number of Z voxels, geometry dependent!
       std::string statMode = "mean"; //!< Use mean or median for calculateMap()
       
               

};




class PCAResult {
   public:
      PCAResult() {} ;
     
     ~PCAResult() {} ;
     
     TVector3 centroid;
     std::pair<TVector3,TVector3> endPoints;
     float length;
     TVector3 eVals;
     std::vector<TVector3> eVecs;
     
     void doPCA(const std::vector<elecInfo> &points);



};

class SCECalib{
   
   public:
     SCECalib() {Intialize(); loadTruthMap(); loadTruthMap(true);}
     
     ~SCECalib(){outputFile->Write(); outputFile->Close();}
               
     std::vector<trackInfo> getLaserTrackSet();
     
     std::vector<trackInfo> getCosmicTrackSet(bool isCalibrated);
     
     double doCoordTransformX(double inputX) const;

     double doCoordTransformY(double inputY) const;

     double doCoordTransformZ(double inputZ) const;
     
     std::vector<distortionMap>  doCosmicCosmicCalib(const std::vector<calibTrackInfo> &cosmicCalibTracks);
     
     std::vector<calibTrackInfo> makeCalibTracks(const std::vector<trackInfo> &tracks);
     
     //int doCalibration(const std::vector<trackInfo> &laserTracks, const std::vector<trackInfo> &cosmicTracks);
     
     std::vector<double> findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB) const;
     
     std::vector<double> findDistortedClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB) const;
     
     std::vector<double> getParabolaParameters(const std::vector<elecInfo> &parabola_points_track) const;
     
     std::vector<distortionMap> doCalibFaces(const std::vector<calibTrackInfo> &cosmicTracks, int minTrackPoints, int numTrackSegPoints);
     
     int setStepsToRun(calibSteps steps) {stepsToRun = steps; return 1;}
     
     int calculate2DMaps(std::vector<distortionMap> faceCalibDistortions);
     
     void loadTruthMap(std::vector<distortionMap>& inputMaps);
     

   private:
     void Intialize();
     
     double getTruthOffset(std::vector<double> sVec, axisType comp, bool isFwd = false);
     
     double getCalibOffset(std::vector<double> sVec, axisType comp, int face);
     
     void loadTruthMap(bool isFwd = false);          
     
     int cosmicTruthMode;          
     
     std::string inputFileLaser;
     std::string inputFileCosmic;
     
     
     
     calibSteps stepsToRun;

     TFile *outputFile;

     double Lx;
     double Ly;
     double Lz;


     bool isMC;
     bool doBulk;
     bool faceCalibrated;

     double relAngleCut;
     double maxXdist;
     double maxYdist;
     double maxZdist;
     double maxDistFactor;
     double distScale;
     
     double SCEFactor;

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
     
     distortionMap xTruthMap;
     distortionMap yTruthMap; 
     distortionMap zTruthMap;
     
     distortionMap xFwdTruthMap;
     distortionMap yFwdTruthMap; 
     distortionMap zFwdTruthMap;
     
     std::vector< std::vector<float> > xTopMap;
     std::vector< std::vector<float> > yTopMap;
     std::vector< std::vector<float> > zTopMap;
     
     std::vector< std::vector<float> > xBottomMap;
     std::vector< std::vector<float> > yBottomMap;
     std::vector< std::vector<float> > zBottomMap;
     
     std::vector< std::vector<float> > xUpstreamMap;
     std::vector< std::vector<float> > yUpstreamMap;
     std::vector< std::vector<float> > zUpstreamMap;
     
     std::vector< std::vector<float> > xDownstreamMap;
     std::vector< std::vector<float> > yDownstreamMap;
     std::vector< std::vector<float> > zDownstreamMap;
     
     std::vector< std::vector<float> > xCathodeMap;
     std::vector< std::vector<float> > yCathodeMap;
     std::vector< std::vector<float> > zCathodeMap;
     
     
     
     
            
     unsigned int randSeed;
     
     

};


