#!/usr/bin/python

import numpy as np
import scipy
from scipy.spatial import Delaunay
import scipy.interpolate.interpnd as interpnd
import ROOT
from ROOT import TH1, TH3, TFile, TAxis, gROOT, TCanvas
from array import array

Universes = 1

#I can't figure out how to get an array of bins from pyroot, so I need to hard code them for now
xBins = [ -5.1208000, 5.1208000, 15.362400, 25.604000, 35.845600, 46.087200, 56.328800, 66.570400, 76.812000, 87.053600, 97.295200, 107.53680, 117.77840, 128.02000, 138.26160, 148.50320, 158.74480, 168.98640, 179.22800, 189.46960, 199.71120, 209.95280, 220.19440, 230.43600, 240.67760, 250.91920, 261.0 ]
yBins = [ -120.90000, -111.60000, -102.30000, -93.000000, -83.700000, -74.400000, -65.100000, -55.800000, -46.500000, -37.200000, -27.900000, -18.600000, -9.3000000, 0.0000000, 9.3000000, 18.600000, 27.900000, 37.200000, 46.500000, 55.800000, 65.100000, 74.400000, 83.700000, 93.000000, 102.30000, 111.60000, 120.9 ]
zBins = [-5.1840000, 5.1840000, 15.552000, 25.920000, 36.288000, 46.656000, 57.024000, 67.392000, 77.760000, 88.128000, 98.496000, 108.86400, 119.23200, 129.60000, 139.96800, 150.33600, 160.70400, 171.07200, 181.44000, 191.80800, 202.17600, 212.54400, 222.91200, 233.28000, 243.64800, 254.01600, 264.38400, 274.75200, 285.12000, 295.48800, 305.85600, 316.22400, 326.59200, 336.96000, 347.32800, 357.69600, 368.06400, 378.43200, 388.80000, 399.16800, 409.53600, 419.90400, 430.27200, 440.64000, 451.00800, 461.37600, 471.74400, 482.11200, 492.48000, 502.84800, 513.21600, 523.58400, 533.95200, 544.32000, 554.68800, 565.05600, 575.42400, 585.79200, 596.16000, 606.52800, 616.89600, 627.26400, 637.63200, 648.00000, 658.36800, 668.73600, 679.10400, 689.47200, 699.84000, 710.20800, 720.57600, 730.94400, 741.31200, 751.68000, 762.04800, 772.41600, 782.78400, 793.15200, 803.52000, 813.88800, 824.25600, 834.62400, 844.99200, 855.36000, 865.72800, 876.09600, 886.46400, 896.83200, 907.20000, 917.56800, 927.93600, 938.30400, 948.67200, 959.04000, 969.40800, 979.77600, 990.14400, 1000.5120, 1010.8800, 1021.2480, 1031.6160, 1037 ]


#gridXFile = open("MergedFwdMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed.csv")
inputMapFile =  ROOT.TFile.Open("VariedBackwardMapsVolumedSmoothed.root")
histogramFile = ROOT.TFile.Open("VariedForwardMapsVolumedSmoothed.root", "RECREATE")


for uni in range(Universes):
  print "Universe: %d" % uni
  inputMapFile.cd()
  if(Universes == 1):
    xHistName = "Reco_Displacement_X"
    yHistName = "Reco_Displacement_Y"
    zHistName = "Reco_Displacement_Z"

  else:
    xHistName = "Reco_Displacement_X_%d" % uni
    yHistName = "Reco_Displacement_Y_%d" % uni
    zHistName = "Reco_Displacement_Z_%d" % uni

  h_combinedDX = inputMapFile.Get(xHistName)
  h_combinedDY = inputMapFile.Get(yHistName)
  h_combinedDZ = inputMapFile.Get(zHistName)

  nBinsX = h_combinedDX.GetNbinsX()
  nBinsY = h_combinedDX.GetNbinsY()
  nBinsZ = h_combinedDX.GetNbinsZ()

  x_reco = 0.0
  y_reco = 0.0
  z_reco = 0.0


  points = []
  grid = []
  dX = []
  dY = []
  dZ = []

  for i in range(nBinsX):
    if(i != nBinsX):
        x_reco = (xBins[i] + xBins[i+1])/2.0
    else:
        x_reco = 256.0

    for j in range(nBinsY):
        if(j != nBinsY):
           y_reco = (yBins[j] + yBins[j+1])/2.0
        else:
           y_reco = 116.5
        for k in range(nBinsZ):
          if(k != nBinsZ):
               z_reco = (zBins[k] + zBins[k+1])/2.0
          else:
               z_reco = 1037.0

          grid.append((x_reco, y_reco, z_reco))
          xDef = h_combinedDX.GetBinContent(i+1,j+1,k+1)
          yDef = h_combinedDY.GetBinContent(i+1,j+1,k+1)
          zDef = h_combinedDZ.GetBinContent(i+1,j+1,k+1)
#print "%f %f %f %f %f %f" % (x_reco, y_reco, z_reco, xDef, yDef, zDef)
          points.append((x_reco+xDef, y_reco+yDef, z_reco+zDef))
          dX.append(-1.0*xDef)
          dY.append(-1.0*yDef)
          dZ.append(-1.0*zDef)

  pointCloud  = np.asarray(points, dtype = float)
  regularGrid = np.asarray(grid,   dtype = float)

  tri = Delaunay(pointCloud)

  regularDX = interpnd.LinearNDInterpolator(tri,dX, 10e10)(regularGrid)
  regularDY = interpnd.LinearNDInterpolator(tri,dY, 10e10)(regularGrid)
  regularDZ = interpnd.LinearNDInterpolator(tri,dZ, 10e10)(regularGrid)

  print array('d',xBins)[0]

#histogramFile = ROOT.TFile.Open("Forward.root", "RECREATE")
  histogramFile.cd()
  if(Universes == 1):
    outHistNameX = "combined_fwd_dX"
    outHistNameY = "combined_fwd_dY"
    outHistNameZ = "combined_fwd_dZ"
  else:
    outHistNameX = "combined_fwd_dX_%d" % uni
    outHistNameY = "combined_fwd_dY_%d" % uni
    outHistNameZ = "combined_fwd_dZ_%d" % uni
  
  hist_fwdDX = ROOT.TH3F(outHistNameX, "", 26, -5.1208, 261.161, 26, -120.9, 120.9, 101, -5.184, 1041.98 )
  hist_fwdDY = ROOT.TH3F(outHistNameY, "", 26, -5.1208, 261.161, 26, -120.9, 120.9, 101, -5.184, 1041.98 )
  hist_fwdDZ = ROOT.TH3F(outHistNameZ, "", 26, -5.1208, 261.161, 26, -120.9, 120.9, 101, -5.184, 1041.98 )


  for index in range(len(regularDX)):
     hist_fwdDX.Fill(regularGrid[index][0], regularGrid[index][1], regularGrid[index][2], regularDX[index])
     hist_fwdDY.Fill(regularGrid[index][0], regularGrid[index][1], regularGrid[index][2], regularDY[index])
     hist_fwdDZ.Fill(regularGrid[index][0], regularGrid[index][1], regularGrid[index][2], regularDZ[index])

  hist_fwdDX.Write()
  hist_fwdDY.Write()
  hist_fwdDZ.Write()

histogramFile.Close()

