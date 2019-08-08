import ROOT, math, sys, os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
from mpl_toolkits.mplot3d import Axes3D
####################################################################################################
def read_file(filename):
    with open(filename) as file:
        return file.readlines()

def make_dir(dir_name):
	if(os.path.exists(dir_name) != True):
		os.system('mkdir ' + dir_name)
def pythagoras(x1,x2,y1,y2,z1,z2):
	norm = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
	return norm



print(sys.argv[1])
study_dir = "Track_Plots/"
#p = ROOT.TFile.Open("results_alpha_%d_ConeAng_%d_evr_0_%d_theta_%d.root"%(alpha,phi,evr,theta))
file = ROOT.TFile.Open(sys.argv[1])

font = {'family': 'sans',
        'color':  'black',
        'weight': 'normal',
        'size': 13,
	'rotation' : 270,
        }
label_size = 13
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
#mpl.rcParams['ztick.labelsize'] = label_size
'''
xMax = [2.8, 2.8, 2.8, 2.8, 2.8, 2.8]
xMin = [0.9, 0.2, 1.1, 1.1, 1.5, 2.0]

zMin = [1.8, 5.2, 5.0, 5.2, 6.25, 3.25]
zMax = [3.0, 7.4, 9.0, 6.4, 8.25, 6.75]

yMin = [-1.2, -1.2, -1.5, -1.5, -1.0, -1.4]
yMax = [1.2, 0.5, 1.5, 0.5, 1.0, 1.2]
'''

xMax = [2.8, 2.8, 2.8, 2.8, 2.8, 2.8]
xMin = [2.0, 2.2, 1.9, 2.0, 1.8, 1.8]

zMin = [4.2, 5.2, 5.0, 2.9, 3.0, 7.0]
zMax = [6.1, 7.0, 7.0, 4.1, 5.5, 8.5]

yMin = [-1.2, -1.2, -1.5, -1.5, -1.5, -1.1]
yMax = [0.8, -0.5, 1.0, 0.9, 0.25, 0.85]

Red = '#FF0000'
LightRed = '#F1948A'
Blue = '#0000FF'
LightBlue = '#7FB3D5'
markerSize = 3 	

for plot in range(6):
  
  #if(plot != 5):
  #s  continue
  
  h_3DOrigPlotZero = file.Get( ("tracksOrig%d" % plot) )
  h_3DCorrPlotZero = file.Get( ("tracksCorr%d" % plot) )

  nBinsX = h_3DOrigPlotZero.GetNbinsX()
  nBinsY = h_3DOrigPlotZero.GetNbinsY()
  nBinsZ = h_3DOrigPlotZero.GetNbinsZ()

  print "Plotting : %d" % plot

  xValuesOrig = []
  yValuesOrig = []
  zValuesOrig = []

  xValuesCorr = []
  yValuesCorr = []
  zValuesCorr = []

  xValuesOrigZero = []
  yValuesOrigZero = []
  zValuesOrigZero = []

  xValuesCorrZero = []
  yValuesCorrZero = []
  zValuesCorrZero = []


  offset = 0.0

  for i in range(nBinsX):
    for j in range(nBinsY):
        for k in range(nBinsZ):
            if(h_3DOrigPlotZero.GetBinContent(i,j,k) > 0.0): 
                    '''
		    xValuesOrig.append(h_3DOrigPlotZero.GetYaxis().GetBinCenter(i))
                    yValuesOrig.append(h_3DOrigPlotZero.GetZaxis().GetBinCenter(j))
                    zValuesOrig.append(h_3DOrigPlotZero.GetXaxis().GetBinCenter(k))
		    '''
		    xValuesOrig.append(h_3DOrigPlotZero.GetYaxis().GetBinCenter(j))
                    yValuesOrig.append(h_3DOrigPlotZero.GetZaxis().GetBinCenter(k))
                    zValuesOrig.append(h_3DOrigPlotZero.GetXaxis().GetBinCenter(i))
		    
		    xValuesOrigZero.append(xMax[plot])
		    yValuesOrigZero.append(yMin[plot])
		    zValuesOrigZero.append(zMax[plot])
		    
            
	    if(h_3DCorrPlotZero.GetBinContent(i,j,k) > 0.0):
                    '''
		    xValuesCorr.append(h_3DCorrPlotZero.GetYaxis().GetBinCenter(i))
                    yValuesCorr.append(h_3DCorrPlotZero.GetZaxis().GetBinCenter(j))
                    zValuesCorr.append(h_3DCorrPlotZero.GetXaxis().GetBinCenter(k))
		    '''
		    xValuesCorr.append(h_3DCorrPlotZero.GetYaxis().GetBinCenter(j))
                    yValuesCorr.append(h_3DCorrPlotZero.GetZaxis().GetBinCenter(k))
                    zValuesCorr.append(h_3DCorrPlotZero.GetXaxis().GetBinCenter(i))
		    
		    xValuesCorrZero.append(xMax[plot])
		    yValuesCorrZero.append(yMin[plot])
		    zValuesCorrZero.append(zMax[plot]-offset)


  fig= plt.figure(figsize=(12,8))

  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(zValuesOrig, xValuesOrig, yValuesOrig, s=markerSize,     c=Red, marker='o')
  ax.scatter(zValuesOrig, xValuesOrigZero, yValuesOrig, s=markerSize, c=LightRed, marker='o')
  ax.scatter(zValuesOrigZero, xValuesOrig, yValuesOrig,s=markerSize, c=LightRed, marker='o')
  ax.scatter(zValuesOrig, xValuesOrig, yValuesOrigZero, s=markerSize, c=LightRed, marker='o')

  ax.scatter(zValuesCorr, xValuesCorr, yValuesCorr,  s=markerSize,   c=Blue, marker='o')
  ax.scatter(zValuesCorr, xValuesCorrZero, yValuesCorr, s=markerSize, c=LightBlue, marker='o')
  ax.scatter(zValuesCorrZero, xValuesCorr, yValuesCorr, s=markerSize, c=LightBlue, marker='o')
  ax.scatter(zValuesCorr, xValuesCorr, yValuesCorrZero, s=markerSize, c=LightBlue, marker='o')

  ax.view_init(elev=24,azim=-143)
    
  ax.set_xlabel('Z [m]', fontdict=font, labelpad=20.4)
  ax.set_ylabel('X [m]', fontdict=font, labelpad=20.4)
  ax.set_zlabel('Y [m]', fontdict=font, labelpad=20.4, rotation=90)

  ax.set_zlim3d(yMin[plot], yMax[plot])
  ax.set_xlim3d(zMin[plot], zMax[plot])
  ax.set_ylim3d(xMin[plot], xMax[plot])

  ax.xaxis.pane.fill = False
  ax.yaxis.pane.fill = False
  ax.zaxis.pane.fill = False

  ax.xaxis.pane.set_edgecolor('w')
  ax.yaxis.pane.set_edgecolor('w')
  ax.zaxis.pane.set_edgecolor('w')


    


  plt.savefig(("Figure%d.png" % plot))
#plt.set_size_inches(18.5, 10.5)
#plt.close()

"""
	XYZplots = plt.figure()
	XYZplots.suptitle(Track_info, fontsize=16)
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	XYplot.scatter(x_ord ,y_ord, s = nor_size, c = 'k',linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	XYplot.scatter(x_un_rcl ,y_un_rcl, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	
	#XYplot.scatter(px ,py, s = nor_size, c = 'g',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	ZYplot.scatter(z_ord, y_ord, s = nor_size, c = 'k',linewidth=0)
	ZYplot.plot(z_ord,y_ord)
	ZYplot.scatter(z_un_rcl ,y_un_rcl, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	#ZYplot.scatter(pz ,py, s = nor_size, c = 'g' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	XZplot.scatter(z_ord, x_ord, s = nor_size, c = 'k' ,linewidth=0)
	XZplot.plot(z_ord, x_ord)
	XZplot.scatter(z_un_rcl ,x_un_rcl, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	#XZplot.scatter(pz ,px, s = nor_size, c = 'g',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')
	
	File_name = 'CLOSE_CLUSTERS_Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	XYZplots.set_size_inches(fig_size_x,fig_size_y)
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()


	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	XYplot.scatter(x_ord ,y_ord, s = size_sp_vert, c = col_sp_vert,linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	XYplot.scatter(x_un ,y_un, s = nor_size, c = 'm', linewidth=0)
	
	#XYplot.scatter(px ,py, s = nor_size, c = 'g',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	ZYplot.scatter(z_ord, y_ord, s = size_sp_vert, c = col_sp_vert,linewidth=0)
	ZYplot.plot(z_ord,y_ord)
	ZYplot.scatter(z_un ,y_un, s = nor_size, c = 'm', linewidth=0)
	#ZYplot.scatter(pz ,py, s = nor_size, c = 'g' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	XZplot.scatter(z_ord, x_ord, s = size_sp_vert, c = col_sp_vert,linewidth=0)
	XZplot.plot(z_ord, x_ord)
	XZplot.scatter(z_un ,x_un, s = nor_size, c = 'm',linewidth=0)
	#XZplot.scatter(pz ,px, s = nor_size, c = 'g',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')

	File_name = 'SPVERTEX_TRUEVERTEX_vs_FPVERTEX_Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	XYZplots.set_size_inches(fig_size_x,fig_size_y)
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()

#########################################################
	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	XYplot.scatter(x_ord ,y_ord, s = size_vert, c = col_vert,linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	XYplot.scatter(x_un ,y_un, s = nor_size, c = 'm', linewidth=0)
	
	#XYplot.scatter(px ,py, s = nor_size, c = 'g',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	ZYplot.scatter(z_ord, y_ord, s = size_vert, c = col_vert,linewidth=0)
	ZYplot.plot(z_ord,y_ord)
	ZYplot.scatter(z_un ,y_un, s = nor_size, c = 'm', linewidth=0)
	#ZYplot.scatter(pz ,py, s = nor_size, c = 'g' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	XZplot.scatter(z_ord, x_ord, s = size_vert, c = col_vert,linewidth=0)
	XZplot.plot(z_ord, x_ord)
	XZplot.scatter(z_un ,x_un, s = nor_size, c = 'm',linewidth=0)
	#XZplot.scatter(pz ,px, s = nor_size, c = 'g',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')

	File_name = 'TRUEVERTEX_vs_FPVERTEX_Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()
#########################################################
	
	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	XYplot.scatter(x_un ,y_un, s = nor_size, c = 'm' ,linewidth=0)
	XYplot.scatter(x_ord ,y_ord, s = nor_size, c = 'k',linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	
	
	#XYplot.scatter(px ,py, s = nor_size, c = 'g',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	ZYplot.scatter(z_un ,y_un, s = nor_size, c = 'm' ,linewidth=0)
	ZYplot.scatter(z_ord, y_ord, s = nor_size, c = 'k',linewidth=0)
	ZYplot.plot(z_ord,y_ord)

	#ZYplot.scatter(pz ,py, s = nor_size, c = 'g' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	XZplot.scatter(z_un ,x_un, s = nor_size, c = 'm' ,linewidth=0)
	XZplot.scatter(z_ord, x_ord, s = nor_size, c = 'k' ,linewidth=0)
	XZplot.plot(z_ord, x_ord)

	#XZplot.scatter(pz ,px, s = nor_size, c = 'g',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	File_name = 'ORD_UNORD_Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')

	XYZplots.set_size_inches(fig_size_x,fig_size_y)
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()
	################################################################
	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	XYplot.scatter(x_ord ,y_ord, s = nor_size, c = 'k',linewidth=0)
	XYplot.plot(x_ord ,y_ord)
	XYplot.scatter(x_un ,y_un, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	
	#XYplot.scatter(px ,py, s = nor_size, c = 'g',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)
	ZYplot.scatter(z_ord, y_ord, s = nor_size, c = 'k',linewidth=0)
	ZYplot.plot(z_ord,y_ord)
	ZYplot.scatter(z_un ,y_un, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	#ZYplot.scatter(pz ,py, s = nor_size, c = 'g' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	XZplot.scatter(z_ord, x_ord, s = nor_size, c = 'k' ,linewidth=0)
	XZplot.plot(z_ord, x_ord)
	XZplot.scatter(z_un ,x_un, s = nor_size, c = close,linewidth=0,cmap='coolwarm')
	#XZplot.scatter(pz ,px, s = nor_size, c = 'g',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')
	File_name = 'Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	
	XYZplots.set_size_inches(fig_size_x,fig_size_y)
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()

	
	##########################################################
	XYZplots = plt.figure()
	#############XY Plot###############
	XYplot = XYZplots.add_subplot(221)
	
	XYplot.plot(x_rcl ,y_rcl,color = 'k')
	XYplot.scatter(x_rcl ,y_rcl, s = size, c = col,linewidth=0)
	#XYplot.scatter(x_un ,y_un, s = nor_size, c = 'm',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	ZYplot = XYZplots.add_subplot(222)


	ZYplot.plot(z_rcl,y_rcl, color = 'k')
	ZYplot.scatter(z_rcl, y_rcl, s = size, c = col,linewidth=0)
	#ZYplot.scatter(z_un ,y_un, s = nor_size, c = 'm' ,linewidth=0)
	ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	XZplot = XYZplots.add_subplot(223)
	
	XZplot.plot(z_rcl, x_rcl,color = 'k')
	XZplot.scatter(z_rcl, x_rcl, s = size, c = col ,linewidth=0)
	#XZplot.scatter(z_un ,x_un, s = nor_size, c = 'm',linewidth=0)
	XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")
	
	dqdx_plot = XYZplots.add_subplot(224)
	dqdx_plot.scatter(rr ,dqdx, s = 5, c = 'k',linewidth=10)
	dqdx_plot.set_xlabel('Residual Range [cm]')
	dqdx_plot.set_ylabel(r'dQ/dx [$e^-$/cm]')
	dqdx_plot.set_title(r'dQ/dx vs Resiudal Range')

	File_name = 'RCL_Run_%d_Event_%d_Cluster_%d'%(pt[0],pt[1],pt[2])
	
	XYZplots.set_size_inches(fig_size_x,fig_size_y)
	XYZplots.savefig(study_dir + File_name + ".pdf")
	XYZplots.clf()

	##########################################################
	XYplots = plt.figure()
	XYplot = XYplots.add_subplot(111)
	XYplot.plot(x_rcl ,y_rcl,color = 'k')
	XYplot.scatter(x_rcl ,y_rcl, s = size, c = col,linewidth=0)
	#XYplot.scatter(x_un ,y_un, s = nor_size, c = 'm',linewidth=0)
	#XYplot.scatter(x_un,y_un, c = 'm',linewidth=0)
	#XYplot.set_title("XY plot")
	XYplot.set_xlabel("x [cm]")
	XYplot.set_ylabel("y [cm]")
	
	File_name = 'RCL_Run_%d_Event_%d_Cluster_%d_XY_PLOT'%(pt[0],pt[1],pt[2])

	XYplots.set_size_inches(fig_size_x,fig_size_y)
	XYplots.savefig(study_dir + File_name + ".pdf")
	XYplots.clf()

	###########################
	ZYplots = plt.figure()
	ZYplot = ZYplots.add_subplot(111)

	ZYplot.plot(z_rcl,y_rcl, color = 'k')
	ZYplot.scatter(z_rcl, y_rcl, s = size, c = col,linewidth=0)
	#ZYplot.scatter(z_un ,y_un, s = nor_size, c = 'm' ,linewidth=0)
	#ZYplot.set_title("ZY plot")
	ZYplot.set_xlabel("z [cm]")
	ZYplot.set_ylabel("y [cm]")

	File_name = 'RCL_Run_%d_Event_%d_Cluster_%d_ZY_PLOT'%(pt[0],pt[1],pt[2])

	ZYplots.set_size_inches(fig_size_x,fig_size_y)
	ZYplots.savefig(study_dir + File_name + ".pdf")
	ZYplots.clf()

	###########################
	XZplots = plt.figure()

	XZplot = XZplots.add_subplot(111)
	
	XZplot.plot(z_rcl, x_rcl,color = 'k')
	XZplot.scatter(z_rcl, x_rcl, s = size, c = col ,linewidth=0)
	#XZplot.scatter(z_un ,x_un, s = nor_size, c = 'm',linewidth=0)
	#XZplot.set_title("XZ plot")
	XZplot.set_xlabel("z [cm]")
	XZplot.set_ylabel("x [cm]")

	File_name = 'RCL_Run_%d_Event_%d_Cluster_%d_XZ_PLOT'%(pt[0],pt[1],pt[2])

	XZplots.set_size_inches(fig_size_x,fig_size_y)
	XZplots.savefig(study_dir + File_name + ".pdf")
	XZplots.clf()
"""
	
sys.exit()
