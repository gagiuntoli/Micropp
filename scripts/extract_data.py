from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import vtk
import vtk2numpy as vn
import os
import sys

file_in = "../micropp_1_20.vtu"
print "check for file "+file_in,
try:
  f = open(file_in,'rb')
  print " ok"
except IOError:
  print " error, not found"
  sys.exit()

VTU  = XMLUnstructuredGridReader(FileName=file_in)
KEYs = VTU.GetPointDataInformation().keys(); print KEYs;
KEYs = VTU.GetCellDataInformation().keys(); print KEYs;

#n_pvtu_files = 6
#displ_x = np.zeros( (n_pvtu_files,1) )
#displ_y = np.zeros( (n_pvtu_files,1) )
#force_x = np.zeros( (n_pvtu_files,1) )
#force_y = np.zeros( (n_pvtu_files,1) )
#file_matrix = np.zeros((n_pvtu_files,4))
#
#j = 0
#output_names = ["direct.dat", "homog_us.dat", "homog_tp.dat", "homog_ts.dat"]
#
#for result_path in ["direct", "homog/us", "homog/tp", "homog/ts"]: 
#
#    for force_path in ["force_x", "force_y"]: 
#    
#      for i in range(0,n_pvtu_files):
#      
#        file_in = result_path+"/"+force_path+"/macro_t_"+str(i)+".pvtu"
#        print "check for file "+file_in,
#        try:
#           f = open(file_in,'rb')
#           print " ok"
#        except IOError:
#          print " error, not found"
#          sys.exit()
#      
#        PVTU  = XMLPartitionedUnstructuredGridReader(FileName=file_in)
#        KEYs  = PVTU.GetPointDataInformation().keys()
#        KEYs  = PVTU.GetCellDataInformation().keys()
#      
#        PlotOverLine1 = PlotOverLine(Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source")
#        PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
#        PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
#        PlotOverLine1.Source.Resolution = 1000
#        PlotOverLine1.UpdatePipeline()
#      
#        stress_pvtu = vn.getPtsData(PlotOverLine1, "stress")
#        displ_pvtu = vn.getPtsData(PlotOverLine1, "displ")
#        leng_pvtu = vn.getPtsData(PlotOverLine1, "arc_length")
#    
#        if force_path == "force_x":
#           force_x[i] = np.trapz(stress_pvtu[:,0],leng_pvtu);
#           displ_x[i] = displ_pvtu[0,0]
#        elif force_path == "force_y":
#           force_y[i] = np.trapz(stress_pvtu[:,2],leng_pvtu);
#           displ_y[i] = displ_pvtu[0,1]
#
#
#    for i in range(0,n_pvtu_files):
#      file_matrix[i] = [ displ_x[i], displ_y[i], force_x[i], abs(force_y[i]) ]
#    
#    np.savetxt(output_names[j], file_matrix)
#    j=j+1
