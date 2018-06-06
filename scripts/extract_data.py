from paraview.simple import *
import numpy as np
import glob

n_pvtu_files = tifCounter = len(glob.glob1("./","*.vtu"))
displ_y = np.zeros((n_pvtu_files,1))
force_y = np.zeros((n_pvtu_files,1))
alldata = np.zeros((n_pvtu_files,2))

for i in range(0, n_pvtu_files):
#for i in range(54, 55):

    file_in = "micropp_1_"+str(i)+".vtu"
    try:
      f = open(file_in,'rb')
    except IOError:
      print " error, not found"
      sys.exit()

    print "processing : ",file_in
    
    VTU  = XMLUnstructuredGridReader(FileName=file_in)
    KEYs = VTU.GetPointDataInformation().keys(); #print KEYs;
    KEYs = VTU.GetCellDataInformation().keys();  #print KEYs;
    
    slice1 = Slice(Input=VTU)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.SliceType.Origin = [0.0, 0.999, 0.0]
    slice1.SliceType.Normal = [0.0, 1.0  , 0.0]
    
    IntegrateVariables1 = IntegrateVariables(Input=slice1)
    DataSliceFile = paraview.servermanager.Fetch(IntegrateVariables1)
    #print DataSliceFile 
    
    area       = DataSliceFile.GetCellData().GetArray('Area').GetValue(0)
    sigyy_wall = DataSliceFile.GetCellData().GetArray('stress').GetValue(1)
    displ_y[i] = DataSliceFile.GetPointData().GetArray('displ').GetValue(1)
    force_y[i] = sigyy_wall*area
    #print area, sigyy_wall, displ_y[i], force_y[i]

#----------

for i in range(0,n_pvtu_files):
    alldata[i] = [displ_y[i], force_y[i]]

file_out = "f_vs_d.dat"
np.savetxt(file_out, alldata)
