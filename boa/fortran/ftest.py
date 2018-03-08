"""
NAM: ftest.py
DES: test python wrapped f90 methods
HIS: FB040502
     FB040509 adjusted to new f90 methods and call strategy. added timestamps
USE: execfile('fortran/ftest.py')

"""

#from numarray import array
from Numeric import array
from fortran import f90
import time

# converting all back to Numeric if necessary

Flag   = (BoaB.currData.Data_Flag_p)
Data   = (BoaB.currData.Data_Red_p)
Mean   = (BoaB.currData.Cha_Mean)
Med    = (BoaB.currData.Cha_Med)
Rms    = (BoaB.currData.Cha_Rms)
Mean_s = (BoaB.currData.Cha_Mean_s)
Med_s  = (BoaB.currData.Cha_Med_s)
Rms_s  = (BoaB.currData.Cha_Rms_s)

Nph  = 2
Nch  = 0

Data_Bac_p = Data
Flag_Bac_p = Flag

plotc([Nch+1])
plotp(Nph)
BoaB.Graphic.PlotData.setPlotLimits(['*','*',-2000,10000])
cl()
set(ciData,1)
set(lwBox,5)
show(ciData)
signal()

#-- DEMONSTRATE USE OF get_subscan_index:
BoaB.currData.Subscan_Index = f90.baseline.get_subscan_index(BoaB.currData.Subscan_Index)
f90.data.subscan_index = BoaB.currData.Subscan_Index
print "Subscan_Index=",BoaB.currData.Subscan_Index
print f90.data.subscan_index 

#-- DEMONSTRATE USE OF flag_pc:  (_pc means phase+channel)
print f90.chanana.flag_pc.__doc__
Flag = f90.chanana.flag_pc(Flag,Data,Nph,Nch, 2500,3500)
print "new flag:",Flag[:,Nph,Nch]

#-- DEMONSTRATE USE OF flag:
print f90.chanana.flag.__doc__
Flag[:,Nph,Nch],count=f90.chanana.flag(Flag[:,Nph,Nch],Data[:,Nph,Nch], 2600,3400)
print count," new flag:",Flag[:,Nph,Nch]

BoaB.currData.Data_Flag_p = Flag

#-- DEMONSTRATE USE OF compress:
compress_array,nmax = f90.f1.compress(Data[:,Nph,Nch],Flag[:,Nph,Nch], 0)
print "nmax",nmax
print compress_array[0:nmax]

#-- DEMONSTRATE USE OF compress:
#input_array = Numeric.array(range(5),'f')
#flag_array  = Numeric.array([0,1,0,1,0],'i')
#compress_array,nmax = f90.f1.compress(input_array,flag_array,1)
#compress_array = compress_array[0:nmax]
#print "nmax=",nmax," compressed array=",compress_array

#-- DEMONSTRATE USE OF mean_rms:
t0 = time.clock()
Mean,Med,Rms = f90.chanana.mean_rms(Data,Flag,Mean,Med,Rms)
t1 = time.clock()
print "TIME(mean_rms)  "+str(t1-t0)
print "mean=",Mean[:,Nch]
print "median=",Med[:,Nch]
print "rms=",Rms[:,Nch]
print "mean=",Mean[Nph,:]

BoaB.currData.Cha_Mean = Mean  
BoaB.currData.Cha_Med = Med
BoaB.currData.Cha_Rms = Rms

#-- DEMONSTRATE USE OF addpoly:
#Poly = [20.,-20.,5.,10.,-5.,5.,-5.]
Poly = [10.,0.,10.,20.]
DataAdd = Data
for iph in range(3): DataAdd = f90.baseline.addpoly(DataAdd,Poly,Mean,Rms,iph,Nch) 

BoaB.currData.Data_Red_p = DataAdd
set(ciData,1)
set(symbolPoint,3)
set(chPoint,20)
signal()

#-- DEMONSTRATE USE OF mean_rms_s:
t0 = time.clock()
DataMean = array(DataAdd)
DataMean,Mean_s,Med_s,Rms_s = f90.chanana.mean_rms_s(DataMean,Flag,Mean_s,Med_s,Rms_s)
t1 = time.clock()
print "TIME(mean_rms_s)  "+str(t1-t0)

print "mean=",Mean_s[Nph,Nch,:]
print "median=",Med_s[Nph,Nch,:]
print " rms=",Rms_s[Nph,Nch,:]
BoaB.currData.Data_Red_p = DataMean
set(ciData,6)
signal()

Poly = array(range(2),'f')
DataFit = array(DataAdd)
DataFit,Poly = f90.baseline.polyfit(DataFit,Poly,Flag,Mean,Rms,Nph,Nch) 
BoaB.currData.Data_Red_p = DataFit
set(ciData,2)
set(symbolPoint,1)
signal()

Poly = array(range(3),'f')
DataFit = array(DataAdd)
DataFit,Poly = f90.baseline.polyfit(DataFit,Poly,Flag,Mean,Rms,Nph,Nch) 
BoaB.currData.Data_Red_p = DataFit
set(ciData,3)
signal()

Poly = array(range(4),'f')
DataFit = array(DataAdd)
DataFit,Poly = f90.baseline.polyfit(DataFit,Poly,Flag,Mean,Rms,Nph,Nch) 
BoaB.currData.Data_Red_p = DataFit
set(ciData,4)
signal()

BoaB.currData.Data_Red_p = DataAdd-DataFit
set(ciData,5)
signal()

#print Fit[:,Nph,Nch]
#print Poly

#print f90.data.data_red_p[:,Nph,Nch]
#BoaB.currData.Data_Red_p = f90.data.data_red_p
#print BoaB.currData.Data_Red_p[:,Nph,Nch]

#f90.test2.foo(1,2)
#f90.data.data_red_p[0:2,0,1] = 1000.
#print f90.data.data_red_p[:,1,0]
#print f90.data.data_red_p[0:2,1,0]
#print f90.test1.__doc__ 
#f90.test1.trythis()

# restore initial data:
BoaB.currData.Data_Red_p  = Data_Bac_p
BoaB.currData.Data_Flag_p = Flag_Bac_p

#BoaB.currData.Data_Red_p = DataAdd


#-------------- test reimage:

from Numeric import *
from fortran import f90
im = zeros((11,11),'f')
data = array([10,20,30],'f')
xy  = array([ [-50,0,20],[-30,0,20] ],'f')
dxy = array([ [11,6,11] , [11,6,11] ],'f')
zerooffset = array([-50,-50],'f')
xyscale =  array([0.1,0.1],'f')
im = f90.f1.reimage(im,data,xy,dxy,zerooffset,xyscale)
print im






