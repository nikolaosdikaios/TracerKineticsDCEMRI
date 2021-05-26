#Test for simulated enhancement
#Developed on Python 2.7.13 (v2.7.13:a06454b1afa1, Dec 17 2016, 20:42:59) [MSC v.1500 32 bit (Intel)] on win32
#2012 Nikolaos Dikaios (n.dikaios@gmail.com)
import math
import scipy.io
import numpy as np
#import mat4py as m4p

#Read 50 dynamic coronal free-breathing T1 weighted images as a matlab (.mat) file
mat_t1w_2D_images = scipy.io.loadmat('t1w_2D_images.mat')
t1w_2D_images = mat_t1w_2D_images['t1w_2D_images']
dimensions= t1w_2D_images.shape
nx=dimensions[0]
ny=dimensions[1]
ng=dimensions[2]

#Assign an index to the different organs and make a mask from the first framr
# Background         1
# Left Heart         2
# Portal vein        3       
# Right Heart        4
# Aorta              5
# Stomach            6
# fat                7
# Liver              8
# Bowel              9
mat_mask_organs = scipy.io.loadmat('mask_organs.mat')
mask_organs = mat_mask_organs['mask_organs']


#----------------------------------------------------
#Acquisition parameters
TR = 2.3046*1e-3; #Repetition time in seconds  
FA = 10*np.pi/180; # flip angle (conversion to radians)
#----------------------------------------------------

#Based on the segmented organ regions - assign tracer kinetics values found in the literature to each organ
vp=[] 
Ktrans=[] #1/seconds
ve=[]
gamma=[]
onset_time=[] #seconds
T10=[] #seconds
for i in range(nx):
    vp.append([])
    Ktrans.append([])
    ve.append([])
    gamma.append([])
    onset_time.append([])
    T10.append([])
    for j in range(ny):
        if ( mask_organs[i][j]==1 or mask_organs[i][j]==10):
            vp[i].append(0)
            Ktrans[i].append(0)
            ve[i].append(0)
            gamma[i].append(0)
            onset_time[i].append(0)
            T10[i].append(0.382)
        elif ( mask_organs[i][j]==2):
            vp[i].append(0.6)
            Ktrans[i].append(0)
            ve[i].append(0)
            gamma[i].append(0)
            onset_time[i].append(10)
            T10[i].append(1.471)
        elif ( mask_organs[i][j]==3):
            vp[i].append(0)
            Ktrans[i].append(1.15/60)
            ve[i].append(0.41)
            gamma[i].append(0.74)
            onset_time[i].append(27)
            T10[i].append(1.932)
        elif ( mask_organs[i][j]==4):
            vp[i].append(1)
            Ktrans[i].append(0)
            ve[i].append(0)
            gamma[i].append(0)
            onset_time[i].append(5)
            T10[i].append(1.471)
        elif ( mask_organs[i][j]==5):
            vp[i].append(1)
            Ktrans[i].append(0)
            ve[i].append(0)
            gamma[i].append(0)
            onset_time[i].append(15)
            T10[i].append(1.932)
        elif ( mask_organs[i][j]==6):
            vp[i].append(0.001)
            Ktrans[i].append(0.7/60)
            ve[i].append(0.5)
            gamma[i].append(0)
            onset_time[i].append(27)
            T10[i].append(1.328)
        elif ( mask_organs[i][j]==7):
            vp[i].append(0)
            Ktrans[i].append(0)
            ve[i].append(1)
            gamma[i].append(0)
            onset_time[i].append(1)
            T10[i].append(0.382)
        elif ( mask_organs[i][j]==8):
            vp[i].append(0)
            Ktrans[i].append(0.23*6/60)
            ve[i].append(0.23)
            gamma[i].append(0.25)
            onset_time[i].append(25)
            T10[i].append(0.809)
        elif ( mask_organs[i][j]==9):
            vp[i].append(0.01)
            Ktrans[i].append(0.55/60)
            ve[i].append(0.4)
            gamma[i].append(0)
            onset_time[i].append(40)
            T10[i].append(1.328)
        else:
            vp[i].append(0)
            Ktrans[i].append(0)
            ve[i].append(0)
            gamma[i].append(0)
            onset_time[i].append(0)
            T10[i].append(0)

#Save tracer kinetic maps in a matlab (.mat) format
scipy.io.savemat('T10.mat',mdict={'T10': T10})
scipy.io.savemat('vp.mat',mdict={'vp': vp})
scipy.io.savemat('Ktrans.mat',mdict={'Ktrans': Ktrans})
scipy.io.savemat('ve.mat',mdict={'ve': ve})
scipy.io.savemat('gamma.mat',mdict={'gamma': gamma})
scipy.io.savemat('onset_time.mat',mdict={'onset_time': onset_time})
#--------------------------------------------------------------

# Set simulated DCE acquisition parameters

# The time of acquisition of each "gate" is in the list named time.

##### Example 1: The  time resolution of the dynamic acquisition is constant
nt = 150 # Number of time frames
time_resolution = 1 # 1 image acquired every 1 seconds
acq_time = nt*time_resolution # Acquisition length (sec)
time= range(1,acq_time+1,time_resolution) #Define the timing of each acquired image

#### Example 2: The  time resolution of the dynamic acquisition is irregular
##nt=10
##acq_time=63.6
##time= [1.2, 1.8, 2.3, 3.1, 5.4 ,6.1, 8.9, 10.3, 42.7, 63.6]


#--- Arterial and Dual input function based on a mathematical model from Orton et.al. (2008)---
from Orton08_AIF import Orton08_AIF
aif= [0]*(int(acq_time));
ab = 4.9 
mb = 22.8  
ag = 1.36
mg = 0.171
i=0
for t in time:
    aif[i]= Orton08_AIF(ab,mb,ag,mg,t)
    i=i+1

dif= [0]*(int(acq_time));
ab = 1.69; 
mb = 11.8; 
ag = 2.33; 
mg = 0.145;
tp= 0.09*60;
i=0
for t in time:
    dif[i]= Orton08_AIF(ab,mb,ag,mg,t-tp)
    i=i+1
    
#--------------------------------------------------------------

##from SimulateSignalEnhancement import DCEsignal
##aa=DCEsignal(0.01, 0.55/60, 0.6, 0.23, 23, 1000, aif, dif, time, TR, 4.51, FA, 1,'Toft');
##print(aa)
#--- Simulate Signal Enhancement images ---
from SimulateSignalEnhancement import DCEsignal
r1 = 4.51; # T10 relaxivity (s^-1 nM^-1)
new = [[[0 for k in range(nt)] for j in range(ny)] for i in xrange(nx)]
S0= [[0 for j in range(ny)] for i in range(nx)]
for i in range(nx):
    for j in range(ny):
        S0[i][j] = t1w_2D_images[i][j][0]* (1-np.exp(-TR/T10[i][j])*np.cos(FA)) / ((1-np.exp(-TR/T10[i][j]))*np.sin(FA))

        if (mask_organs[i][j] == 8 or mask_organs[i][j] == 3):
            new[i][j][:]=DCEsignal(vp[i][j], Ktrans[i][j], ve[i][j], gamma[i][j], onset_time[i][j], S0[i][j], aif, dif, time, TR, r1, FA, T10[i][j],'Orton');
        else:
            new[i][j][:]=DCEsignal(vp[i][j], Ktrans[i][j], ve[i][j], gamma[i][j], onset_time[i][j], S0[i][j], aif, dif, time, TR, r1, FA, T10[i][j],'Toft');

simulated_dce_images= new
scipy.io.savemat('S0.mat',mdict={'S0': S0})
scipy.io.savemat('simulated_dce_images.mat',mdict={'simulated_dce_images': simulated_dce_images})

