#Simulate signal enhancement based on the modified Toft model or the dual Input Orton
import numpy as np
from vec_orton import vec_orton
from vec_modToft import vec_modToft
from operator import sub
from itertools import imap

def DCEsignal(vp, Ktrans, ve, gamma, onset_time, S0, aif, dif, time, TR, rc1, fa, T10, tk_model):
# generate DCE signal from extended Toft or Orton model parameters, transfer rate
# Ktrans are in (1/secs), onset time, T10 and TR are in secs and fa is in radians
# 2012 Nikolaos Dikaios (n.dikaios@gmail.com)
    if T10 != 0:    
        R10=float(1.0/T10)
        sina= np.sin(fa); cosa= np.cos(fa);
        constant= S0*(1-np.exp(-TR*R10))*sina/(1-cosa*np.exp(-TR*R10));
        nt=len(time)
        St= [0]*nt
        Rt= [0]*nt
        enh_time=[]
        if tk_model == 'Toft':
            tk= [vp, Ktrans, ve, onset_time]; i=0;
            for t in time:
                j=0
                if (t >= onset_time and onset_time != 0):
                    enh_time.append(t-onset_time+1)
                    Rt = R10+rc1*vec_modToft(enh_time,tk,aif) #omitted T10 % St(t) = k(4)*(sin(phi)/(1-cos(phi)))*TR*Rt(t); %
                    St[i] = S0*(1-np.exp(-TR*Rt))*sina/(1-cosa*np.exp(-TR*Rt))
                    j=j+1
                else:
                    St[i]= constant
                    
                i=i+1

        elif tk_model == 'Orton':
            tk= [Ktrans, ve, gamma, onset_time]; i=0; 
            for t in time:
                j=0
                if (t >= onset_time and onset_time != 0):
                    enh_time.append(t-onset_time+1)
                    Rt = R10+rc1*vec_orton(enh_time,tk,aif,dif)
                    St[i] = S0*(1-np.exp(-TR*Rt))*sina/(1-cosa*np.exp(-TR*Rt))
                    j=j+1
                else:
                    St[i]= constant;
                    
                i=i+1;

    else:
        St= [S0]*nt
        

    return St;

