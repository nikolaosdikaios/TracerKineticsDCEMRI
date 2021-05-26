import numpy as np
 
def vec_modToft(etime,tracer_kinetics,aif):
# 2012 Nikolaos Dikaios (n.dikaios@gmail.com)
# PS Toft, Journal Magnetic Resonance Imaging. 1997; 7: 91-101
    t=etime[-1]
    index= int(etime.index(t))
    if tracer_kinetics[2] > 0:
        u=np.asarray(etime[:(index+1)],dtype=float)
        tt=np.asarray(['%.1f' % member for member in [t]*(index+1)],dtype=float)
        Rdelay= np.asarray(np.exp(-tracer_kinetics[1]*np.subtract(tt,u)/tracer_kinetics[2]))[np.newaxis]
        Ca= np.asarray([ aif[i] for i in range(index+1) ])[np.newaxis]
        Ct = float(tracer_kinetics[0]*aif[index]+tracer_kinetics[1]*np.dot(Rdelay,Ca.T))
    else:
        Ct = float(tracer_kinetics[0]*aif[index])


    return Ct;
