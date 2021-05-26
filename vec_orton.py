import numpy as np

def vec_orton(etime,tracer_kinetics, aif, dif):
# 2012 Nikolaos Dikaios (n.dikaios@gmail.com)
# MR Orton et al, Phys Med Biol. 2009;54(7):2197-215.
    t=etime[-1]
    index= int(etime.index(t))
    if (tracer_kinetics[1] > 0):
        u=np.asarray(etime[:(index+1)],dtype=float)
        tt=np.asarray(['%.1f' % member for member in [t]*(index+1)],dtype=float)
        Rdelay= np.asarray(np.exp(-tracer_kinetics[0]*np.subtract(tt,u)/tracer_kinetics[1]))[np.newaxis]
        Ca=  np.asarray([ aif[i] for i in range(index+1)])[np.newaxis]
        Cp= np.asarray([ dif[i] for i in range(index+1)])[np.newaxis]
        Cdif= np.asarray(tracer_kinetics[2]*Ca + (1-tracer_kinetics[2])*Cp)
        Ct = float(tracer_kinetics[0]*np.dot(Rdelay,Cdif.T))
    else:
        Ct = 0.0

    return Ct;
