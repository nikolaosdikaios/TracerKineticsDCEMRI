import math;

def Orton08_AIF(ab,mb,ag,mg,t0):
#  This function reproduces the models described in Orton MR, d'Arcy JA, Walker-Samuel S, Hawkes DJ, Atkinson D, Collins DJ, Leach MO.
#  Computationally efficient vascular input function models for quantitative kinetic modelling using DCE-MRI.
#  Phys Med Biol. 2008 Mar 7;53(5):1225-39

# input values as defined in Orton et al 2008, page 1231
#     ab = 4.9; # (2.84 in Orton et al. 2008)
#     mb = 25; #20;
#     ag = 1.36; 
#     mg = 0.171;
#     t  = (0:0.1:60)/60; # must be in minutes
    
    t=t0/60.0; #Needs the 60.0 to understand that this is operation between floats
    Cp = 0.0
    tb = 2.0*math.pi/mb;                             # duration of arterial bolus

    if (t>0) & (t<tb): # indices of arterial bolus
        Cp = ab * (1 - math.cos(mb*t))        # cosine bolus
        Cp = Cp + ab * ag * func(t,mg,mb)
    else:
        pass

    if t >= tb:                             # indices post-bolus
        Cp= ab * ag * func(tb,mg,mb) * math.exp(-mg * (t-tb))
    else:
        pass

    return Cp;

# Function used as decribed in Orton et al. 2008 Eq.(A.2)
def func(th,ah,mh):
    f = (1 - math.exp(-ah*th))/ah - (ah*math.cos(mh*th) + mh*math.sin(mh*th) - ah*math.exp(-ah*th))/(ah**2+mh**2);
    return f;
