import numpy as np
import math
import sys

def Get_lamination_pars(thicknesses, angles):
    n = len(thicknesses)
    
    if n!=len(angles):
        sys.exit("ERROR: thickness and angle arrays must be of equal size!")
    
    h = sum(thicknesses)    
    angles = [math.radians(angle) for angle in angles]
    
    xiA = np.zeros([4])
    xiB = np.zeros([4])
    xiD = np.zeros([4])
    
    z = -0.5*h + 0.5*thicknesses[0]
    
    for ply in range(0,n):
        for i in range(1,5):
            if i==1:
                cs = np.cos(2*angles[ply])
            elif i==2:
                cs = np.cos(4*angles[ply])
            elif i==3:
                cs = np.sin(2*angles[ply])
            elif i==4:
                cs = np.sin(4*angles[ply])
                
            xiA[i-1] = xiA[i-1] + cs*thicknesses[ply]
            xiB[i-1] = xiB[i-1] + cs*thicknesses[ply]*z
            xiD[i-1] = xiD[i-1] + cs*(thicknesses[ply]*z**2 + (thicknesses[ply]**3)/12.0)
    
        if ply<n-1:
            z = z + 0.5*thicknesses[ply] + 0.5*thicknesses[ply+1]
    
    for i in range(0,4):
        xiA[i] = xiA[i] * 1.0 / h
        xiB[i] = xiB[i] * 4.0 / (h**2)
        xiD[i] = xiD[i] * 12.0 / (h**3)
    
    return xiA, xiB, xiD, h



#"Execution"
#thicknesses = [0.125,
#               0.125,
#               0.125,
#               0.125,
#               0.125,
#               0.125,               
#               0.225,               
#               0.125]
#angles = [0.0,
#          90.0,
#          45.0,
#         -45.0,
#         -45.0,
#          45.0,
#          90.0,
#          0.0]

#print Get_lamination_pars(thicknesses, angles)
    
