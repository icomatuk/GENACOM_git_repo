import numpy as np
import math
import sys

from Get_lamination_pars import Get_lamination_pars


#def Get_ABD(thicknesses, angles):


thicknesses = [0.125,
               0.125,
               0.125,
               0.125,
               0.125,
               0.125,               
               0.125,               
               0.125]
angles = [0.0,
          90.0,
          45.0,
         -45.0,
         -45.0,
          45.0,
          90.0,
          0.0]

xiA, xiB, xiD, h = Get_lamination_pars(thicknesses, angles)

'--------------------------------------------------------------'

E11 = 150000
E22 = 20000
G12 = 5000
v12 = 0.3

Q11 = (E11**2)/(E11-E22*v12**2)
Q22 = (E11*E22)/(E11-E22*v12**2)
Q12 = v12*Q22
Q66 = G12

U1 = (1/8.0)*(3*Q11+3*Q22+2*Q12+4*Q66)
U2 = (1/2.0)*(Q11-Q22)
U3 = (1/8.0)*(Q11+Q22-2*Q12-4*Q66)
U4 = (1/8.0)*(Q11+Q22+6*Q12-4*Q66)
U5 = (1/8.0)*(Q11+Q22-2*Q12+4*Q66)

U = np.array([[U1],
              [U2],
              [U3],
              [U4],
              [U5]])

xiA_mat = h*np.array([[1, xiA[0], xiA[1], 0, 0],
                     [1, -xiA[0], xiA[1], 0, 0],
                     [0, 0, -xiA[1], 1, 0],
                     [0, 0, -xiA[1], 0, 1],
                     [0, 0.5*xiA[2], xiA[3], 0, 0],
                     [0, 0.5*xiA[2], -xiA[3], 0, 0]])

xiB_mat = ((h**2)/4.0)*np.array([[0, xiB[0], xiB[1], 0, 0],
                                 [0, -xiB[0], xiB[1], 0, 0],
                                 [0, 0, -xiB[1], 0, 0],
                                 [0, 0, -xiB[1], 0, 0],
                                 [0, 0.5*xiB[2], xiB[3], 0, 0],
                                 [0, 0.5*xiB[2], -xiB[3], 0, 0]])
    
xiD_mat = ((h**3)/12.0)*np.array([[1, xiD[0], xiD[1], 0, 0],
                                  [1, -xiD[0], xiD[1], 0, 0],
                                  [0, 0, -xiD[1], 1, 0],
                                  [0, 0, -xiD[1], 0, 1],
                                  [0, 0.5*xiD[2], xiD[3], 0, 0],
                                  [0, 0.5*xiD[2], -xiD[3], 0, 0]])

A = np.matmul(xiA_mat,U)
B = np.matmul(xiB_mat,U)
D = np.matmul(xiD_mat,U)

'--------------------------------------------------------------'

print A
print B
print D















