import numpy as np
from pylab import *
# US Standard atmosphere 1976 model for chosen point 74
P = 39271 #Pa
rho = 0.568607 #kg/m3
T = 240.6 #K

M = 1.4

Vsound = 20.05*np.sqrt(T)
Veq = M*Vsound
Q = 1/2*rho*pow(Veq, 2)
S = 34 # m2
m = 8400 #kg
g0 = 9.81

C = 0.52
f_delta = 0.9
f = 0.608
l_t = 3/2*5.24 #m
X = f*l_t - C*l_t
Y = f_delta*l_t - C*l_t

# Initialization
alpha, Fpx = 0, 0
alphaList = []
alphaList.append(alpha)
# Other parameters
Cz_deltam = 0.59
delta_m0 = 0
k = 0.4
C_mq = -0.36
C_x0 = 0.033
C_zalpha = 2.58
alpha_0 = 0.008


Cz_eq = 1/(Q*S)*(m*g0-Fpx*np.sin(alpha))
Cx_eq = C_x0 + k*pow(Cz_eq, 2)
Cx_deltam = 2*k*Cz_eq*Cz_deltam
delta_meq = delta_m0 -(Cx_eq*np.sin(alpha) + Cz_eq*np.cos(alpha))/(Cx_deltam*np.sin(alpha)+Cz_deltam*np.cos(alpha))
dela_meq = delta_meq*(X/(Y-X))

alpha = alpha_0 + Cz_eq/C_zalpha - Cz_deltam/C_zalpha*dela_meq
alphaList.append(alpha)
Fpx = Q*S*Cx_eq/(np.cos(alpha))
i = 1
epsilon = 0.00000000000000000000001
while abs(alphaList[i] - alphaList[i-1])>epsilon:
    Cz_eq = 1/(Q*S)*(m*g0-Fpx*np.sin(alpha))
    Cx_eq = C_x0 + k*pow(Cz_eq, 2)
    Cx_deltam = 2*k*Cz_eq*Cz_deltam
    delta_meq = delta_m0 -(Cx_eq*np.sin(alpha) + Cz_eq*np.cos(alpha))/(Cx_deltam*np.sin(alpha)+Cz_deltam*np.cos(alpha))
    dela_meq = delta_meq*(X/(Y-X))

    alpha = alpha_0 + Cz_eq/C_zalpha - Cz_deltam/C_zalpha*dela_meq
    alphaList.append(alpha)
    Fpx = Q*S*Cx_eq/(np.cos(alpha))
    i +=1

x = np.linspace(0, len(alphaList), len(alphaList))
print(Fpx)
print(Cz_eq)
print(alpha)

