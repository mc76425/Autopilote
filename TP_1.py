# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 09:20:18 2020

@author: 33643
"""
import numpy as np
import scipy.linalg as sal

#Atmosphere characteristics
gamma = 1.4
r = 287.053 
rho = 0.568607 #kg/m3
P = 39271 #Pa
T = 240.6 #K
vs = np.sqrt(gamma*r*T) #m/s
g0 = 9.81 #m/s2

#Plane specification
S = 34 #m2
lref = 5.24 #m
lt = (3/2)*lref #m
rg = 2.65 #m
m = 8400 #kg
c = 0.52 #centering

#Flight parameters
M = 1.4
Veq = vs*M #m/s
h = 24000 #ft
Czdm = 0.59
dm0 = 0
fd = 0.9
k = 0.4
Cmq = -0.36
Cx0 = 0.033
Cza = 2.58
a0 = 0.008
f = 0.608

Czeq = 0.0439
gammaeq = 0
aeq = 0.030
Feq = 62916.1 #N

Q = (1/2)*rho*(Veq**2) 
Cxeq = Cx0 + k*(Czeq**2)
Cxa = 2*k*Cza*Czeq
Cxdm = 2*k*Czdm*Czeq
Iyy = m*(rg**2)

X = -lt*(f-c)
Y = -lt*(fd-c)

Cma = (X/lref)*(Cxa*np.sin(aeq) + Cza*np.cos(aeq))
Cmdm = (Y/lref)*(Cxdm*np.sin(aeq) + Czdm*np.cos(aeq))

Xv = (2*Q*S*Cxeq)/(m*Veq)
Zv = (2*Q*S*Czeq)/(m*Veq)
Xg = (g0*np.cos(gammaeq))/(Veq)
Xa = (Feq*np.sin(aeq)+(Q*S*Cxa))/(m*Veq)
Za = (Feq*np.cos(aeq)+(Q*S*Cza))/(m*Veq)
ma = (Q*S*lref*Cma)/(Iyy)
mq = (Q*S*lref**2*Cmq)/(Veq*Iyy)
Zdm = (Q*S*Czdm)/(m*Veq)
mdm = (Q*S*lref*Cmdm)/(Iyy)

A = np.zeros([6,6])
B = np.zeros([6,1])
C = np.eye(6)
D = np.zeros([6,1])

A[0,0] = -Xv
A[0,1] = -Xg
A[0,2] = -Xa
A[1,0] = Zv
A[1,2] = Za
A[2,0] = -Zv
A[2,2] = -Za
A[2,3] = 1
A[3,2] = ma
A[3,3] = mq
A[4,3] = 1
A[5,1] = Veq

B[1] = Zdm
B[2] = -Zdm
B[3] = mdm

np.set_printoptions(4,suppress=True)
print("A : \n",A)
print("B : \n",B)
print("C : \n",C)
print("D : \n",D)

print(sal.eig(A))