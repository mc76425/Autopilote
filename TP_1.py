# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 09:20:18 2020

@author: 33643
"""
import numpy as np
import scipy.linalg as sal
import control.matlab as coma
import control as co
import matplotlib.pyplot as plt
from sisopy31 import *

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

sys = coma.ss(A,B,C,D)
coma.damp(sys)
print("\n")

############################################################ Phugoid

Ap = A[0:2,0:2]
Bp = B[0:2,0:1]
Cpv = np.array([1,0])
Cpg = np.array([0,1])
sysp_v = coma.ss(Ap,Bp,Cpv,0)
sysp_g = coma.ss(Ap,Bp,Cpg,0)
print("Eigenvalues of Phugoid mode\n")
[print(round(i,5)) for i in sal.eig(Ap)[0]]

print("\n")
TFp_V_dm = coma.ss2tf(sysp_v)
print("TF V/delta_m = ",TFp_V_dm,"\n")
plt.figure(2)
Yp_v, Tp_v = coma.step(TFp_V_dm, np.arange(0,500,0.1))


print("\n")
TFp_G_dm = coma.ss2tf(sysp_g)
print("TF G/delta_m = ",TFp_G_dm,"\n")
Yp_g, Tp_g = coma.step(TFp_G_dm, np.arange(0,500,0.1))

plt.plot(Tp_v, Yp_v, 'b', Tp_g, Yp_g, 'r')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Phugoid step response")
plt.legend(["V", "gamma"])
plt.show()

print("Step info for V : \n")
a,b,c = step_info(Tp_v, Yp_v)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

print("Step info for gamma : \n")
a,b,c = step_info(Tp_g, Yp_g)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

############################################################ Short period

Asp = A[2:4,2:4]
Bsp = B[2:4,0:1]
Cspa = np.array([1,0])
Cspq = np.array([0,1])
syssp_a = coma.ss(Asp,Bsp,Cspa,0)
syssp_q = coma.ss(Asp,Bsp,Cspq,0)
print("Eigenvalues of Short period mode\n")
[print(round(i,5)) for i in sal.eig(Asp)[0]]

print("\n")
TFsp_a_dm = coma.tf(syssp_a)
print("TF V/delta_m = ",TFsp_a_dm,"\n")
plt.figure(3)
Ysp_a, Tsp_a = coma.step(TFsp_a_dm, np.arange(0,10,0.01))

print("\n")
TFsp_q_dm = coma.tf(syssp_q)
print("TF G/delta_m = ",TFsp_q_dm,"\n")
Ysp_q, Tsp_q = coma.step(TFsp_q_dm, np.arange(0,10,0.01))

plt.plot(Tsp_a, Ysp_a, 'b', Tsp_q, Ysp_q, 'r')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Short-period step response")
plt.legend(["alpha", "q"])
plt.show()

print("Step info for alpha : \n")
a,b,c = step_info(Tsp_a, Ysp_a)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

print("Step info for q : \n")
a,b,c = step_info(Tsp_q, Ysp_q)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

############################################################ New state space vector

As = A[1:6,1:6]
Bs = B[1:6,0:1]
Cs_q = np.array([0,0,1,0,0])
Ds = np.zeros([5,1])
syss = coma.ss(As,Bs,Cs_q,0)
TFs = coma.tf(syss)
print(coma.pole(syss))
#sisotool(-TFs)
syss_CL = co.feedback(0.1544*syss,-1)
print(syss_CL)
TF_CL = coma.ss2tf(syss_CL)
print(TF_CL)
print(coma.damp(syss_CL))

plt.figure(4)
Y_cl, T_cl = coma.step(TF_CL, np.arange(0,5,0.01))

plt.plot(T_cl, Y_cl,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Closed loop step response")
plt.show()