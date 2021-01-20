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

As = A[1:6,1:6]
Bs = B[1:6,0:1]
Ds = np.zeros([5,1])

Cs_g = np.array([1,0,0,0,0])
Cs_a = np.array([0,1,0,0,0])
Cs_q = np.array([0,0,1,0,0])
Cs_t = np.array([0,0,0,1,0])
Cs_z = np.array([0,0,0,0,1])

TFs_g = coma.tf(coma.ss(As,Bs,Cs_g,0))
TFs_a = coma.tf(coma.ss(As,Bs,Cs_a,0))
TFs_q = coma.tf(coma.ss(As,Bs,Cs_q,0))
TFs_t = coma.tf(coma.ss(As,Bs,Cs_t,0))
TFs_z = coma.tf(coma.ss(As,Bs,Cs_z,0))

###############################################################q feedback loop
print(coma.pole(TFs_q))
#sisotool(-TFs_q)
Kr = -0.1544

#TF_fb_q = coma.feedback(co.series(coma.tf(Kr,1),TFs_q),sign=-1) #Not working
Aq = As - Kr*Bs*Cs_q
Bq = Kr*Bs

TF_fbq = coma.ss2tf(coma.ss(Aq,Bq,Cs_q,0))
print(coma.damp(TF_fbq))

plt.figure(2)
Y_fbq, T_fbq = coma.step(TF_fbq, T=np.arange(0,10,0.01))
plt.plot(T_fbq, Y_fbq,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("q closed loop step response")
plt.show()

################################################################washout filter


###########################################################gamma feedback loop
#TF_ol_g = coma.series(coma.feedback(coma.tf(Kr,1),TFs_q,sign=-1), TFs_g) #Not working
TF_fbq_g = coma.ss2tf(coma.ss(Aq,Bq,Cs_g,0))
# sisotool(TF_fbq_g)
Kg = 8.11037

#TF_fb_g = coma.feedback(coma.series(coma.tf(Kg,1), TF_ol_g), 1, sign=-1) #Not working

Ag = Aq - Kg*Bq*Cs_g
Bg = Kg*Bq

TF_fbg = coma.ss2tf(coma.ss(Ag,Bg,Cs_g,0))
print(coma.damp(TF_fbg))

plt.figure(3)
Y_fbg, T_fbg = coma.step(TF_fbg, T=np.arange(0,10,0.01))
plt.plot(T_fbg, Y_fbg,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("gamma closed loop step response")
plt.show()

###############################################################z feedback loop
TF_fbg_z = coma.ss2tf(coma.ss(Ag,Bg,Cs_z,0))
sisotool(TF_fbg_z)
Kz = 0.00122

Az = Ag - Kz*Bg*Cs_z
Bz = Kz*Bg

TF_fbz = coma.ss2tf(coma.ss(Az,Bz,Cs_z,0))
print(coma.damp(TF_fbz))

plt.figure(4)
Y_fbz, T_fbz = coma.step(TF_fbz, T=np.arange(0,10,0.01))
plt.plot(T_fbz, Y_fbz,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("altitude closed loop step response")
plt.show()