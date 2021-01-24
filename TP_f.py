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

#######################################################################Phugoid

print('\n-----------------------------Phugoid-------------------------------')
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
plt.legend(["V (m/s)", "gamma (rad)"])
plt.show()

print("Step info for V : \n")
a,b,c = step_info(Tp_v, Yp_v)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

print("Step info for gamma : \n")
a,b,c = step_info(Tp_g, Yp_g)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))
print('\n')

##################################################################Short period

print('\n-----------------------------Short period-------------------------------')
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
print("TF a/delta_m = ",TFsp_a_dm,"\n")
plt.figure(3)
Ysp_a, Tsp_a = coma.step(TFsp_a_dm, np.arange(0,10,0.01))

print("\n")
TFsp_q_dm = coma.tf(syssp_q)
print("TF q/delta_m = ",TFsp_q_dm,"\n")
Ysp_q, Tsp_q = coma.step(TFsp_q_dm, np.arange(0,10,0.01))

plt.plot(Tsp_a, Ysp_a, 'b', Tsp_q, Ysp_q, 'r')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Short-period step response")
plt.legend(["alpha (rad)", "q (rad/s)"])
plt.show()

print("Step info for alpha : \n")
a,b,c = step_info(Tsp_a, Ysp_a)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))

print("Step info for q : \n")
a,b,c = step_info(Tsp_q, Ysp_q)
print("Overshoot (%) : {} Rise time (s) : {} Settling time : {}".format(a,b,c))
print('\n')

##############################################################Simplified model

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
print('\n-----------------------------q feedback-------------------------------')
print(coma.pole(TFs_q))
sisotool(-TFs_q)
Kr = -0.1544

#TF_fb_q = coma.feedback(co.series(coma.tf(Kr,1),TFs_q),sign=-1) #Not working
Aq = As - Kr*Bs*Cs_q
Bq = Kr*Bs

TF_fbq = coma.ss2tf(coma.ss(Aq,Bq,Cs_q,0))
print(coma.damp(TF_fbq))
print(TF_fbq)
print('A : ',Aq,'\nB : ',Bq,'C : ',Cs_q)

plt.figure(4)
Y_fbq, T_fbq = coma.step(TF_fbq, T=np.arange(0,10,0.01))
plt.plot(T_fbq, Y_fbq,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (rad/s)")
plt.title("q closed loop step response")
plt.show()
print('\n')

################################################################washout filter
tau = 0.5
TF_wf = coma.tf([tau,0], [tau, 1])

sys_ol_a = TFs_a
sys_cl_a = coma.series(coma.tf(1/Kr,1),coma.ss2tf(coma.ss(Aq,Bq,Cs_a,0)))
sys_clwf_a = coma.series(coma.tf(1/Kr,1),coma.feedback(coma.tf(Kr,1),coma.series(TFs_q,TF_wf),sign=-1),TFs_a)

plt.figure(5)
y_ol_a, t = coma.step(sys_ol_a,T=np.arange(0,10,0.01))
y_cl_a, t = coma.step(sys_cl_a,T=np.arange(0,10,0.01))
y_clwf_a, t = coma.step(sys_clwf_a,T=np.arange(0,10,0.01))
plt.plot(t, y_ol_a,'b',t,y_cl_a,'r',t,y_clwf_a,'g')
plt.title('alpha step reponses')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (rad/s)")
plt.legend(('Open-loop','Feedback-loop','Feedback-loop with washout filter'))
plt.plot()

###########################################################gamma feedback loop
print('\n---------------------------gamma feedback-----------------------------')
#TF_ol_g = coma.series(coma.feedback(coma.tf(Kr,1),TFs_q,sign=-1), TFs_g) #Not working
TF_fbq_g = coma.ss2tf(coma.ss(Aq,Bq,Cs_g,0))
# sisotool(TF_fbq_g)
Kg = 8.11037

#TF_fb_g = coma.feedback(coma.series(coma.tf(Kg,1), TF_ol_g), 1, sign=-1) #Not working

Ag = Aq - Kg*Bq*Cs_g
Bg = Kg*Bq

TF_fbg = coma.ss2tf(coma.ss(Ag,Bg,Cs_g,0))
print(coma.damp(TF_fbg))
print(TF_fbg)
print('A : ',Ag,'\nB : ',Bg,'C : ',Cs_g)

plt.figure(6)
Y_fbg, T_fbg = coma.step(TF_fbg, T=np.arange(0,10,0.01))
plt.plot(T_fbg, Y_fbg,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (rad)")
plt.title("gamma closed loop step response")
plt.show()
print('\n')

###############################################################z feedback loop
print('\n-----------------------------z feedback-------------------------------')
TF_fbg_z = coma.ss2tf(coma.ss(Ag,Bg,Cs_z,0))
# sisotool(TF_fbg_z)
Kz = 0.00122

Az = Ag - Kz*Bg*Cs_z
Bz = Kz*Bg

TF_fbz = coma.ss2tf(coma.ss(Az,Bz,Cs_z,0))
print(coma.damp(TF_fbz))
print(TF_fbz)
print('A : ',Az,'\nB : ',Bz,'C : ',Cs_z)

plt.figure(7)
Y_fbz, T_fbz = coma.step(TF_fbz, T=np.arange(0,10,0.01))
plt.plot(T_fbz, Y_fbz,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (ft)")
plt.title("altitude closed loop step response")
plt.show()
print('\n')

##########################################################################Satu

print('\n-----------------------------Saturation-------------------------------')
Dnz = 2.8*9.81
a_max = aeq + (aeq-a0)*Dnz
print("alpha max : ",round(np.rad2deg(a_max),2),' deg')

a_i = 0
alpha = a_i
gamma_i = 0
gamma_min = 0
gamma_max = np.pi
i = 0
Tg = np.arange(0,50,0.05)

while (abs(alpha - a_max) > 10**(-4) and i < 10**4):
    gamma_i = (gamma_min + gamma_max)/2
    ug = np.ones_like(Tg)*gamma_i
    y, tg, x = coma.lsim(coma.ss(Ag,Bg,Cs_a,0), U = ug, T=Tg)
    alpha = max(y)
    if(alpha < a_max):
        gamma_min = (gamma_min + gamma_max)/2
    elif(alpha > a_max):
        gamma_max = (gamma_min + gamma_max)/2
    i += 1
gamma_Max = gamma_i
print("gamma max : ",round(np.rad2deg(gamma_Max),2),' deg')
ugg = np.ones_like(Tg)*gamma_Max
ygg, tg, xgg = coma.lsim(coma.ss(Ag,Bg,Cs_a,0), U = ugg, T=Tg)
plt.figure(8)
plt.plot(tg, ygg,'b')
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (rad)")
plt.title("Step response of gamma max")
plt.show()
print('\n')

##########################################################################Simu

t = np.arange(0,300,0.01)
n = int(np.shape(t)[0]/3)
u1 = np.linspace(10000,24000,num=n)
u2 = np.linspace(24000,24000,num=n)
u3 = np.linspace(24000,10000,num=n)
u = np.concatenate((u1,u2,u3))
    
x0 = np.transpose(np.array([0,0,0,0,10000]))

plt.figure(9)
yout, tt, xout = coma.lsim(coma.ss(Az,Bz,np.eye(5),0), U=u, T=t, X0=x0)
plt.plot(tt, yout[:,-1], 'r', tt, u,'b')
plt.xlabel("Time (s)")
plt.ylabel("altitude (ft)")
plt.title("Altitude control simulation A-C-D")
plt.legend(('aircraft altitude','altitude command'))
plt.show()


fig1, axs1 = plt.subplots(2, 2)
axs1[0, 0].plot(tt, yout[:,0])
axs1[0, 0].set_title("gamma")
axs1[1, 0].plot(tt, yout[:,1])
axs1[1, 0].set_title("alpha")
axs1[0, 1].plot(tt, yout[:,2])
axs1[0, 1].set_title("q")
axs1[1, 1].plot(tt, yout[:,3])
axs1[1, 1].set_title("theta")
fig1.tight_layout()

tauf = 10
h0 = 600
tflare = np.arange(0,120,0.01)
uc = np.linspace(600,600,num=int(np.shape(tflare)[0]/4))
uflare = np.concatenate((uc,np.array([ h0*np.exp(-(v-30)/tauf) for v in tflare if v >= 30])))
x0flare = np.transpose(np.array([0,0,0,0,600]))

plt.figure(11)
youtf, ttf, xoutf = coma.lsim(coma.ss(Az,Bz,np.eye(5),0), U=uflare, T=tflare, X0=x0flare)
plt.plot(ttf, youtf[:,-1], 'r', ttf, uflare,'b')
plt.xlabel("Time (s)")
plt.ylabel("altitude (ft)")
plt.title("Flare altitude control simulation")
plt.legend(('aircraft altitude','altitude command'))
plt.show()


fig2, axs2 = plt.subplots(2, 2)
axs2[0, 0].plot(ttf, youtf[:,0])
axs2[0, 0].set_title("gamma")
axs2[1, 0].plot(ttf, youtf[:,1])
axs2[1, 0].set_title("alpha")
axs2[0, 1].plot(ttf, youtf[:,2])
axs2[0, 1].set_title("q")
axs2[1, 1].plot(ttf, youtf[:,3])
axs2[1, 1].set_title("theta")
fig2.tight_layout()