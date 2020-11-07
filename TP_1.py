# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 09:20:18 2020

@author: 33643
"""
import numpy as np

#Atmosphere characteristics
gamma = 1.4
r = 287.053 
rho = 0.568607 #kg/m3
P = 39271 #Pa
T = 240.6 #K
c = np.sqrt(gamma*r*T) #m/s

#Plane specification
S = 34 #m2
lref = 5.24 #m
lt = (3/2)*lref #m
rg = 2.65 #m
m = 8400 #kg
c = 0.52 #centering

#Flight parameters
M = 1.4
Veq = c*M #m/s
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

Czeq = 0.0449

Q = (1/2)*rho*(Veq**2) 
Cxeq = Cx0 + k*(Czeq**2)

Xv = (2*Q*S*Cxeq)/(m*Veq)
Zv = (2*Q*S*Czeq)/(m*Veq)