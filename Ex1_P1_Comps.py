#!/home/jmartin/anaconda2/envs/py36/bin/python
# coding: utf-8

# In[1]:


# Implement a Lax-Friedrichs solver
# Environment configuration
import numpy as np
import matplotlib.pyplot as plt


# In[2]:

xBins = 500
xLength = 1.0
deltaX = xLength/xBins
x = np.arange(0,xBins)
x = x*deltaX
print(x)
print(deltaX)

# In[3]:


def uNextFromm(u,deltaT):
    # Fromm scheme
    # print("In Fromm")
    dTdX = deltaT/deltaX
    up1 = np.roll(u,-1)
    up2 = np.roll(u,-2)
    um1 = np.roll(u,1)
    um2 = np.roll(u,2)
    fp12 = u - 0.25*(um1 - up1) - 0.5*dTdX*(u - um1)
    fm12 = um1 - 0.25*(um2 - u) - 0.5*dTdX*(um1 - um2)
    temp = u + dTdX*(fm12-fp12)
    return temp

def uNextBW(u,deltaT):
    #  Beam-Warming scheme
    dTdX = deltaT/deltaX
    up1 = np.roll(u,-1)
    up2 = np.roll(u,-2)
    um1 = np.roll(u,1)
    um2 = np.roll(u,2)
    fp12 = 0.5*(3.0*u - um1) - 0.5*dTdX*(u - um1)
    fm12 = 0.5*(3.0*um1-um2) - 0.5*dTdX*(um1-um2)
    temp = u + dTdX*(fm12-fp12)
    return temp

def uNextLaxF(u,deltaT):
    # Lax-Friedrichs Scheme
    up1 = np.roll(u,-1)
    um1 = np.roll(u,1)
    temp = 0.5*(up1+um1) - (0.5*deltaT/deltaX)*(up1-um1)
    return temp

def uNextLaxW(u,deltaT):
    # Lax-Wendroff scheme
    # print('In LaxW')
    dTdX = deltaT/deltaX
    up1 = np.roll(u,-1)
    um1 = np.roll(u,1)
    temp = u - 0.5*dTdX*(up1-um1) + 0.5*dTdX*dTdX*(up1 - 2.0*u + um1)
    return temp

def uNextUpwind(u,deltaT):
    # Upwind differencing scheme
    um1 = np.roll(u,1)
    temp = u + (deltaT/deltaX)*(um1 - u)
    return temp

###################################################
# Implement the initial condition
def IC(vec,bins):
    halfPoint = int(bins/2)
    vec[:halfPoint] = 1.0
    vec[halfPoint:] = 0.0
    return vec

def IC_GAUSS(vec,bins):
    vec = np.exp(-(1./(2*0.05*0.05))*(x - 0.5)*(x-0.5))
    return vec
    


# Implement the initial conditiond
u = np.zeros(x.shape)
u = IC_GAUSS(u,xBins)
uLaxF = uLaxW = uBW = uFromm = uUpwind = u
tFinal = 5.
dt = 0.7*deltaX

plt.figure(figsize=(9,9))
time = 0

while( time < tFinal):
    # uLaxF = uNextLaxF(uLaxF,dt)
    # uLaxW = uNextLaxW(uLaxW,dt)
    # uUpwind = uNextUpwind(uUpwind,dt)
    uBW = uNextBW(uBW,dt)
    # uFromm = uNextFromm(uFromm,dt)
    time = time + dt
    plt.cla()
    plt.ylim(-0.3,1.3)
    plt.title(f"time = {time:.2f}")
    # plt.plot(x,uLaxF,'b.-',label='LaxF')
    # plt.plot(x,uLaxW,'c.-',label='LaxW')
    # plt.plot(x,uUpwind,'k.-',label='Upwind')
    plt.plot(x,uBW,'g.-',label='BW')
    # plt.plot(x,uFromm,'r.-',label='Fromm')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show(block=False)
    plt.pause(0.01)
    
