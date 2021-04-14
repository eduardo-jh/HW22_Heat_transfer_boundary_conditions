#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW22 - Question 9. Cooling of a chocolate bar
Prof. Waller solution

Created on Mon Apr 12 13:20:01 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

T_air = 20  # ambient temperature in C
L = 0.006  # thickness of solid in m
tmax = 400  # simulation period in sec
nx = 20  # number of space steps
nt = 2000  # number of time steps
b_initial = 40  # initial solid temperature
h = 60  # aerodynamic conductivity in W/m/K
k = 0.19  # thermal conductivity for the solid W/m
Cp = 1170  # specific heat of the solid J /kg/ C
roe = 1200  # density of the solid kg /m3
f = 0.25

alpha = k / (Cp * roe)  # thermal diffusivity
dt = tmax/nt  # time step in sec
dx = L/nx  # space step in m

print('dx =', dx, 'dt =', dt)

T = np.zeros((nx+1,nt+1))  # creating the solution matrix
x = np.linspace(0,L,nx+1)  # the space vector
t = np.linspace(0,tmax,nt+1)  # the time vector

A = (alpha*dt)/(dx**2)
B = 1 - 2*(alpha * dt)/ (dx**2)
print('B=', B)

T[:,0] = b_initial  # the initial boundary condition

conditionleft = 2  # type of the left boundary condition
conditionright = 3  # type pf the right boundary condition

if conditionleft == 1:
    T[0,1:nt+1] = 10
if conditionright == 1:
    T[nx,1:nt+1] = 10

for n in range (nt):

    # Neumann condition dT/dz = 0
    if conditionleft == 2:
        T[0,n] = T[2,n]
    if conditionright == 2:
        T[nx,n] = T[nx-1,n]
        
    if conditionleft == 3:  # I put both forms of the Robbins condition here. It will use the second one.
        # T[0,n] = (dx*h/k*T_air+T[1,n])/(1+dx*h/k)  # first order Taylor series for dT/dz
        T[0,n] = (4*T[1,n]-T[2,n]+2*dx*(h/k*T_air))/(3+2*dx*h/k)  # second order Taylor series for dT/dz
    if conditionright == 3:
        # T[nx,n] = (4*T[nx-1,n]-T[nx-2,n]+2*dx*(h/k*T_air))/(3+2*dx*h/k)  # second order Taylor series for dT/dz
        T[nx,n] = (dx*h/k*T_air+T[nx-1,n])/(1+dx*h/k)  # first order Taylor series for dT/dz
    
    # Dirichlet boundary condition
    for i in range (1,nx):
        T[i,n+1] = A * T[i+1,n] + B * T[i,n] + A * T[i-1,n]+f*dt

plotlayers = [10] 

plt.figure(0)
for layer in plotlayers:
    plt.plot(t,T[layer,:], label = 'T(t) at specific layer ({0})'.format(layer)) 
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.xlabel('Time (sec)')
plt.ylabel('Temperature (°C)')
plt.savefig('q9_chocolate_cooling_time_%dl_%dr.png' % (conditionleft, conditionright), dpi=300, bbox_inches='tight')
plt.show()

plottime = [50, 100, 200]

plt.figure(1)
for t in plottime:
    plt.plot(x,T[:,t], label = 'temps after {0} seconds'.format(t)) 
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('Thickness (m)')
plt.ylabel('Temperature (°C)')
plt.savefig('q9_chocolate_cooling_space_%dl_%dr.png' % (conditionleft, conditionright), dpi=300, bbox_inches='tight')
plt.show()