#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:23:52 2024

@author: pier
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import BL

# ------------------------------------------------------------------------ INPUT

# pitch law

A0 = 14                                                                         # harmonic pitch law average value [deg]
A1 = 10                                                                         # harmonic pitch law amplitude [deg]

k = 0.077                                                                       # reduced frequency [-]

# inflow conditions

M = 0.1                                                                         # freestream Mach number
T0 = 298.15                                                                     # freestream temperature [K]

# airfoil properties

chord = 0.457                                                                   # airfoil chord [m]
x_AC = 0.25                                                                     # airfoil aerodynamic center [x/c] - depends on Mach

# numerical set-up

N_rev = 10                                                                      # number of pitching cycles
N = 180                                                                         # number of timesteps per cycle

formulation = 'compressible'                                                    # formulation of the attached flow module: incompressible | compressible
fMode = 'raw'                                                                   # handling of the f function: 'fit' use Leishman exponential fitting | 'raw' use data from static polars
vortexModule = 'on'                                                             # activate LEV module: 'on' | 'off'
timeConstantsMod = 'on'                                                         # activate modification of time constants due to LEV: 'on' | 'off'
secondaryVortex ='on'                                                           # activate secondary vortex shedding: 'on' | 'off'

# input data

file_validation = 'reference data/14+10_k0077_M01.txt'                          # reference data
file_constants = 'polarData/S809_constants.txt'                              # BL model constants
file_polar = 'polarData/S809_Re1000k.txt'                                       # airfoil polar data

# ------------------------------------------------------------------------



## initialization

# input validation data

validationData = np.loadtxt(file_validation) 

AOA_val = validationData[:,0] 
CL_val = validationData[:,1] 
CD_val = validationData[:,2] 
CM_val = validationData[:,3] 

CN_val =  CL_val * np.cos(np.deg2rad(AOA_val)) + CD_val * np.sin(np.deg2rad(AOA_val)) 
CC_val =  CL_val * np.sin(np.deg2rad(AOA_val)) - CD_val * np.cos(np.deg2rad(AOA_val)) 

# input calibration data

A = pd.read_table(file_constants, delimiter='\s+', header=None)

calibrationData = np.array(A.iloc[:,1]) 

# input polar data

polarData = np.loadtxt(file_polar) 

AOA = polarData[:,0] 
CL = polarData[:,1] 
CD = polarData[:,2] 
CM = polarData[:,3]  

CN =  CL * np.cos(np.deg2rad(AOA)) + CD * np.sin(np.deg2rad(AOA)) 
CC =  CL * np.sin(np.deg2rad(AOA)) - CD * np.cos(np.deg2rad(AOA)) 

# derived parameters

a = np.sqrt(1.4*287*T0)                                                        # freestream speed of sound [m/s] 
V = M*a                                                                        # freestream velocity [m/s]

omega = (2*k*V)/chord                                                          # pitching pulsation [1/rad]
T = 2*np.pi/omega                                                              # pitching period [s]


## Pitch simulation

# generate pitch law

t = np.linspace(0, N_rev*T, N_rev*N) 
aoa_f= np.deg2rad(A0) + np.deg2rad(A1) * np.sin(omega*t) 
aoa_rate = np.deg2rad(A1) * omega * np.cos(omega*t) 
theta_rate = aoa_rate 
h_rate = aoa_rate * 0.0 

# initialize vectors

cl = np.zeros(len(t)) 
cd = np.zeros(len(t)) 
cn = np.zeros(len(t)) 
ct = np.zeros(len(t)) 
cm = np.zeros(len(t)) 
f_lag = np.zeros(len(t)) 
tv = np.zeros(len(t)) 
comp = np.zeros((len(t),8)) 
bl = np.zeros((len(t),6)) 
state = np.zeros(28) 

dt = t[1]-t[0] 

# run pitching simulation

for i in range(len(t)):

    cn[i], ct[i], cl[i], cd[i], cm[i], f_lag[i], tv[i], comp[i,:], bl[i,:], state = BL.BL(aoa_f[i], aoa_rate[i], theta_rate[i], h_rate[i], V, M, dt, chord, x_AC, calibrationData, polarData, formulation, fMode, timeConstantsMod, vortexModule, secondaryVortex, state) 


## plot data

range_plot = t >= (N_rev - 1) * T

# figure 1 - simulation history

fig1, axs1 = plt.subplots(nrows=8, ncols=1, sharex=True, layout="constrained") 

fig1.set_figheight(16)
fig1.set_figwidth(4)

xlimits = ((N_rev-1), N_rev)

axs1[0].plot(t/T, np.rad2deg(aoa_f),'-k')
axs1[0].set_title('AoA')
axs1[0].set_ylim(A0-A1-5, A0+A1+5)

axs1[1].plot(t/T, aoa_rate,'-k')
axs1[1].set_title('aoarate')
axs1[1].set_ylim(np.min(aoa_rate), np.max(aoa_rate))

axs1[2].plot(t/T, f_lag,'-k')
axs1[2].set_title('f')
axs1[2].set_ylim(0, 1)

axs1[3].plot(t/T, tv/calibrationData[32],'-k')
axs1[3].set_title('tv/Tvl')
axs1[3].set_ylim(0, 3)

axs1[4].plot(t/T, bl[:,0],'-k')
axs1[4].set_title('Tf')
axs1[4].set_ylim(0, 20)

axs1[5].plot(t/T, bl[:,1],'-k')
axs1[5].set_title('Tv')
axs1[5].set_ylim(0, 10)

axs1[6].plot(t/T, cn,'-k')
axs1[6].set_title('CN')
axs1[6].set_ylim([0, 2]) 

axs1[7].plot(t/T, ct,'-k')
axs1[7].set_title('CC')
axs1[7].set_xlabel('t/T [-]')
axs1[7].set_ylim(-0.1, 0.5)
axs1[7].set_xlim(xlimits)

for ax in axs1:
    ax.grid(axis='x', linestyle='--')


# figure 2 - CN, CC

aoa_plot = np.rad2deg(aoa_f[range_plot]) 

cl = cl[range_plot] 
cd = cd[range_plot] 
cm = cm[range_plot] 
cn = cn[range_plot] 
ct = ct[range_plot] 

fig2, axs2 = plt.subplots(nrows=1, ncols=2, layout="constrained") 

fig2.set_figheight(3.5)
fig2.set_figwidth(8)


xlimits = (-5, 30)

axs2[0].plot(AOA, CN, '-o', aoa_plot, cn, '-k', AOA_val, CN_val, '-r')
axs2[0].set_title('CN')

axs2[1].plot(AOA, CC,'-o', aoa_plot, ct, '-k', AOA_val, CC_val, '-r')
axs2[1].set_title('CC')


for ax in axs2:
    ax.set_xlim(xlimits)
    ax.set_xlabel('AoA [deg]')
    ax.grid(linestyle='--')


# figure 3 - CL, CD, CM

fig3, axs3 = plt.subplots(nrows=3, ncols=1, sharex=True, layout="constrained") 

fig3.set_figheight(16)
fig3.set_figwidth(4)


xlimits = (-5, 30)

axs3[0].plot(AOA, CL, '-o', aoa_plot, cl, '-k', AOA_val, CL_val, '-r')
axs3[0].set_title('CL')

axs3[1].plot(AOA, CD, '-o', aoa_plot, cd, '-k', AOA_val, CD_val, '-r')
axs3[1].set_title('CD')

axs3[2].plot(AOA, CM, '-o', aoa_plot, cm, '-k', AOA_val, CM_val, '-r')
axs3[2].set_title('CM')
axs3[2].set_xlabel('AoA [deg]')
axs3[2].set_xlim(xlimits)

axs3[2].legend(['static data', 'BL model','EXP'])

for ax in axs3:
    ax.grid(linestyle='--')

## output to file

P = np.column_stack((aoa_plot, cn, ct, cd, cm))

pd.DataFrame(data=P, columns=['AOA','CL','CD','CM', 'CC']).to_csv('Output/' +'Results.dat', index=False, sep='\t', mode='w', float_format='%-15.4g' )
