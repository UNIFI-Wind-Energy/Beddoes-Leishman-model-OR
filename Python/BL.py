#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 10:23:30 2024

@author: pier
"""

import numpy as np

def f_BL(alpha, alpha1, S1, S2):
    
    '''
    LEISHMAN SEPARATION POINT  exponential fitting of raw separation point
    data
    '''

    if np.abs(alpha) <= alpha1:
    
        f = 1 - 0.3 *  np.exp((np.abs(alpha)-alpha1)/S1)
    
    else:
    
        f = 0.04 + 0.66 *  np.exp((alpha1-np.abs(alpha))/S2)

    return f



def evaluatePolar(alpha_0, m_CN, CM0, polarData):

    '''
    BL - EVALUATE POLAR Compute static characteristics of the boundary layer
    - separation point f and center of pressure x_cp profiles - from input
    polar data
    '''
    
    # polar data pre-processing
    
    AOA = polarData[:,0]
    CL = polarData[:,1]
    CD = polarData[:,2]
    CM = polarData[:,3]
    
    AOA = np.deg2rad(AOA) 
    
    CN = CL * np.cos(AOA) + CD * np.sin(AOA) 
    
    # compute static separation point f
    
    f = np.zeros(len(AOA))
    
    for i in range(len(AOA)):
    
        f[i] = (2 * np.sqrt( np.abs(CN[i] / ( m_CN*(AOA[i]-alpha_0)+1e-6)) ) - 1)**2
            
        # apply boundaries saturation - f c [0 1]
    
        if f[i] > 1:
    
            f[i] = 1 - 1e-12
        
        elif f[i] < 0:
        
            f[i] = 1e-12
        
    
    # compute static center of pressure xCP
    
    x_CP = (CM-CM0) / CN 
    
    
    return AOA, f, x_CP



def BL_attachedFlow_incompressible(alpha, theta_dot, V, dt, chord, bsc, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, state):

    '''
    ATTACHED FLOW MODULE Computes potential unsteady loads in the
    Beddoes-Leishman dynamic stall model
    '''

    # initialisation
    
    # update state
    
    alpha_prev = state[0]
    theta_dot_prev = state[1]
    V_prev = state[2]
    X1_prev = state[3]
    Y1_prev = state[4]
    w_34_prev = state[5]

    DP_prev = state[14]
    CN_prev = state[15]
    
    # derived quantities
    
    ds = (V+V_prev) * dt / chord                                                   # non-dimensional timestep [-]
    
    Delta_alpha = alpha - alpha_prev
    Delta_theta_dot = theta_dot - theta_dot_prev
    Delta_V = V - V_prev
    
    
    # -------------------------------------------------------------------------- CN
    
    
    # ---------------------circulatory part
    
    
    # relative downwash at 3/4 chord 
    
    w_34 = V * alpha - (bsc - 0.75) * chord * theta_dot
    
    # deficiency functions - shed vorticity effect
    
    Delta_w_34 = w_34 - w_34_prev
    
    X1 = X1_prev * np.exp(-b1*ds) + A1 * Delta_w_34 * np.exp(-0.5*b1*ds)
    Y1 = Y1_prev * np.exp(-b2*ds) + A2 * Delta_w_34 * np.exp(-0.5*b2*ds)
    
    w_34_E = w_34 - X1 - Y1
    
    # effective angle of attack at 3/4 chord -> circulatory loads computation
    
    alphaE = w_34_E / V
    
    CN_C = m_CN * (alphaE - alpha0)
    

    # --------------------- impulsive part
    
    # relative downwash acceleration at 1/2 chord
    
    w_12_dot = V * Delta_alpha/dt + alpha * Delta_V/dt - (bsc - 0.5) * chord * Delta_theta_dot/dt
    
    CN_I = m_CN/4 * chord * w_12_dot/V**2
    
    
    # ------------------------- total load
    
    CN = CN_C + CN_I
    
    
    # -------------------------------------------------------------------------- LE pressure lag
    
    DP = DP_prev * np.exp(-ds/TP) + (CN-CN_prev) * np.exp(-0.5*ds/TP)
    
    CN_lag = CN - DP
    
    alpha_lag = alpha0 + CN_lag/m_CN
    
    
    # -------------------------------------------------------------------------- CC
    
    CC = m_CN * (alphaE - alpha0) * np.tan(alphaE)
    
    
    # -------------------------------------------------------------------------- CM
    
    
    # ------------------------- circulatory part
    
    # angle of attack change
    
    CM_alpha_c = (0.25 - x_AC) * CN_C
    
    # pitch rate (virtual camber) component
    
    CM_q_c = - m_CN/16 * chord/V * theta_dot
    
    
    # ------------------------- impulsive part
    
    CM_I = - 0.25 * CN_I - m_CN/128 * chord**2 / V**2 * Delta_theta_dot/dt
    
    
    
    # update previous timestep
    
    # [  0        1     2  3  4  5  6  7  8  10 11 12 13  14 15]
    # [alpha theta_dot  V  X1 Y1 w_34                     DP CN] 
    
    state[0] = alpha
    state[1] = theta_dot
    state[2] = V
    state[3] = X1
    state[4] = Y1
    state[5] = w_34
    state[14] = DP
    state[15] = CN
    
    return CN_C, CN_I, CN_lag, alpha_lag, CC, CM_alpha_c, CM_q_c, CM_I, ds, state



def BL_attachedFlow(alpha, dthetadt, V, M, dt, chord, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, A3, b3, A4, b4, A5, b5, state):

    '''
    # ATTACHED FLOW MODULE Computes potential unsteady loads in the
    # Beddoes-Leishman dynamic stall model
    '''
    
    # initialisation
    
    # update state
    
    alpha_prev = state[0] 
    q_prev = state[1] 
    Delta_alpha_prev = state[2]
    Delta_q_prev = state[3]
    X1_prev = state[4]
    Y1_prev = state[5]
    D_prev = state[6]
    X2_prev = state[7]
    Y2_prev = state[8]
    Dq_prev = state[9]
    X3_prev = state[10]
    X4_prev = state[11]
    Y4_prev = state[12]
    DMq_prev = state[13]
    DP_prev = state[14]
    CN_prev = state[15]
    
    # derived quantities
    
    ds = 2*V*dt/chord                                                           # non-dimensional timestep [-]
    
    q = dthetadt*chord/V                                                        # non-dimensional pitch rate [-]
    
    a = V/M                                                                     # speed of sound [m/s]
    beta = np.sqrt(1-M**2)                                                      # Glauert compressibility factor [-]
    
    Delta_alpha = alpha - alpha_prev 
    Delta_q = q - q_prev 
    
    
    ## -------------------------------------------------------------------------- CN
    
    # ---------------------circulatory part
    
    # angle of attack change
    
    X1 = X1_prev * np.exp(-b1*beta**2*ds) + A1 * Delta_alpha * np.exp(-0.5*b1*beta**2*ds) 
    
    Y1 = Y1_prev * np.exp(-b2*beta**2*ds) + A2 * Delta_alpha * np.exp(-0.5*b2*beta**2*ds) 
    
    CN_alpha_c = m_CN * (alpha - X1 - Y1 - alpha0) 
    
    # pitch rate change
    
    X2 = X2_prev * np.exp(-b1*beta**2*ds) + A1 * Delta_q * np.exp(-0.5*b1*beta**2*ds) 
    
    Y2 = Y2_prev * np.exp(-b2*beta**2*ds) + A2 * Delta_q * np.exp(-0.5*b2*beta**2*ds) 
    
    qE = q - X2 - Y2 
    
    CN_q_c = 0.5 * m_CN * qE 
    
    # normal load
    
    CN_C = CN_alpha_c + CN_q_c 
    
    # effective angle of attack at 3/4 chord -> computation of CC
    
    alphaE = (alpha - X1 - Y1) + qE/2 
    
    #--------------------- impulsive part
    
    TI = chord/a 
    
    # angle of attack change
    
    K_alpha = 1 / ( 1-M + 0.5 * m_CN * beta**2 * M**2 * (A1*b1+A2*b2) ) 
    
    T_alpha = 0.75 * K_alpha * TI 
    
    D = D_prev * np.exp(-dt/T_alpha) + (Delta_alpha-Delta_alpha_prev)/dt * np.exp(-0.5*dt/T_alpha) 
    
    CN_alpha_i = 4*T_alpha/M * (Delta_alpha/dt - D) 
    
    # pitch rate change
    
    K_q = 1 / ( 1-M + m_CN*beta**2 * M**2 * (A1*b1+A2*b2) ) 
    
    T_q = 0.75 * K_q * TI 
    
    Dq = Dq_prev * np.exp(-dt/T_q) + (Delta_q-Delta_q_prev)/dt * np.exp(-0.5*dt/T_q) 
    
    CN_q_i = T_q/M * (Delta_q/dt - Dq) 
    
    # impulsive load
    
    CN_I = CN_alpha_i + CN_q_i 
    
    # ------------------------- total load
    
    CN = CN_C + CN_I 
    
    
    ## -------------------------------------------------------------------------- LE pressure lag
    
    DP = DP_prev * np.exp(-ds/TP) + (CN-CN_prev) * np.exp(-0.5*ds/TP) 
    
    CN_lag = CN - DP 
    
    alpha_lag = alpha0 + CN_lag/m_CN 
    
    ## -------------------------------------------------------------------------- CC
    
    CC = m_CN * (alphaE - alpha0) * np.tan(alphaE) 
    
    
    ## -------------------------------------------------------------------------- CM
    
    # ------------------------- circulatory part
    
    # angle of attack change
    
    CM_alpha_c = (0.25-x_AC) * CN_alpha_c  
    
    # pitch rate change
    
    X3 = X3_prev * np.exp(-b5*beta**2*ds) + A5 * Delta_q * np.exp(-0.5*b5*beta**2*ds)  
    
    CM_q_c= -m_CN/16 * (q - X3)  
    
    # ------------------------- impulsive part
    
    # angle of attack change
    
    KalphaM = (A3*b4+A4*b3) / (b3*b4*(1-M)) 
    
    X4 = X4_prev * np.exp(-dt/(b3*KalphaM*TI)) + (Delta_alpha-Delta_alpha_prev)/dt * np.exp(-0.5*dt/(b3*KalphaM*TI)) 
    Y4 = Y4_prev * np.exp(-dt/(b4*KalphaM*TI)) + (Delta_alpha-Delta_alpha_prev)/dt * np.exp(-0.5*dt/(b4*KalphaM*TI)) 
    
    CM_alpha_i = - 1/M * ( A3 * (b3*KalphaM*TI) * (Delta_alpha/dt-X4) + A4 * (b4*KalphaM*TI) * (Delta_alpha/dt-Y4) )  
    
    # pitch rate change
    
    KqM = 7 / ( 15*(1-M) + 1.5*m_CN*beta**2 * M**2 *A5*b5 )  
    
    DMq = DMq_prev * np.exp(-dt/(KqM*TI)) + (Delta_q-Delta_q_prev)/dt * np.exp(-0.5*dt/(KqM*TI)) 
    
    CM_q_i = - 7*KqM*TI/(12*M) * (Delta_q/dt - DMq) 
    
    # impulsive load
    
    CM_I = CM_alpha_i + CM_q_i 
    
    
    ## update previous timestep
    
    # [              1   2      3          4   5  6  7 8  9  10 11 12 13  14 15 16]
    
    state[0:16] = [alpha, q, Delta_alpha, Delta_q, X1, Y1, D, X2, Y2, Dq, X3, X4, Y4, DMq, DP, CN] 

    
    return CN_alpha_c, CN_q_c, CN_I, CN_lag, alpha_lag, CC, CM_alpha_c, CM_q_c, CM_I, state



def  BL_computeTimeConstants(alpha, dalphadt, tv, f_lag, CN_lag, alpha10, delta_alpha1, CN1, Tf0, Tvl, dfdt):

    '''
    TIME CONSTANTS Modify boundary layer dynamic response based on the system and LEV status 
    
    '''
    
    # compute dynamic separation point speed Tf
    
    # separation speed is increased if beyond a certain amount of separated flow (Leishman, 1986)
    
    if f_lag >= 0.7:
        Tf = Tf0 
    else:
        Tf = 0.5*Tf0 
    
    # vortex shedding - presence of the LEV modifies boundary layer behaviour
    
    if (abs(CN_lag) >= CN1): 
    
        # case #1: separating flow during vortex shedding - motion of the LEV along 
        # the chord accelerates separation in both pitch-up and pitch down
    
        if (dfdt < 0): 
    
            if tv > 0 and tv < Tvl:
                Tf = 0.5 * Tf0 
            elif tv >= Tvl and tv < 2*Tvl:
                Tf = Tf0 
            else:
                Tf = Tf0 
    
        # case #2: re-attaching flow during vortex shedding - presence of the LEV 
        # in the near-wake slows-down flow re-attachment 
    
        else: 
    
            if tv > 0 and tv < 2*Tvl:
                Tf = 4*Tf0 

    
    # case #3: re-attaching flow after vortex shedding - recovery capacity of
    # the boundary layer is restored
    
    elif np.abs(CN_lag) < CN1 and dfdt >= 0:
    
        Tf = 2*Tf0 
    

    ## compute separation angle alpha1
    
    alpha1 = alpha10 
    
    if np.abs(CN_lag) >= CN1 and alpha*dalphadt < 0:                               # vortex shedding
    
        alpha1 = alpha10 - (1-f_lag) **0.25 * delta_alpha1 
       

    return Tf, alpha1
    


def BL_separatedFlow(f_lag, x_CP, CN_C, CC_pot, eta):

    '''
    # BEDDOES-LEISHMAN - SEPARATED FLOW MODULE Computation of blade loads
    # during TE separation
    '''
    
    ## ------------------------------------------------------------------------ unsteady loads
    
    # normal load
    
    CN_f = CN_C * ( (1+np.sqrt(f_lag))/2 )**2 
    
    # tangential load
    
    CC_f = eta * CC_pot * np.sqrt(f_lag) 
    
    # pitching moment
        
    CM_f = CN_f * x_CP 
       
    return CN_f, CC_f, CM_f



def BL_unsteadyBoundaryLayer(alpha, alpha_lag, ds, alpha1, S1, S2, Tf, K0, K1, K2, m, Tf0, F1, AOA, f_raw, x_CP_raw, fMode, state):

    '''
    BOUNDARY LAYER DYNAMIC Compute unsteady separation point position f_lag
    and center of pressure x_CP
    '''
    ## initialization
    
    f_lag_prev = state[16] 
    Df_prev = state[17] 
    fM_prev = state[18] 
    DfM_prev = state[19] 
    
    alpha_lag_prev = state[20] 
    Dalpha_prev = state[21] 
    
    ## ------------------------------------------------------------------------ unsteady separation point (Kirchhoff theory)
    
    if fMode == 'raw':
    
        f = np.interp(alpha, AOA, f_raw) 
    
        f_lag = np.interp(alpha_lag, AOA, f_raw) 
    
    elif fMode == 'fit':
    
        f = f_BL(alpha, alpha1, S1, S2) 
    
        f_lag = f_BL(alpha_lag, alpha1, S1, S2) 
    
    
    Df = Df_prev * np.exp(-ds/Tf) + (f_lag-f_lag_prev) * np.exp(-0.5*ds/Tf) 
    
    f_lag2 = f_lag - Df 
    
    
    ## ------------------------------------------------------------------------ unsteady center of pressure (Kirchhoff theory)
    
    fM = f 
    
    DfM = DfM_prev * np.exp(-ds/(F1*Tf0)) + (fM-fM_prev) * np.exp(-0.5*ds/(F1*Tf0)) 
    
    fM_lag = fM - DfM 
    
    # if raw x_CP data is used, delayed angle of attack instead of separation
    # point f'' is used for data query - time constant Tf is the same
    
    Dalpha = Dalpha_prev * np.exp(-ds/Tf) + (alpha_lag-alpha_lag_prev) * np.exp(-0.5*ds/Tf) 
    
    alpha_lag2 = alpha_lag - Dalpha 
    
    
    if fMode == 'raw':
    
        x_CP = np.interp(alpha_lag2, AOA, x_CP_raw) 
    
    elif fMode == 'fit':
    
        if f_lag2 >= fM_lag:
    
            x_CP = K0 + K1 * (1-f_lag2) + K2 * np.sin(np.pi * f_lag2**m) 
    
        else:
    
            x_CP = K0 + K1 * (1-fM_lag) + K2 * np.sin(np.pi * fM_lag**m) 
            
    
    ## update previous timestep
    
    state[16:22] = [f_lag, Df, fM, DfM, alpha_lag, Dalpha] 
    
    
    return f_lag2, x_CP, f, fM, state



def BL_vortexShedding(alpha, dalphadt, CN_C, CN_f, dt, ds, tv, Tv0, Tvl, timeConstantsMod, state):

    '''
    BEDDOES-LEISHMAN VORTEX SHEDDING MODULE Compute unsteady loads due to
    primary (LEV) and secondary vortex shedding
    '''    
    
    Cv_prev = state[22] 
    CN_v_prev = state[23] 
    
    ## compute load response time to LEV Tv
    
    Tv = Tv0 
    
    if timeConstantsMod == 'on':
    
        if alpha*dalphadt >= 0: # pitch-up
    
            if tv > 0 and tv <= Tvl:
                Tv = Tv0 
            elif tv > Tvl and tv <= 2*Tvl:
                Tv = 0.5*Tv0 
            else:
                Tv = Tv0 

    
        else: # pitch-down
    
            if tv > 0 and tv <= Tvl:
                Tv = 0.5*Tv0 
            elif tv > Tvl and tv <= 2*Tvl:
                Tv = 0.25*Tv0 
            else:
                Tv = 0.25*Tv0 

    
    ## ------------------------------------------------------------------------ LEV strength
    
    Cv = CN_C - CN_f 
    
    dCvdt = (Cv-Cv_prev)/dt 
    
    
    ## ------------------------------------------------------------------------ load computation
    
    # normal loads
    
    if tv > 0 and tv <= Tvl and alpha*dCvdt >= 0:
        
        CN_v = CN_v_prev *  np.exp(-ds/Tv) + (Cv-Cv_prev) *  np.exp(-0.5*ds/Tv) 
        
    else:
        
        CN_v = CN_v_prev *  np.exp(-ds/Tv) 
        
    # pitching moment
    
    x_CP = 0.25 * ( 1-np.cos(np.pi*tv/Tvl) ) 
    
    CM_v = - x_CP * CN_v 
        
    # update previous timestep
    
    state[22:24]= [Cv, CN_v] 
    
    return CN_v, CM_v, Tv, state




def BL(alpha, dalphadt, dthetadt, V, M, dt, chord, bsc, x_AC, calibrationData, polarData, formulation, fMode, timeConstantsMod, vortexModule, secondaryVortex, state):

    '''
    BEDDOES-LEISHMAN (OR) Original model - Indicial formulation - v2.6
    
    Closed-loop version
    Secondary vortex shedding
    Missing def for alpha, q computation under generic pitching conditions (i.e. arbitrary pitching point and/or plunging) 
    
    References:
     - Leishman J. G. and Beddoes T. S. A Generalised Model for Airfoil Unsteady Aerodynamic Behaviour and Dynamic Stall Using the Indicial Method
       Proceedings of the 42nd Annual Forum of the American Helicopter Society, Washington, DC, 1986
     - Leishman J. G. and Crouse G. L. Jr. State-Space Model for Unsteady Airfoil Behaviour and Dynamic Stall. Proceedings of the
       AIAA/AHS/ASME Structural Dynamics and Materials Conference, Mobile, Alabama, April 1989. AIAA paper 89-1319
     - Leishman J. G. and Beddoes T. S. A Semi-Empirical Model for Dynamic Stall. Journal of the American Helicopter Society, 34:3–17, 1989
       Chantharasenawong C., Nonlinear Aeroelastic Behaviour of Aerofoils Under Dynamic Stall, PhD thesis, Imperial College London, 2007
    '''

    ## Initialisation
    
    CN_v = 0 
    CM_v = 0
    
    f_lag_prev = state[24] 
    dfdt = state[25] 
    tv = state[26] 
    f_LEV = state[27] 
    
    
    ## constants from polar data
    
    # attached flow
    
    A1 = calibrationData[0] 	   							                    # constant of the indicial response circulatory CN - aoa [-] 
    b1 = calibrationData[1] 								                    # constant of the indicial response circulatory CN - aoa [-] 
    A2 = calibrationData[2]  							                        # constant of the indicial response circulatory CN - aoa [-] 
    b2 = calibrationData[3] 								                    # constant of the indicial response circulatory CN - aoa [-] 
    A3 = calibrationData[4] 							                        # constant of the indicial response impulsive CM - alpha [-] 
    b3 = calibrationData[5] 								                    # constant of the indicial response impulsive CM - alpha [-] 
    A4 =calibrationData[6]   							                        # constant of the indicial response impulsive CM - alpha [-] 
    b4 = calibrationData[7] 								                    # constant of the indicial response impulsive CM - alpha [-] 
    A5 = calibrationData[8]   							                        # constant of the indicial response impulsive CM - q [-] 
    b5 = calibrationData[9]  							                        # constant of the indicial response impulsive CM - q [-] 
    
    CD0 = calibrationData[10]                                                   # drag coefficient @ zero-lift [-]
    CM0 = calibrationData[11]                                                   # pitching moment coefficient @ zero-lift [-]
    alpha0 = calibrationData[12]                                                # zero-lift angle of attack [rad] [-]
    m_CN = calibrationData[13]                                                  # normal force coefficient slope [1/rad] 
    TP = calibrationData[14]                                                    # time constant related to LE pressure dynamics [-]
    eta = calibrationData[15]                                                   # chordwise force recovery factor [-] 
    
    # separated flow
    
    F1 = calibrationData[16]                                                    # constant regulating the speed difference between separation point and center of pressure [-] 
    delta_alpha1 = calibrationData[17]                                          # constant regulating the shift of the separation point towards the LE after shedding [-]
    
    alpha10 = calibrationData[18]                                               # Leishman fitting - separation point - positive breakpoint angle [rad] - needs to be taken as absolute value
    S1 = calibrationData[19]                                                    # Leishman fitting - separation point - positive constant [rad]
    S2 = calibrationData[20]                                                    # Leishman fitting - separation point - positive constant [rad]
    
    alpha20 = calibrationData[21]                                               # Leishman fitting - separation point - negative breakpoint angle [rad] - needs to be taken as absolute value
    S3 = calibrationData[22]                                                    # Leishman fitting - separation point - negative constant [rad]
    S4 = calibrationData[23]                                                    # Leishman fitting - separation point - negative constant [rad]
    
    K0 = calibrationData[24]                                                    # Leishman fitting - center of pressure - constant [-]
    K1 =calibrationData[25]                                                     # Leishman fitting - center of pressure - constant [-]
    K2 = calibrationData[26]                                                    # Leishman fitting - center of pressure - constant [-] 
    m = calibrationData[27]                                                     # Leishman fitting - center of pressure - constant [-] 
    
    CN1 = calibrationData[28]                                                   # positive critical normal load [-] - needs to be taken as absolute value
    CN2 = calibrationData[29]                                                   # negative critical normal load [-] - needs to be taken as absolute value
    
    Tf0 = calibrationData[30]                                                   # time constant related to the unsteady boundary layer dynamics [-]
    
    # vortex shedding
    
    Tv0 = calibrationData[31]                                                   # time constant associated to LEV decay [-]
    Tvl = calibrationData[32]                                                   # characteristic time in semi-chordsrequired to the LEV to go from LE to TE [-]
    Str = calibrationData[33]                                                   # LEV Strouhal number [-]
    Df = calibrationData[34]                                                    # constant in the computation of CC during vortex shedding [-]
    k_CC = calibrationData[35]                                                  # constant in the computation of CC during vortex shedding [-]
    
    ## derived quantities
    
    ds = 2*V*dt/chord                                                           # non-dimensional timestep [-]
    
    # select constants based on the sign of alpha
    
    if alpha >= alpha0:
    
        alpha1 = alpha10 
    
    else:
    
        alpha10 = np.abs(alpha20) 
        alpha1 = alpha10 
        S1 = S3 
        S2 = S4 
        CN1 = np.abs(CN2) 
    
    Tf = Tf0 
    Tv = Tv0 
    
    
    ## apply BL model
    
    # attached flow module
    
    if formulation == 'compressible':
    
        CN_alpha_c, CN_q_c, CN_I, CN_lag, alpha_lag, CC_pot, CM_alpha_c, CM_q_c, CM_I, state = BL_attachedFlow(alpha, dthetadt, V, M, dt, chord, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, A3, b3, A4, b4, A5, b5, state) 
    
        CN_C = CN_alpha_c + CN_q_c
    
    elif formulation == 'incompressible':
    
        CN_C, CN_I, CN_lag, alpha_lag, CC_pot, CM_alpha_c, CM_q_c, CM_I, ds, state = BL_attachedFlow_incompressible(alpha, dthetadt, V, dt, chord, bsc, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, state) 
     
        
    CM_C = CM_alpha_c + CM_q_c 
    
    # compute characteristics of unsteady boundary layer
    
    AOA, f_raw, x_CP_raw = evaluatePolar(alpha0, m_CN, CM0, polarData) 
    
    if timeConstantsMod == 'on':
    
        Tf, alpha1 = BL_computeTimeConstants(alpha, dalphadt, tv, f_lag_prev, CN_lag, alpha10, delta_alpha1, CN1, Tf0, Tvl, dfdt) 
    
    
    f_lag, x_CP, f, fM, state = BL_unsteadyBoundaryLayer(alpha, alpha_lag, ds, alpha1, S1, S2, Tf, K0, K1, K2, m, Tf0, F1, AOA, f_raw, x_CP_raw, fMode, state) 
    
    
    # separated flow module
    
    [CN_f, CC_f, CM_f] = BL_separatedFlow(f_lag, x_CP, CN_C, CC_pot, eta) 
    
    # vortex shedding module
    
    Tst = 2*(1-f_lag)/Str                                                       # secondary vortex characteristic shedding time
    
    if vortexModule == 'on' and np.abs(CN_lag) >= CN1:
    
        if secondaryVortex == 'on' and f_LEV==0 or tv > Tvl + Tst:
            tv = 0 
            f_LEV = 1 
            print('VortexShedding!')

    
        CN_v, CM_v, Tv, state = BL_vortexShedding(alpha, dalphadt, CN_C, CN_f, dt, ds, tv, Tv0, Tvl, timeConstantsMod, state) 
    
        CC_f = k_CC + CC_f * f_lag**( Df*(np.abs(CN_lag)-CN1) + (f_lag-f) ) 
    
        tv = tv + ds 
    
    elif np.abs(CN_lag) < CN1:
    
        tv = 0 
        f_LEV = 0 
    
    
    ## compute airfoil loads
    
    # normal and tangential loads
    
    CN = CN_I + CN_f + CN_v 
    CC = CC_f 
    
    # lift, drag and pitching moment
    
    CL = CN * np.cos(alpha) + CC * np.sin(alpha)  
    CD = CD0 + CN * np.sin(alpha) - CC * np.cos(alpha) 
    CM = CM0 + CM_I + (CM_q_c + CM_f) + CM_v 
    
    
    ## update status
    
    dfdt = f_lag - f_lag_prev 
    tv_output = tv 
    
    state[24:28] = [f_lag, dfdt, tv, f_LEV] 
    
    
    ## output quantities of interest
    
    comp = [CN_C, CN_I, CM0, CM_C, CM_I, CM_f, CM_v]                                # load components
    
    bl = [Tf, Tv, alpha1, f, fM, x_CP]                                              # unsteady boundary layer properties
    
    
    return CN, CC, CL, CD, CM, f_lag, tv_output, comp, bl, state
    
