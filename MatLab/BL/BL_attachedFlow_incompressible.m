function [CN_C, CN_I, CN_lag, alpha_lag, CC, CM_alpha_c, CM_q_c, CM_I, ds, state] = BL_attachedFlow_incompressible(alpha, theta_dot, V, dt, chord, bsc, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, state)

% ATTACHED FLOW MODULE Computes potential unsteady loads in the
% Beddoes-Leishman dynamic stall model

%% initialisation

% update state

alpha_prev = state(1);
theta_dot_prev = state(2);
V_prev = state(3);
X1_prev = state(4);
Y1_prev = state(5);
w_34_prev = state(6);

DP_prev = state(15);
CN_prev = state(16);

% derived quantities

ds = (V+V_prev)*dt/chord;                                                          % non-dimensional timestep [-]

Delta_alpha = alpha - alpha_prev;
Delta_theta_dot = theta_dot - theta_dot_prev;
Delta_V = V - V_prev;


%% -------------------------------------------------------------------------- CN

% ---------------------circulatory part

% relative downwash at 3/4 chord 
    
w_34 = V * alpha - (bsc - 0.75) * chord * theta_dot;

% deficiency functions - shed vorticity effect

Delta_w_34 = w_34 - w_34_prev;

X1 = X1_prev * exp(-b1*ds) + A1 * Delta_w_34 * exp(-0.5*b1*ds);
Y1 = Y1_prev * exp(-b2*ds) + A2 * Delta_w_34 * exp(-0.5*b2*ds);

w_34_E = w_34 - X1 - Y1;

% effective angle of attack at 3/4 chord -> circulatory loads computation

alphaE = w_34_E / V;

CN_C = m_CN * (alphaE - alpha0);

%--------------------- impulsive part

% relative downwash acceleration at 1/2 chord

w_12_dot = V * Delta_alpha/dt + alpha * Delta_V/dt - (bsc - 0.5) * chord * Delta_theta_dot/dt;

CN_I = m_CN/4 * chord * w_12_dot/V^2;

% ------------------------- total load

CN = CN_C + CN_I;


%% -------------------------------------------------------------------------- LE pressure lag

DP = DP_prev *  exp(-ds/TP) + (CN-CN_prev) * exp(-0.5*ds/TP);

CN_lag = CN - DP;

alpha_lag = alpha0 + CN_lag/m_CN;


%% -------------------------------------------------------------------------- CC
    
CC = m_CN * (alphaE - alpha0) * tan(alphaE);


%% -------------------------------------------------------------------------- CM

% ------------------------- circulatory part

% angle of attack change

CM_alpha_c = (0.25 - x_AC) * CN_C; 

% pitch rate (virtual camber) component

CM_q_c = - m_CN/16 * chord/V * theta_dot;

% ------------------------- impulsive part

CM_I = - 0.25 * CN_I - m_CN/128 * chord^2 / V^2 * Delta_theta_dot/dt;


%% update previous timestep

% [  1        2     3  4  5  6  7  8  9  10 11 12 13  15 16]
% [alpha theta_dot  V  X1 Y1 w_34                     DP CN]

state(1) = alpha;
state(2) = theta_dot;
state(3) = V;
state(4) = X1;
state(5) = Y1;
state(6) = w_34;
state(15) = DP;
state(16) = CN;

end

