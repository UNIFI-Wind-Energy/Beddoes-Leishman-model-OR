function [CN_alpha_c, CN_q_c, CN_I, alphaE, CN_lag, alpha_lag, CM_alpha_c, CM_q_c, CM_I, state] = BL_attachedFlow_incompressible(alpha, dthetadt, dhdt, V, dt, chord, x_AC, alpha0, m_CN, TP, state)

% ATTACHED FLOW MODULE Computes potential unsteady loads in the
% Beddoes-Leishman dynamic stall model

%% initialisation

% update state

alpha_prev = state(1);
q_prev = state(2);
dhdt_prev = state(3);
X1_prev = state(5);
Y1_prev = state(6);
X2_prev = state(8);
Y2_prev = state(9);

DP_prev = state(15);
CN_prev = state(16);

% calibration constants

A1 = 0.165;
b1 = 0.0455;
A2 = 0.335;
b2 = 0.3;

% derived quantities

ds = 2*V*dt/chord;                                                          % non-dimensional timestep [-]

q = dthetadt*chord/V;                                                       % non-dimensional pitch rate [-]

Delta_alpha = alpha - alpha_prev;
Delta_q = q - q_prev;
Delta_h = dhdt - dhdt_prev;


%% -------------------------------------------------------------------------- CN

% ---------------------circulatory part

% angle of attack change

X1 = X1_prev * exp(-b1*ds) + A1 * Delta_alpha * exp(-0.5*b1*ds);

Y1 = Y1_prev * exp(-b2*ds) + A2 * Delta_alpha * exp(-0.5*b2*ds);

CN_alpha_c = m_CN * (alpha - X1 - Y1 - alpha0);

% pitch rate change

X2 = X2_prev * exp(-b1*ds) + A1 * Delta_q * exp(-0.5*b1*ds);

Y2 = Y2_prev * exp(-b2*ds) + A2 * Delta_q * exp(-0.5*b2*ds);

qE = q - X2 - Y2;

CN_q_c = 0.5 * m_CN * qE;

% normal load

CN_C = CN_alpha_c + CN_q_c;

% effective angle of attack at 3/4 chord -> computation of CC

alphaE = (alpha - X1 - Y1) + qE/2;

%--------------------- impulsive part

CN_I = m_CN/2 * (Delta_alpha/ds + 1/V * Delta_h/ds + 0.25 * Delta_q/ds);

% ------------------------- total load

CN = CN_C + CN_I;


%% -------------------------------------------------------------------------- LE pressure lag

DP = DP_prev *  exp(-ds/TP) + (CN-CN_prev) * exp(-0.5*ds/TP);

CN_lag = CN - DP;

alpha_lag = alpha0 + CN_lag/m_CN;


%% -------------------------------------------------------------------------- CM

% ------------------------- circulatory part

% angle of attack change

CM_alpha_c = (0.25-x_AC) * CN_alpha_c; 

% pitch rate change

CM_q_c = -m_CN/8 * q/2; 

% ------------------------- impulsive part

CM_I = -m_CN/8 * ( Delta_alpha/ds + 2/V * Delta_h/ds + 5/8 * Delta_q/ds );


%% update previous timestep

% [              1   2      3          4   5  6  7 8  9  10 11 12 13  14 15 16]
% [            alpha q Delta_alpha Delta_q X1 Y1 D X2 Y2 Dq X3 X4 Y4 DMq DP CN];

state(1) = alpha;
state(2) = q;
state(3) = dhdt;
state(5) = X1;
state(6) = Y1;
state(8) = X2;
state(9) = Y2;
state(15) = DP;
state(16) = CN;

end

