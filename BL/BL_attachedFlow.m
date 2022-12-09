function [CN_alpha_c, CN_q_c, CN_I, alphaE, CN_lag, alpha_lag, CM_alpha_c, CM_q_c, CM_I, state] = BL_attachedFlow(alpha, dalphadt, V, M, dt, chord, x_AC, alpha0, m_CN, TP, state)

% ATTACHED FLOW MODULE Computes potential unsteady loads in the
% Beddoes-Leishman dynamic stall model

%% initialisation

% update state

alpha_prev = state(1);
q_prev = state(2);
Delta_alpha_prev = state(3);
Delta_q_prev = state(4);
X1_prev = state(5);
Y1_prev = state(6);
D_prev = state(7);
X2_prev = state(8);
Y2_prev = state(9);
Dq_prev = state(10);
X3_prev = state(11);
X4_prev = state(12);
Y4_prev = state(13);
DMq_prev = state(14);
DP_prev = state(15);
CN_prev = state(16);

% calibration constants

A1 = 0.3;
b1 = 0.14;
A2 = 0.7;
b2 = 0.53;

A3 = 1.5;
b3 = 0.25;
A4 = -0.5;
b4 = 0.1;
A5 = 1;
b5 = 5.0;

% derived quantities

ds = 2*V*dt/chord;                                                          % non-dimensional timestep [-]

q = dalphadt*chord/V;                                                       % non-dimensional pitch rate [-]

a = V/M;                                                                    % speed of sound [m/s]
beta = sqrt(1-M^2);                                                         % Glauert compressibility factor [-]

Delta_alpha = alpha - alpha_prev;
Delta_q = q - q_prev;


%% -------------------------------------------------------------------------- CN

% ---------------------circulatory part

% angle of attack change

X1 = X1_prev * exp(-b1*beta^2*ds) + A1 * Delta_alpha * exp(-0.5*b1*beta^2*ds);

Y1 = Y1_prev * exp(-b2*beta^2*ds) + A2 * Delta_alpha * exp(-0.5*b2*beta^2*ds);

CN_alpha_c = m_CN * (alpha - X1 - Y1 - alpha0);

% pitch rate change

X2 = X2_prev * exp(-b1*beta^2*ds) + A1 * Delta_q * exp(-0.5*b1*beta^2*ds);

Y2 = Y2_prev * exp(-b2*beta^2*ds) + A2 * Delta_q * exp(-0.5*b2*beta^2*ds);

qE = q - X2 - Y2;

CN_q_c = 0.5 * m_CN * qE;

% normal load

CN_C = CN_alpha_c + CN_q_c;

% effective angle of attack at 3/4 chord -> computation of CC

alphaE = (alpha - X1 - Y1) + qE/2;

%--------------------- impulsive part

TI = chord/a;

% angle of attack change

K_alpha = 1 / ( 1-M + 0.5*m_CN*beta^2*M^2 * (A1*b1+A2*b2) );

T_alpha = 0.75 * K_alpha * TI;

D = D_prev * exp(-dt/T_alpha) + (Delta_alpha-Delta_alpha_prev)/dt * exp(-0.5*dt/T_alpha);

CN_alpha_i = 4*T_alpha/M * (dalphadt - D);

% pitch rate change

K_q = 1 / ( 1-M + m_CN*beta^2*M^2 * (A1*b1+A2*b2) );

T_q = 0.75 * K_q * TI;

Dq = Dq_prev * exp(-dt/T_q) + (Delta_q-Delta_q_prev)/dt * exp(-0.5*dt/T_q);

CN_q_i = T_q/M * (Delta_q/dt - Dq);

% impulsive load

CN_I = CN_alpha_i + CN_q_i;

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

X3 = X3_prev * exp(-b5*beta^2*ds) + A5 * Delta_q * exp(-0.5*b5*beta^2*ds); 

CM_q_c= -m_CN/16 * (q - X3); 

% ------------------------- impulsive part

% angle of attack change

KalphaM = (A3*b4+A4*b3) / (b3*b4*(1-M));

X4 = X4_prev * exp(-dt/(b3*KalphaM*TI)) + A3 * Delta_alpha * exp(-0.5*dt/(b3*KalphaM*TI));
Y4 = Y4_prev * exp(-dt/(b4*KalphaM*TI)) + A4 * Delta_alpha * exp(-0.5*dt/(b4*KalphaM*TI));

CM_alpha_i = - 1/M * (X4+Y4); 

% pitch rate change

KqM = 7 / ( 15*(1-M) + 1.5*m_CN*beta^2*M^2*A5*b5 ); 

DMq = DMq_prev * exp(-dt/(KqM*TI)) + (Delta_q-Delta_q_prev)/dt * exp(-0.5*dt/(KqM*TI));

CM_q_i = - 7*KqM*TI/(12*M) * (Delta_q/dt - DMq);

% impulsive load

CM_I = CM_alpha_i + CM_q_i;


%% update previous timestep

% [              1   2      3          4   5  6  7 8  9  10 11 12 13  14 15 16]

state(1:16) = [alpha q Delta_alpha Delta_q X1 Y1 D X2 Y2 Dq X3 X4 Y4 DMq DP CN];

end

