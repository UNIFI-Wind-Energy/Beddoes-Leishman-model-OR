function [CN, CC, CL, CD, CM, f_lag, tv_output, comp, bl, state] = BL(alpha, dalphadt, dthetadt, V, M, dt, chord, x_AC, calibrationData, polarData, fMode, timeConstantsMod, vortexModule, secondaryVortex, state)

% BEDDOES-LEISHMAN (OR) Original model - Indicial formulation - v2.3
%
% Closed-loop version
% Secondary vortex shedding
% Missing function for alpha, q computation under generic pitching conditions (i.e. arbitrary pitching point and/or plunging) 
%
% References:
% - Leishman J. G. and Beddoes T. S. A Generalised Model for Airfoil Unsteady Aerodynamic Behaviour and Dynamic Stall Using the Indicial Method
%   Proceedings of the 42nd Annual Forum of the American Helicopter Society, Washington, DC, 1986
% - Leishman J. G. and Crouse G. L. Jr. State-Space Model for Unsteady Airfoil Behaviour and Dynamic Stall. Proceedings of the
%   AIAA/AHS/ASME Structural Dynamics and Materials Conference, Mobile, Alabama, April 1989. AIAA paper 89-1319
% - Leishman J. G. and Beddoes T. S. A Semi-Empirical Model for Dynamic Stall. Journal of the American Helicopter Society, 34:3–17, 1989
%   Chantharasenawong C., Nonlinear Aeroelastic Behaviour of Aerofoils Under Dynamic Stall, PhD thesis, Imperial College London, 2007


%% Initialisation

CN_v = 0;
CM_v = 0;

f_lag_prev = state(25);
dfdt = state(26);
tv = state(27);
f_LEV = state(28);


%% constants from polar data

% attached flow

CD0 = calibrationData(1);                                                   % drag coefficient @ zero-lift [-]
CM0 = calibrationData(2);                                                   % pitching moment coefficient @ zero-lift [-]
alpha0 = calibrationData(3);                                                % zero-lift angle of attack [rad] [-]
m_CN = calibrationData(4);                                                  % normal force coefficient slope [1/rad] 
TP = calibrationData(5);                                                    % time constant related to LE pressure dynamics [-]
eta = calibrationData(6);                                                   % chordwise force recovery factor [-] 

% separated flow

F1 = calibrationData(7);                                                    % constant regulating the speed difference between separation point and center of pressure [-] 
delta_alpha1 = calibrationData(8);                                          % constant regulating the shift of the separation point towards the LE after shedding [-]

alpha10 = calibrationData(9);                                               % Leishman fitting - separation point - positive breakpoint angle [rad] - needs to be taken as absolute value
S1 = calibrationData(10);                                                   % Leishman fitting - separation point - positive constant [rad]
S2 = calibrationData(11);                                                   % Leishman fitting - separation point - positive constant [rad]

alpha20 = calibrationData(12);                                              % Leishman fitting - separation point - negative breakpoint angle [rad] - needs to be taken as absolute value
S3 = calibrationData(13);                                                   % Leishman fitting - separation point - negative constant [rad]
S4 = calibrationData(14);                                                   % Leishman fitting - separation point - negative constant [rad]

K0 = calibrationData(15);                                                   % Leishman fitting - center of pressure - constant [-]
K1 =calibrationData(16);                                                    % Leishman fitting - center of pressure - constant [-]
K2 = calibrationData(17);                                                   % Leishman fitting - center of pressure - constant [-] 
m = calibrationData(18);                                                    % Leishman fitting - center of pressure - constant [-] 

CN1 = calibrationData(19);                                                  % positive critical normal load [-] - needs to be taken as absolute value
CN2 = calibrationData(20);                                                  % negative critical normal load [-] - needs to be taken as absolute value

Tf0 = calibrationData(21);                                                  % time constant related to the unsteady boundary layer dynamics [-]

% vortex shedding

Tv0 = calibrationData(22);                                                  % time constant associated to LEV decay [-]
Tvl = calibrationData(23);                                                  % characteristic time in semi-chordsrequired to the LEV to go from LE to TE [-]
Str = calibrationData(24);                                                  % LEV Strouhal number [-]
Df = calibrationData(25);                                                   % constant in the computation of CC during vortex shedding [-]


%% derived quantities

ds = 2*V*dt/chord;                                                          % non-dimensional timestep [-]

% select constants based on the sign of alpha

if(alpha >= alpha0)

    alpha1 = alpha10;

else

    alpha10 = abs(alpha20);
    alpha1 = alpha10;
    S1 = S3;
    S2 = S4;
    CN1 = abs(CN2);

end

Tf = Tf0;
Tv = Tv0;


%% apply BL model

% attached flow module

[CN_alpha_c, CN_q_c, CN_I, alphaE, CN_lag, alpha_lag, CM_alpha_c, CM_q_c, CM_I, state] = BL_attachedFlow(alpha, dthetadt, V, M, dt, chord, x_AC, alpha0, m_CN, TP, state);

CN_C = CN_alpha_c + CN_q_c;
CM_C = CM_alpha_c + CM_q_c;

% compute characteristics of unsteady boundary layer

[AOA, f_raw, x_CP_raw] = evaluatePolar(alpha0, m_CN, CM0, polarData);

if(strcmp(timeConstantsMod,'on'))

    [Tf, alpha1] = BL_computeTimeConstants(alpha, dalphadt, tv, f_lag_prev, CN_lag, alpha10, delta_alpha1, CN1, Tf0, Tvl, dfdt);

end

[f_lag, x_CP, f, fM, state] = BL_unsteadyBoundaryLayer(alpha, alpha_lag, ds, alpha1, S1, S2, Tf, K0, K1, K2, m, Tf0, F1, AOA, f_raw, x_CP_raw, fMode, state);

% separated flow module

[CN_f, CC_f, CM_f] = BL_separatedFlow(f_lag, x_CP, alpha0, alphaE, CN_C, m_CN, eta);

% vortex shedding module

Tst = 2*(1-f_lag)/Str;                                                      % secondary vortex characteristic shedding time

if (strcmp(vortexModule,'on') && abs(CN_lag) >= CN1)

    if(strcmp(secondaryVortex,'on') && (f_LEV==0 || tv > Tvl + Tst))
        tv = 0;
        f_LEV = 1;
        disp('VortexShedding!')
    end

    [CN_v, CM_v, Tv, state] = BL_vortexShedding(alpha, dalphadt, CN_C, CN_f, dt, ds, tv, Tv0, Tvl, timeConstantsMod, state);

    CC_f = CC_f * f_lag^(Df*(abs(CN_lag)-CN1));

    tv = tv + ds;

elseif (abs(CN_lag) < CN1)

    tv=0;
    f_LEV=0;

end


%% compute airfoil loads

% normal and tangential loads

CN = CN_I + CN_f + CN_v;
CC = CC_f;

% lift, drag and pitching moment

CL = CN * cos(alpha) + CC * sin(alpha); 
CD = CD0 + CN * sin(alpha) - CC * cos(alpha);
CM = CM0 + CM_I + (CM_q_c + CM_f) + CM_v;


%% update status

dfdt = f_lag - f_lag_prev;
tv_output = tv;

state(25:28) = [f_lag dfdt tv f_LEV];


%% output quantities of interest

comp = [CN_C CN_I CN_f CN_v CM0 CM_C CM_I CM_f CM_v];                      % load components

bl = [Tf Tv alpha1 f fM x_CP];                                             % unsteady boundary layer properties


end

