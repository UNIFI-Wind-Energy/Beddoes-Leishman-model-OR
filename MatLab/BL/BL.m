function [CN, CC, CL, CD, CM, f_lag, tv_output, comp, bl, state] = BL(alpha, dalphadt, dthetadt, V, M, dt, chord, bsc, x_AC, calibrationData, polarData, formulation, fMode, timeConstantsMod, vortexModule, secondaryVortex, state)

% BEDDOES-LEISHMAN (OR) Original model - Indicial formulation - v2.6
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

A1 = calibrationData(1); 	   							                    % constant of the indicial response circulatory CN - aoa [-] 
b1 = calibrationData(2); 								                    % constant of the indicial response circulatory CN - aoa [-] 
A2 = calibrationData(3);  							                        % constant of the indicial response circulatory CN - aoa [-] 
b2 = calibrationData(4); 								                    % constant of the indicial response circulatory CN - aoa [-] 
A3 = calibrationData(5);							                        % constant of the indicial response impulsive CM - alpha [-] 
b3 = calibrationData(6);								                    % constant of the indicial response impulsive CM - alpha [-] 
A4 =calibrationData(7);   							                        % constant of the indicial response impulsive CM - alpha [-] 
b4 = calibrationData(8);								                    % constant of the indicial response impulsive CM - alpha [-] 
A5 = calibrationData(9);   							                        % constant of the indicial response impulsive CM - q [-] 
b5 = calibrationData(10);  							                        % constant of the indicial response impulsive CM - q [-] 

CD0 = calibrationData(11);                                                   % drag coefficient @ zero-lift [-]
CM0 = calibrationData(12);                                                   % pitching moment coefficient @ zero-lift [-]
alpha0 = calibrationData(13);                                                % zero-lift angle of attack [rad] [-]
m_CN = calibrationData(14);                                                  % normal force coefficient slope [1/rad] 
TP = calibrationData(15);                                                    % time constant related to LE pressure dynamics [-]
eta = calibrationData(16);                                                   % chordwise force recovery factor [-] 

% separated flow

F1 = calibrationData(17);                                                    % constant regulating the speed difference between separation point and center of pressure [-] 
delta_alpha1 = calibrationData(18);                                          % constant regulating the shift of the separation point towards the LE after shedding [-]

alpha10 = calibrationData(19);                                               % Leishman fitting - separation point - positive breakpoint angle [rad] - needs to be taken as absolute value
S1 = calibrationData(20);                                                   % Leishman fitting - separation point - positive constant [rad]
S2 = calibrationData(21);                                                   % Leishman fitting - separation point - positive constant [rad]

alpha20 = calibrationData(22);                                              % Leishman fitting - separation point - negative breakpoint angle [rad] - needs to be taken as absolute value
S3 = calibrationData(23);                                                   % Leishman fitting - separation point - negative constant [rad]
S4 = calibrationData(24);                                                   % Leishman fitting - separation point - negative constant [rad]

K0 = calibrationData(25);                                                   % Leishman fitting - center of pressure - constant [-]
K1 =calibrationData(26);                                                    % Leishman fitting - center of pressure - constant [-]
K2 = calibrationData(27);                                                   % Leishman fitting - center of pressure - constant [-] 
m = calibrationData(28);                                                    % Leishman fitting - center of pressure - constant [-] 

CN1 = calibrationData(29);                                                  % positive critical normal load [-] - needs to be taken as absolute value
CN2 = calibrationData(30);                                                  % negative critical normal load [-] - needs to be taken as absolute value

Tf0 = calibrationData(31);                                                  % time constant related to the unsteady boundary layer dynamics [-]

% vortex shedding

Tv0 = calibrationData(32);                                                  % time constant associated to LEV decay [-]
Tvl = calibrationData(33);                                                  % characteristic time in semi-chordsrequired to the LEV to go from LE to TE [-]
Str = calibrationData(34);                                                  % LEV Strouhal number [-]
Df = calibrationData(35);                                                   % constant in the computation of CC during vortex shedding [-]
k_CC = calibrationData(36);                                                   % constant in the computation of CC during vortex shedding [-]

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

if(strcmp(formulation,'compressible'))

    [CN_alpha_c, CN_q_c, CN_I, CN_lag, alpha_lag, CC_pot, CM_alpha_c, CM_q_c, CM_I, state] = BL_attachedFlow(alpha, V, M, dt, chord, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, A3, b3, A4, b4, A5, b5, state);

    CN_C = CN_alpha_c + CN_q_c;

elseif(strcmp(formulation,'incompressible'))

    [CN_C, CN_I, CN_lag, alpha_lag, CC_pot, CM_alpha_c, CM_q_c, CM_I, ds, state] = BL_attachedFlow_incompressible(alpha, dthetadt, V, dt, chord, bsc, x_AC, alpha0, m_CN, TP, A1, b1, A2, b2, state);

end

CM_C = CM_alpha_c + CM_q_c;

% compute characteristics of unsteady boundary layer

[AOA, f_raw, x_CP_raw] = evaluatePolar(alpha0, m_CN, CM0, polarData);

if(strcmp(timeConstantsMod,'on'))

    [Tf, alpha1] = BL_computeTimeConstants(alpha, dalphadt, tv, f_lag_prev, CN_lag, alpha10, delta_alpha1, CN1, Tf0, Tvl, dfdt);

end

[f_lag, x_CP, f, fM, state] = BL_unsteadyBoundaryLayer(alpha, alpha_lag, ds, alpha1, S1, S2, Tf, K0, K1, K2, m, Tf0, F1, AOA, f_raw, x_CP_raw, fMode, state);

% separated flow module

[CN_f, CC_f, CM_f] = BL_separatedFlow(f_lag, x_CP, CN_C, CC_pot, eta);

% vortex shedding module

Tst = 2*(1-f_lag)/Str;                                                      % secondary vortex characteristic shedding time

if (strcmp(vortexModule,'on') && abs(CN_lag) >= CN1)

    if(strcmp(secondaryVortex,'on') && (f_LEV==0 || tv > Tvl + Tst))
        tv = 0;
        f_LEV = 1;
        disp('VortexShedding!')
    end

    [CN_v, CM_v, Tv, state] = BL_vortexShedding(alpha, dalphadt, CN_C, CN_f, dt, ds, tv, Tv0, Tvl, timeConstantsMod, state);

    CC_f = k_CC + CC_f * f_lag^( Df*(abs(CN_lag)-CN1) + (f_lag-f) );

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

