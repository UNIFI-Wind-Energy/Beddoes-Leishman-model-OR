function [f_lag2, x_CP, f, fM, state] = BL_unsteadyBoundaryLayer(alpha, alpha_lag, ds, alpha1, S1, S2, Tf, K0, K1, K2, m, Tf0, F1, AOA, f_raw, x_CP_raw, fMode, state)

% BOUNDARY LAYER DYNAMIC Compute unsteady separation point position f_lag
% and center of pressure x_CP

%% initialization

f_lag_prev = state(17);
Df_prev = state(18);
fM_prev = state(19);
DfM_prev = state(20);

alpha_lag_prev = state(21);
Dalpha_prev = state(22);

%% ------------------------------------------------------------------------ unsteady separation point (Kirchhoff theory)

if (strcmp(fMode, 'raw'))

    f = interp1(AOA, f_raw, alpha);

    f_lag = interp1(AOA, f_raw, alpha_lag);

elseif (strcmp(fMode, 'fit'))

    f = f_BL(alpha, alpha1, S1, S2);

    f_lag = f_BL(alpha_lag, alpha1, S1, S2);

end

Df = Df_prev * exp(-ds/Tf) + (f_lag-f_lag_prev) * exp(-0.5*ds/Tf);

f_lag2 = f_lag - Df;


%% ------------------------------------------------------------------------ unsteady center of pressure (Kirchhoff theory)

fM = f;

DfM = DfM_prev * exp(-ds/(F1*Tf0)) + (fM-fM_prev) * exp(-0.5*ds/(F1*Tf0));

fM_lag = fM - DfM;

% if raw x_CP data is used, delayed angle of attack instead of separation
% point f'' is used for data query - time constant Tf is the same

Dalpha = Dalpha_prev * exp(-ds/Tf) + (alpha_lag-alpha_lag_prev) * exp(-0.5*ds/Tf);

alpha_lag2 = alpha_lag - Dalpha;


if (strcmp(fMode, 'raw'))

    x_CP = interp1(AOA, x_CP_raw, alpha_lag2);


elseif (strcmp(fMode, 'fit'))

    if (f_lag2 >= fM_lag)

        x_CP = K0 + K1 * (1-f_lag2) + K2 * sin(pi*f_lag2^m);

    else

        x_CP = K0 + K1 * (1-fM_lag) + K2 * sin(pi*fM_lag^m);

    end

end

%% update previous timestep

state(17:22) = [f_lag Df fM DfM alpha_lag Dalpha];

end