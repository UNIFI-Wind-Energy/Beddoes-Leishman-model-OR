function [CN_v, CM_v, Tv, state] = BL_vortexShedding(alpha, dalphadt, CN_C, CN_f, dt, ds, tv, Tv0, Tvl, timeConstantsMod, state)

% BEDDOES-LEISHMAN VORTEX SHEDDING MODULE Compute unsteady loads due to
% primary (LEV) and secondary vortex shedding

Cv_prev = state(23);
CN_v_prev = state(24);

%% compute load response time to LEV Tv

Tv = Tv0;

if(strcmp(timeConstantsMod,'on'))

    if (alpha*dalphadt >= 0) % pitch-up

        if(tv > 0 && tv <= Tvl)
            Tv = Tv0;
        elseif (tv > Tvl && tv <= 2*Tvl)
            Tv = 0.5*Tv0;
        else
            Tv = Tv0;
        end

    else % pitch-down

        if(tv > 0 && tv <= Tvl)
            Tv = 0.5*Tv0;
        elseif (tv > Tvl && tv <= 2*Tvl)
            Tv = 0.25*Tv0;
        else
            Tv = 0.25*Tv0;
        end

    end

end

%% ------------------------------------------------------------------------ LEV strength

Cv = CN_C - CN_f;

dCvdt = (Cv-Cv_prev)/dt;

%% ------------------------------------------------------------------------ load computation

% normal loads

if (tv > 0 && tv <= Tvl && alpha*dCvdt >= 0)
    
    CN_v = CN_v_prev * exp(-ds/Tv) + (Cv-Cv_prev) * exp(-0.5*ds/Tv);
    
else
    
    CN_v = CN_v_prev * exp(-ds/Tv);
    
end

% pitching moment

x_CP = 0.25 * ( 1-cos(pi*tv/Tvl) );

CM_v = - x_CP * CN_v;
    
% update previous timestep

state(23:24)= [Cv CN_v];

end

