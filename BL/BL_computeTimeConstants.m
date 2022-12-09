function [Tf, alpha1] = BL_computeTimeConstants(alpha, dalphadt, tv, f_lag, CN_lag, alpha10, delta_alpha1, CN1, Tf0, Tvl, dfdt)

% TIME CONSTANTS Modify boundary layer dynamic response based on the system and LEV status 

%% compute dynamic separation point speed Tf

% separation speed is increased if beyond a certain amount of separated flow (Leishman, 1986)

if (f_lag >= 0.7)
    Tf = Tf0;
else
    Tf = 0.5*Tf0;
end

% vortex shedding - presence of the LEV modifies boundary layer behaviour

if (abs(CN_lag) >= CN1) 

    % case #1: separating flow during vortex shedding - motion of the LEV along 
    % the chord accelerates separation in both pitch-up and pitch down

    if (dfdt < 0) 

        if(tv > 0 && tv < Tvl)
            Tf = 0.5 * Tf0;
        elseif (tv >= Tvl && tv < 2*Tvl)
            Tf = Tf0;
        else
            Tf = Tf0;
        end

    % case #2: re-attaching flow during vortex shedding - presence of the LEV 
    % in the near-wake slows-down flow re-attachment 

    else 

        if(tv > 0 && tv < 2*Tvl)
            Tf = 4*Tf0;
        end

    end

% case #3: re-attaching flow after vortex shedding - recovery capacity of
% the boundary layer is restored

elseif(abs(CN_lag) < CN1 && dfdt >= 0)

    Tf = 2*Tf0;

end



%% compute separation angle alpha1

alpha1 = alpha10;

if (abs(CN_lag) >= CN1 && alpha*dalphadt < 0)                               % vortex shedding

    alpha1 = alpha10 - (1-f_lag)^0.25 * delta_alpha1;

end


end

