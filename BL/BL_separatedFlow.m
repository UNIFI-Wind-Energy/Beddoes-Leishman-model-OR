function [CN_f, CC_f, CM_f] = BL_separatedFlow(f_lag, x_CP, alpha0, alphaE, CN_C, m_CN, eta)

% BEDDOES-LEISHMAN - SEPARATED FLOW MODULE Computation of blade loads
% during TE separation

%% ------------------------------------------------------------------------ unsteady loads

% normal load

CN_f = CN_C * ( (1+sqrt(f_lag))/2 )^2;

% tangential load

CC_f = eta * m_CN * (alphaE-alpha0) * tan(alphaE+alpha0) * sqrt(f_lag);

% pitching moment
    
CM_f = CN_f * x_CP;
   
end

