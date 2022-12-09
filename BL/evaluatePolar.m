function [AOA, f, x_CP] = evaluatePolar(alpha_0, m_CN, CM0, polarData)

% BL - EVALUATE POLAR Compute static characteristics of the boundary layer
% - separation point f and center of pressure x_cp profiles - from input
% polar data

%% polar data pre-processing

AOA = polarData(:,1);
CL = polarData(:,2);
CD = polarData(:,3);
CM = polarData(:,4);

AOA = deg2rad(AOA);

CN = CL .* cos(AOA) + CD .* sin(AOA);

%% compute static separation point f

f = zeros(length(AOA),1);

for i=1:length(AOA)

    f(i) = (2 * sqrt( CN(i)./(m_CN*(AOA(i)-alpha_0)-0.0001) ) - 1).^2;

    % apply boundaries saturation - f c [0 1]

    if (f(i) > 1)

        f(i) = 1 - 1e-12;
    
    elseif (f(i) < 0)
    
        f(i) = 1e-12;
    
    end

end

%% compute static center of pressure xCP

x_CP = (CM-CM0) ./ CN;

end