function [f] = f_BL(alpha, alpha1, S1, S2)

% LEISHMAN SEPARATION POINT Exponential fitting of raw separation point
% data

if (abs(alpha) <= alpha1)

    f = 1 - 0.3 * exp( (abs(alpha)-alpha1)/S1 );

else

    f= 0.04 + 0.66 * exp( (alpha1-abs(alpha))/S2 );

end


end

