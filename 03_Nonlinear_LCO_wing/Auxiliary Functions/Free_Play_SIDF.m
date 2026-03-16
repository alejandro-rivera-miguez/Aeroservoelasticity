%% =======================================================================
% FUNCTION: Free_Play_SIDF
% DESCRIPTION: Implements the Sinusoidal Input Describing Function (SIDF) 
%              for an actuator with a free-play non-linearity of range 
%              2*delta and stiffness k.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================
function NA = Free_Play_SIDF(A, delta, k)
    NA = nan(size(A));
    for i = 1:length(A)
        if A(i) <= delta
            NA(i) = 0;
        else
            NA(i) = k * (1 - 2/pi * (asin(delta/A(i)) + delta/A(i) * sqrt(1 - (delta/A(i))^2)));
        end
    end
end