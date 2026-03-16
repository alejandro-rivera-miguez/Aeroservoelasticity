%% =======================================================================
% FUNCTION: Free_Play_dSIDFdA
% DESCRIPTION: Implements the derivative with respect to amplitude (A) of 
%              the Sinusoidal Input Describing Function for a free-play 
%              non-linearity of range 2*delta and stiffness k.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================
function dNA = Free_Play_dSIDFdA(A, delta, k)
    dNA = nan(size(A));
    for i = 1:length(A)
        if A(i) <= delta
            dNA(i) = 0;
        else
            dNA(i) = 4 * k * delta / pi / (A(i)^2) * sqrt(1 - (delta/A(i))^2);
        end
    end
end