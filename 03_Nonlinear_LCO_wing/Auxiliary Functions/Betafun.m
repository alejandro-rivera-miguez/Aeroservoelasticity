%% =======================================================================
% FUNCTION: Betafun
% DESCRIPTION: Returns the aileron/flap deflection (beta) value given 
%              the corresponding structural mode vector 'q'.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
% =======================================================================
function val = Betafun(y, q)
    yL = 3.35;         % [m]
    yU = 3.35 + 1.8;   % [m]
    val = nan(size(y));
    for i = 1:length(y)
        if yL <= y(i) && y(i) <= yU
            val(i) = q(end);
        else
            val(i) = 0;
        end
    end
end