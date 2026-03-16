%% =======================================================================
% FUNCTION: AeroIntegrand
% DESCRIPTION: Wrapper function acting as an integrand for computing 
%              aerodynamic modal projection integrals. Handles vector 
%              inputs to be compatible with MATLAB's 'integral' command.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================
function val = AeroIntegrand(y, MyModei, MyModej, matrix)
    val = nan(size(y));
    for i = 1:length(y)
        % Account for the negative sign in psiwj using a transformation matrix
        val(i) = MyModei(y(i))' * matrix * [-1, 0, 0; 0, 1, 0; 0, 0, 1] * MyModej(y(i)); 
    end
end