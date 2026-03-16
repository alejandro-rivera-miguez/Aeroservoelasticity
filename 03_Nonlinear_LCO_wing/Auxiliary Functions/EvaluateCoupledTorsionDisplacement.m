%% =======================================================================
% FUNCTION: EvaluateCoupledTorsionDisplacement
% DESCRIPTION: Evaluates the torsional displacement at point 'y' associated 
%              with a mode given by the coefficients in vector 'q'.
%              Uses FEM linear Lagrange shape functions.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
% =======================================================================
function w = EvaluateCoupledTorsionDisplacement(y, q, deltay, Connect)
    [Ne, ~] = size(Connect);
    w = nan(size(y));
    for k = 1:length(y)
        % Identify element
        EID = floor(y(k)/deltay) + 1; 
        if EID == Ne + 1 % Catch fringe case where y=L
            EID = EID - 1;
        end
        % Extract relevant degrees of freedom
        NIDs = Connect(EID, 5:7); 
        
        % Evaluate shape functions
        xi = (y(k) - (EID - 1)*deltay) / deltay; 
        psi1 = 1 - 3*xi + 2*xi^2;
        psi2 = 4*xi - 4*xi^2;
        psi3 = 2*xi^2 - xi;
        
        w(k) = [psi1, psi2, psi3] * q(NIDs); 
    end
end