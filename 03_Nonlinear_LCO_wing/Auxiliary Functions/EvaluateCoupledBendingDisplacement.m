%% =======================================================================
% FUNCTION: EvaluateCoupledBendingDisplacement
% DESCRIPTION: Evaluates the bending displacement at point 'y' associated 
%              with a mode given by the coefficients in vector 'q'.
%              Uses FEM cubic Hermite shape functions.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
% =======================================================================
function w = EvaluateCoupledBendingDisplacement(y, q, deltay, Connect)
    [Ne, ~] = size(Connect);
    w = nan(size(y));
    for k = 1:length(y)
        % Identify element
        EID = floor(y(k)/deltay) + 1; 
        if EID == Ne + 1 % Catch fringe case where y=L
            EID = EID - 1;
        end
        % Extract relevant degrees of freedom
        NIDs = Connect(EID, 1:4); 
        
        % Evaluate shape functions
        xi = (y(k) - (EID - 1)*deltay) / deltay; 
        psi1 = 1 - 3*xi^2 + 2*xi^3;
        psi2 = xi - 2*xi^2 + xi^3;
        psi3 = 3*xi^2 - 2*xi^3;
        psi4 = -xi^2 + xi^3;
        
        w(k) = [psi1, psi2*deltay, psi3, psi4*deltay] * q(NIDs); 
    end
end