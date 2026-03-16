%% =======================================================================
% FUNCTION: TheodorsenCoefficients
% DESCRIPTION: Computes the 13 geometric T-coefficients for the unsteady 
%              aerodynamic Theodorsen model of a trailing-edge flap.
%
% INPUTS:
%   c - Non-dimensional position of flap hinge in reference [-1, +1]
%   e - Non-dimensional position of the elastic axis in reference [-1, +1]
%
% OUTPUT:
%   CoefficientsVector - Vector [13x1] of T1 through T13 coefficients
%% =======================================================================
function CoefficientsVector = TheodorsenCoefficients(c, e)
    
    mu = acos(c);
    
    T1  = -(1/3) * sqrt(1 - c^2) * (2 + c^2) + c * mu;
    T2  = c * (1 - c^2) - sqrt(1 - c^2) * (1 + c^2) * mu + c * mu^2;
    T3  = -(1/8 + c^2) * mu^2 + (1/4) * c * sqrt(1 - c^2) * mu * (7 + 2*c^2) ...
          - (1/8) * (1 - c^2) * (5*c^2 + 4);
    T4  = -mu * c * sqrt(1 - c^2);
    T5  = -(1 - c^2) - mu^2 + 2 * c * sqrt(1 - c^2) * mu;
    
    T6  = 0; % T6 is not used in standard 3-DOF formulation, set to zero
    
    T7  = -(1/8 + c^2) * mu + (1/8) * c * sqrt(1 - c^2) * (7 + 2*c^2);
    T8  = -(1/3) * sqrt(1 - c^2) * (2*c^2 + 1) + c * mu;
    T9  = 0.5 * ((1/3) * (sqrt(1 - c^2))^3 + e * T4);
    T10 = sqrt(1 - c^2) + mu;
    T11 = mu * (1 - 2*c) + sqrt(1 - c^2) * (2 - c);
    T12 = sqrt(1 - c^2) * (2 + c) - mu * (2*c + 1);
    T13 = -0.5 * (T7 + (c - e) * T1);
    
    CoefficientsVector = [T1; T2; T3; T4; T5; T6; T7; T8; T9; T10; T11; T12; T13];
end