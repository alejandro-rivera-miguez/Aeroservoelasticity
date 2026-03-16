%% =======================================================================
% FUNCTION: ExactBeamBendingModes
% DESCRIPTION: Computes the exact analytical bending modes of a uniform 
%              cantilever beam.
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   x  - Vector giving the grid of points where modal displacements 
%        must be output. Range is automatically scaled between 0 and 1.
%   No - Number of modes to compute
%
% OUTPUTS:
%   f  - Vector containing the computed frequencies
%   Nx - Matrix [No x length(x)] containing modal forms along each row
%% =======================================================================
function [f, Nx] = ExactBeamBendingModes(x, No)
    % Length is normalized to 1
    L = 1;
    
    if (x(end) > 1)
        warning('x values should be between 0 and 1. The vector values are rescaled.');
        x = x ./ x(end);
    end
    
    % Analytical frequencies for a Cantilever Uniform Beam
    f_roots = [0.59686; 1.49417; 2.5002; (1/2 * (7:2:((No-4)*2+7)))'];
    f = f_roots.^2 * (pi^2 / L^2);   
    Nx = zeros(length(f), length(x));
    
    for i = 1:length(f)
        beta = sqrt(f(i)); 
        A1 = (-1)^(i-1) * 1/2;
        
        num_eps = (cosh(beta*L) - sinh(beta*L)) + (cos(beta*L) - sin(beta*L));
        den_eps = sinh(beta*L) + sin(beta*L);
        epsilon = num_eps / den_eps;
        
        A2 = A1 * (1 + epsilon);
        Nx(i,:) = A1 * (cosh(beta.*x) - cos(beta.*x)) + A2 * (sin(beta.*x) - sinh(beta.*x));
    end
end