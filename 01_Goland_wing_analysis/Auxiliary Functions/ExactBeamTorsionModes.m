%% =======================================================================
% FUNCTION: ExactBeamTorsionModes
% DESCRIPTION: Computes the exact analytical torsional modes of a beam.
%              Frequencies are computed for GJ/(Jp L^2) = 1.
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft (2021)
%
% INPUTS:
%   type - Selects the type of constraints (e.g., 1 for Constrained-Free)
%   x    - Vector giving the grid of points where modal displacements 
%          must be output. Range is automatically scaled between 0 and 1.
%   No   - Number of modes to compute
%
% OUTPUTS:
%   f    - Vector containing the computed frequencies
%   Nx   - Matrix [No x length(x)] containing modal forms along each row
%% =======================================================================
function [f, Nx] = ExactBeamTorsionModes(type, x, No)
    % Length is normalized to 1
    L = 1;
    
    if (x(end) > 1)
        warning('x values should be between 0 and 1. The vector values are rescaled.');
        x = x ./ x(end);
    end
    
    switch type
        case 1
            % Constrained-Free analytical frequencies
            f = (1:2:(1 + 2*(No - 1))) * pi / (2 * L);   
            Nx = zeros(No, length(x));
            
            for i = 1:No
                beta = f(i); 
                A1 = 1;
                Nx(i,:) = A1 * sin(beta .* x);
            end
            
        otherwise
            disp('Unknown constraint type selected.');
    end
end