%% =======================================================================
% FUNCTION: FEM_beam_torsion
% DESCRIPTION: Computes the torsional modes of a structure.
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   x     - Grid vector of nodes (non-dimensional coordinates, [0, 1])
%   GJ    - Vector of torsional stiffness properties for each element
%   Jp    - Vector of mass (polar inertia) properties for each element
%   L     - Total length of the beam
%   cdofs - List of constrained DOFs
%   No    - Number of modes to output
%
% OUTPUTS:
%   f      - Frequencies from 1 to No
%   Modest - Modal shapes for torsional rotation at nodes
%% =======================================================================
function [f, Modest] = FEM_beam_torsion(x, GJ, Jp, L, cdofs, No)
 
    n_elem = length(x) - 1;
    n_nodes = n_elem + 1;
    n_dofs = n_nodes;
    
    if (x(end) > 1)
        x = x ./ x(end);
        warning('Nodes position is rescaled to be in the interval [0,1].');
    end
    
    if (length(GJ) ~= n_elem) || (length(Jp) ~= n_elem) 
        error('Wrong size of input property vectors.');
    end
    
    K = zeros(n_dofs, n_dofs);
    M = zeros(n_dofs, n_dofs);
    
    % Matrix Assembly
    for i = 1:n_elem
        [Me, Ke] = elem_T_matrices(GJ(i), Jp(i), (x(i+1) - x(i))*L);
        K(i:(i+1), i:(i+1)) = K(i:(i+1), i:(i+1)) + Ke;
        M(i:(i+1), i:(i+1)) = M(i:(i+1), i:(i+1)) + Me;
    end
    
    % Apply constraints
    dofs = 1:n_dofs;
    f_dofs = setdiff(dofs, cdofs); 
    
    Kr = K(f_dofs, f_dofs);
    Mr = M(f_dofs, f_dofs);
    
    % Solve Eigenvalue Problem
    [V, E] = eig(Kr, Mr);
    [Eo, I] = sort(diag(E), 'ascend');   % omega^2 sorted in ascending order
    
    f = sqrt(Eo(1:No));                  % Natural frequencies (rad/s)
    Modes = zeros(length(f), n_dofs);
    Modes(:, f_dofs) = V(:, I(1:No))';
    
    % Scale the modes to unit mass
    Mm = Modes * M * Modes';
    Modes = diag(1 ./ sqrt(diag(Mm))) * Modes;
    Modest = Modes;
end