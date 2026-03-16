%% =======================================================================
% FUNCTION: FEM_beam_bending
% DESCRIPTION: Computes the bending modes of a structure given the mass 
%              and stiffness distributions.
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   x     - Grid vector of nodes (non-dimensional coordinates, [0, 1])
%   EJ    - Vector of bending stiffness properties for each element
%   m     - Vector of mass properties for each element
%   L     - Total length of the beam
%   cdofs - List of constrained DOFs (e.g., [1, 2] sets w and w' of node 1 to 0)
%   No    - Number of modes to output
%
% OUTPUTS:
%   f       - Frequencies from 1 to No
%   Modesw  - Modal shapes for displacements at nodes [No x length(x)]
%   Modeswp - Modal shapes for rotations (w') at nodes [No x length(x)]
%% =======================================================================
function [f, Modesw, Modeswp] = FEM_beam_bending(x, EJ, m, L, cdofs, No)
 
    n_elem = length(x) - 1;
    n_nodes = n_elem + 1;
    n_dofs = n_nodes * 2;
    
    if (x(end) > 1)
        x = x ./ x(end);
        warning('Nodes position is rescaled to be in the interval [0,1].');
    end
    
    if (length(EJ) ~= n_elem) || (length(m) ~= n_elem) 
        error('Wrong size of input property vectors.');
    end
    
    K = zeros(n_dofs, n_dofs);
    M = zeros(n_dofs, n_dofs);
    
    % Matrix Assembly
    for i = 1:n_elem
        [Me, Ke] = elem_B_matrices(EJ(i), m(i), (x(i+1) - x(i))*L);
        idx = (i-1)*2+1 : (i+1)*2;
        K(idx, idx) = K(idx, idx) + Ke;
        M(idx, idx) = M(idx, idx) + Me;
    end
    
    % Apply constraints by isolating free DOFs
    dofs = 1:n_dofs;
    f_dofs = setdiff(dofs, cdofs); 
    
    Kr = K(f_dofs, f_dofs);
    Mr = M(f_dofs, f_dofs);
    
    % Solve Eigenvalue Problem
    [V, E] = eig(Kr, Mr);
    [Eo, I] = sort(diag(E), 'ascend');
    
    f = sqrt(Eo(1:No));
    Modes = zeros(length(f), n_dofs);
    Modes(:, f_dofs) = V(:, I(1:No))';
    
    % Scale the modes to unit mass
    Mm = Modes * M * Modes';
    Modes = diag(1 ./ sqrt(diag(Mm))) * Modes;
    
    Modesw  = Modes(1:No, 1:2:end);
    Modeswp = Modes(1:No, 2:2:end);
end