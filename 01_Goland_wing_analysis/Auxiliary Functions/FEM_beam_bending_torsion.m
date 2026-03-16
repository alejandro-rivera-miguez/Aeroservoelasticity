%% =======================================================================
% FUNCTION: FEM_beam_bending_torsion
% DESCRIPTION: Computes the coupled bending-torsion modes of a structure.
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   x     - Grid vector of nodes (non-dimensional coordinates, [0, 1])
%   EJ    - Vector of bending stiffness properties for each element
%   m     - Vector of mass properties for each element
%   GJ    - Vector of torsional stiffness properties for each element
%   Jp    - Vector of polar moment of inertia properties for each element
%   xt    - Vector of static unbalance (CG offset) for coupling
%   L     - Total length of the beam
%   cdofs - List of constrained DOFs
%   No    - Number of modes to output
%
% OUTPUTS:
%   f     - Coupled frequencies from 1 to No
%   Modes - Modal shapes (displacements and rotations) at nodes
%% =======================================================================
function [f, Modes] = FEM_beam_bending_torsion(x, EJ, m, GJ, Jp, xt, L, cdofs, No)
 
    n_elem = length(x) - 1;
    n_nodes = n_elem + 1;
    n_dofs = n_nodes * 3;
    
    if (x(end) > 1)
        x = x ./ x(end);
        warning('Nodes position is rescaled to be in the interval [0,1].');
    end
    
    if (length(EJ) ~= n_elem) || (length(m) ~= n_elem) || (length(GJ) ~= n_elem) || (length(Jp) ~= n_elem)
        error('Wrong size of input property vectors.');
    end
    
    K = zeros(n_dofs, n_dofs);
    M = zeros(n_dofs, n_dofs);
    
    % Assembly for bending
    for i = 1:n_elem
        [Me, Ke] = elem_B_matrices(EJ(i), m(i), (x(i+1) - x(i))*L);
        id = [(i-1)*3+1, (i-1)*3+2, i*3+1, i*3+2]; 
        K(id, id) = K(id, id) + Ke;
        M(id, id) = M(id, id) + Me;
    end
    
    % Assembly for torsion
    for i = 1:n_elem
        [Me, Ke] = elem_T_matrices(GJ(i), Jp(i), (x(i+1) - x(i))*L);
        id = [3*i, 3*(i+1)];
        K(id, id) = K(id, id) + Ke;
        M(id, id) = M(id, id) + Me;
    end
    
    % Assembly of the inertial coupling term (Static Unbalance)
    for i = 1:n_elem
        Le = (x(i+1) - x(i)) * L;
        Me = m(i) * xt(i) * [ 7*Le/20,   3*Le/20;
                              Le^2/20,   Le^2/30;
                              3*Le/20,   7*Le/20;
                             -Le^2/30,  -Le^2/20 ];
                         
        id_l = [(i-1)*3+1, (i-1)*3+2, i*3+1, i*3+2];
        id_r = [3*i, 3*(i+1)];
        
        M(id_l, id_r) = M(id_l, id_r) + Me;
        M(id_r, id_l) = M(id_r, id_l) + Me';                  
    end
    
    % Apply constraints
    dofs = 1:n_dofs;
    f_dofs = setdiff(dofs, cdofs); 
    
    Kr = K(f_dofs, f_dofs);
    Mr = M(f_dofs, f_dofs);
    
    % Solve Eigenvalue Problem (accounting for formulation signs)
    [V, E] = eig(-Kr, Mr);
    [Eo, I] = sort(diag(E), 'descend');
    
    f = imag(sqrt(Eo(1:No)));
    Modes = zeros(length(f), n_dofs);
    Modes(:, f_dofs) = V(:, I(1:No))';
    
    % Scale the modes to unit mass
    Mm = Modes * M * Modes';
    Modes = diag(1 ./ sqrt(diag(Mm))) * Modes;
end