%% =======================================================================
% FUNCTION: elem_T_matrices
% DESCRIPTION: Mass and stiffness matrices for a FEM torsion element of 
%              length L with constant moment of inertia Jp and constant 
%              stiffness GJ.
%              The DOFs are: theta(1), theta(2).
%              The shape functions are linear functions.
%              (Ref: Cooper Wright Pag 397 - 400)
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   GJ - Torsional stiffness constant of the element
%   Jp - Polar moment of inertia (mass constant) of the element
%   L  - Length of the element
%
% OUTPUTS:
%   M  - Element Mass matrix (2x2)
%   K  - Element Stiffness matrix (2x2)
%% =======================================================================
function [M, K] = elem_T_matrices(GJ, Jp, L) 

    M = (Jp * L / 6) * [ 2,  1;
                         1,  2 ];
    
    K = (GJ / L) * [  1, -1;
                     -1,  1 ];
end