%% =======================================================================
% FUNCTION: elem_B_matrices
% DESCRIPTION: Mass and stiffness matrices for a FEM beam element of length 
%              L with constant mass m and constant stiffness EJ.
%              The DOFs are: w(1), w'(1), w(2), w'(2).
%              The shape functions are cubic polynomials.
%              (Ref: Cooper Wright Pag 397 - 400)
%
% AUTHOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%
% INPUTS:
%   EJ - Bending stiffness constant of the element
%   m  - Mass per unit length constant of the element
%   L  - Length of the element
%
% OUTPUTS:
%   M  - Element Mass matrix (4x4)
%   K  - Element Stiffness matrix (4x4)
%% =======================================================================
function [M, K] = elem_B_matrices(EJ, m, L) 
    
    M = (m * L / 420) * [ 156,     22*L,     54,    -13*L;
                           22*L,    4*L^2,   13*L,   -3*L^2;
                           54,     13*L,    156,    -22*L;
                          -13*L,   -3*L^2,  -22*L,    4*L^2 ];
                      
    K = (EJ / L^3) * [  12,      6*L,    -12,      6*L;
                         6*L,    4*L^2,   -6*L,    2*L^2;
                        -12,    -6*L,     12,     -6*L;
                         6*L,    2*L^2,   -6*L,    4*L^2 ];
end