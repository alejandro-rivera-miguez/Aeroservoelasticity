%% =======================================================================
% SCRIPT: Generate_LCO_Animation
% DESCRIPTION: Example script to trigger the 3D video generation of the 
%              Goland wing undergoing Limit Cycle Oscillations (LCO).
%              It uses the converged results from the continuation solver.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================

addpath("Auxiliary Functions\")

% --- Animation Usage Example ---
% Inputs required by Animate_Goland_Wing:
% 1. qr    : FEM structural matrix (e.g., 44x2 for Nd=2) -> 'qstruct'
% 2. qF    : Complex mode shape vector -> Real part + 1i * Imaginary part
% 3. omega : LCO frequency obtained at the specific evaluation point [rad/s]
% 4. name  : Output video filename

% Select a target point from the amplitude sweep (e.g., point 55)
target_idx = 55;

% Reconstruct the complex mode shape for the selected amplitude
complex_mode_shape = qrlc(:, target_idx) + 1i * qilc(:, target_idx);

% Extract the corresponding frequency
lco_frequency = wlc(target_idx);

% Generate the 3D .avi video
fprintf('\n--- Launching 3D Animation Render ---\n');
Animate_Goland_Wing(qstruct, complex_mode_shape, lco_frequency, 'Goland_LCO_Simulation.avi');S