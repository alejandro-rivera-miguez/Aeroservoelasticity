%% =======================================================================
% MAIN SCRIPT: Workshop 2 - Puma Blade Aeroelasticity
% DESCRIPTION: Computes and compares stability maps for the Puma helicopter 
%              blade in hover using different modal basis and aerodynamic 
%              assumptions (Quasi-steady, Unsteady, 2 vs 6 modes).
%
% AUTHOR: Alejandro Rivera Míguez
%         (Based on course material by G. Quaranta)
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%% =======================================================================
clc;
clear all; close all;

addpath("Auxiliary Functions\")

%% Task 1: Quasi-Steady Analysis
fprintf('\n--- Running Task 1: Quasi-Steady Analysis (2 Modes) ---\n');
QuasiSteady = true;
PumaBladeFEM_1and2(QuasiSteady);

%% Task 2: Unsteady Analysis
fprintf('\n--- Running Task 2: Unsteady Analysis (2 Modes) ---\n');
QuasiSteady = false;
PumaBladeFEM_1and2(QuasiSteady);

%% Task 3: Stability Map (2 Modes)
fprintf('\n--- Running Task 3: Stability Map (2 Modes) ---\n');
% --- 1. CONFIGURATION ---
res_stiff = 100; 
res_cg    = 100; 

% Ranges
link_stiff = linspace(100, 33032, res_stiff); 
cg_target_vec = linspace(24, 41, res_cg);
Data = puma_data(); 

stability_map = zeros(length(cg_target_vec), length(link_stiff));
damping_map = zeros(length(cg_target_vec), length(link_stiff), 2);

% Sweep loop
for i = 1:length(link_stiff)
    fprintf('Processing pitch link stiffness step %d / %d...\n', i, length(link_stiff));
    for j = 1:length(cg_target_vec)
        [is_stable, eig_vals] = PumaBladeFEM_3(link_stiff(i), cg_target_vec(j));
        stability_map(j, i) = is_stable;
        damping_map(j, i, :) = eig_vals;
    end
end

% --- Professional Color Configuration ---
color_unstable = [0.85, 0.32, 0.30]; % Coral Red
color_stable   = [0.40, 0.75, 0.65]; % Pastel Green
color_map = [color_unstable; color_stable];

X_grid = cg_target_vec;
Y_grid = (link_stiff / 33032) * 100;
Z_grid = stability_map'; % 1 = Stable, 0 = Unstable

figure('Name', 'Stability Map (Task 3)', 'Color', 'w', 'Position', [100, 100, 700, 500]);
% 1. Smooth fill guaranteeing both colors
contourf(X_grid, Y_grid, Z_grid, [0, 0.5, 1], 'LineColor', 'none');
hold on;
% 2. Exact separator line at threshold 0.5
contour(X_grid, Y_grid, Z_grid, [0.5, 0.5], 'k', 'LineWidth', 2);
% 3. Apply colors
colormap(color_map);
if ~verLessThan('matlab', '9.12')
    clim([0 1]); 
else
    caxis([0 1]); 
end

% --- Formatting and Aesthetics ---
grid on; grid minor;
set(gca, 'Layer', 'top', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
ax = gca; ax.GridAlpha = 0.4;
xlabel('Blade CoG [\% Chord]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Pitch Link Stiffness [\% Nominal]', 'Interpreter', 'latex', 'FontSize', 14);
cbar = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Unstable', 'Stable'});
set(cbar, 'TickLabelInterpreter', 'latex', 'FontSize', 12);


%% Task 4: Stability Map (6 Modes)
fprintf('\n--- Running Task 4: Stability Map (6 Modes) ---\n');
% --- 1. CONFIGURATION ---
res_stiff2 = 100; 
res_cg2    = 100; 

% Ranges
link_stiff2 = linspace(100, 33032, res_stiff2); 
cg_target_vec2 = linspace(24, 41, res_cg2);

stability_map2 = zeros(length(cg_target_vec2), length(link_stiff2));
damping_map2 = zeros(length(cg_target_vec2), length(link_stiff2), 6);

% Sweep loop 
for i = 1:length(link_stiff2)
    fprintf('Processing pitch link stiffness step %d / %d...\n', i, length(link_stiff2));
    for j = 1:length(cg_target_vec2)
        [is_stable, eig_vals] = PumaBladeFEM_4(link_stiff2(i), cg_target_vec2(j));
        stability_map2(j, i) = is_stable;
        damping_map2(j, i, :) = eig_vals;
    end
end

X_grid2 = cg_target_vec2;
Y_grid2 = (link_stiff2 / 33032) * 100;
Z_grid2 = stability_map2'; 

figure('Name', 'Stability Map (Task 4)', 'Color', 'w', 'Position', [100, 100, 700, 500]);
contourf(X_grid2, Y_grid2, Z_grid2, [0, 0.5, 1], 'LineColor', 'none');
hold on;
contour(X_grid2, Y_grid2, Z_grid2, [0.5, 0.5], 'k', 'LineWidth', 2);
colormap(color_map);
if ~verLessThan('matlab', '9.12'), clim([0 1]); else, caxis([0 1]); end

% --- Formatting and Aesthetics ---
grid on; grid minor;
set(gca, 'Layer', 'top', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
ax = gca; ax.GridAlpha = 0.4;
xlabel('Blade CoG [\% Chord]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Pitch Link Stiffness [\% Nominal]', 'Interpreter', 'latex', 'FontSize', 14);
cbar = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Unstable', 'Stable'});
set(cbar, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% --- COMPARATIVE PLOT (TASK 3 vs TASK 4) ---
figure('Name', 'Stability Comparison: 2 vs 6 Modes', 'Color', 'w', 'Position', [150, 150, 750, 500]);
hold on;
% 1. Draw Task 3 boundary (2 modes) in Dark Blue
[C1, h1] = contour(X_grid, Y_grid, Z_grid, [0.5, 0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', '2 Modes (1st Flap \& Torsion)');
% 2. Draw Task 4 boundary (6 modes) in Dark Orange
[C2, h2] = contour(X_grid2, Y_grid2, Z_grid2, [0.5, 0.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', '6 Modes (Full Basis)');

% --- Formatting and Aesthetics ---
grid on; grid minor;
set(gca, 'Layer', 'top', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
ax = gca; ax.GridAlpha = 0.4;
xlabel('Blade CoG [\% Chord]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Pitch Link Stiffness [\% Nominal]', 'Interpreter', 'latex', 'FontSize', 14);
xlim([24 41]); ylim([0 100]);
legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);


%% Task 5: Exact Unsteady Integration
fprintf('\n--- Running Task 5: Stability Map (Exact Unsteady) ---\n');
% --- 1. CONFIGURATION ---
res_stiff3 = 100; 
res_cg3    = 100; 

link_stiff3 = linspace(100, 33032, res_stiff3); 
cg_target_vec3 = linspace(24, 41, res_cg3);

stability_map3 = zeros(length(cg_target_vec3), length(link_stiff3));
damping_map3 = zeros(length(cg_target_vec3), length(link_stiff3), 6);

for i = 1:length(link_stiff3)
    fprintf('Processing pitch link stiffness step %d / %d...\n', i, length(link_stiff3));
    for j = 1:length(cg_target_vec3)
        [is_stable, eig_vals] = PumaBladeFEM_5(link_stiff3(i), cg_target_vec3(j));
        stability_map3(j, i) = is_stable;
        damping_map3(j, i, :) = eig_vals;
    end
end

X_grid3 = cg_target_vec3;
Y_grid3 = (link_stiff3 / 33032) * 100;
Z_grid3 = stability_map3'; 

figure('Name', 'Stability Map (Task 5)', 'Color', 'w', 'Position', [100, 100, 700, 500]);
contourf(X_grid3, Y_grid3, Z_grid3, [0, 0.5, 1], 'LineColor', 'none');
hold on;
contour(X_grid3, Y_grid3, Z_grid3, [0.5, 0.5], 'k', 'LineWidth', 2);
colormap(color_map);
if ~verLessThan('matlab', '9.12'), clim([0 1]); else, caxis([0 1]); end

% --- Formatting and Aesthetics ---
grid on; grid minor;
set(gca, 'Layer', 'top', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
ax = gca; ax.GridAlpha = 0.4;
xlabel('Blade CoG [\% Chord]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Pitch Link Stiffness [\% Nominal]', 'Interpreter', 'latex', 'FontSize', 14);
cbar = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Unstable', 'Stable'});
set(cbar, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% --- FINAL COMPARATIVE PLOT (ALL TASKS) ---
figure('Name', 'Final Stability Comparison', 'Color', 'w', 'Position', [150, 150, 750, 500]);
hold on;
[C1, h1] = contour(X_grid, Y_grid, Z_grid, [0.5, 0.5], 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', '2 Modes reference section');
[C2, h2] = contour(X_grid2, Y_grid2, Z_grid2, [0.5, 0.5], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', '6 Modes reference section');
[C3, h3] = contour(X_grid3, Y_grid3, Z_grid3, [0.5, 0.5], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 2, 'DisplayName', '6 Modes exact unsteady');

% --- Formatting and Aesthetics ---
grid on; grid minor;
set(gca, 'Layer', 'top', 'FontSize', 12, 'TickLabelInterpreter', 'latex');
ax = gca; ax.GridAlpha = 0.4;
xlabel('Blade CoG [\% Chord]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Pitch Link Stiffness [\% Nominal]', 'Interpreter', 'latex', 'FontSize', 14);
xlim([24 41]); ylim([0 100]);
legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);