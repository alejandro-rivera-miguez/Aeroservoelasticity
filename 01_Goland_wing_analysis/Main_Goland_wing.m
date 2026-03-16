%% =======================================================================
% MAIN SCRIPT: Workshop 1 - Aeroelastic Analysis of the Goland Wing
% DESCRIPTION: Evaluates structural convergence, coupled frequencies,
%              flutter boundaries, parameter sensitivity, tip mass
%              balancing, and compressibility effects.
%
% AUTHOR: Alejandro Rivera Míguez
% PROFESSOR: Giuseppe Quaranta, Politecnico di Milano
% COURSE: Aeroservoelasticity of fixed and rotary wing aircraft
%% =======================================================================

clc;
clear all;
addpath("Auxiliary Functions\")

%% ==========================================
% STYLE CONFIGURATION (RUN AT START)
% ==========================================
% 1. Cleanup
close all; clc;

% 2. LaTeX Texts and Sizes
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 12); % Axes font size
set(groot, 'defaultTextFontSize', 14); % Titles and loose text

% 3. Lines and Graphics
set(groot, 'defaultLineLineWidth', 1.5); % Line width
set(groot, 'defaultAxesBox', 'on');      % Box always on
set(groot, 'defaultAxesXGrid', 'on');    % X Grid on
set(groot, 'defaultAxesYGrid', 'on');    % Y Grid on

% 4. Legend Position 
set(groot, 'defaultLegendLocation', 'northeast');


%% =======================================================================
% TASK 1: Convergence Analysis (Bending & Torsion)
% =======================================================================

% Parameter Definition
EJ = 9.77e6;
m = 35.75;
GJ = 9.876e5;
Jp = 8.65;
L = 6.096;
chord = 1.829;
xt_n = 0.1 * chord;
Cla = 4.2; % Lift slope
cdofs = [1,2];
cdofst = 1;
type = 1;

n_elem0 = 2;
n_elemend = 10;
imax = n_elemend - n_elem0 + 1;
tol = 1e-10;

% --- Bending Loop ---
Modesw = cell(n_elemend - n_elem0, 1);
Modesw_exact = cell(n_elemend - n_elem0, 1);
i = 0;

for n_elem = n_elem0:n_elemend
    i = i + 1;
    x = 0:1/n_elem:1;
    x_av = 1/2 * (x(1:end-1) + x(2:end));
    EJ_v = EJ * x_av.^0;
    m_v = m * x_av.^0;
    
    No = (n_elem + 1) * 2 - 2;
    if No > 10
        No = 10;
    end
    
    [f, Modesw{i}, ~] = FEM_beam_bending(x, EJ_v, m_v, L, cdofs, No);
    [f_exact, Modesw_exact{i}] = ExactBeamBendingModes(x, No);
    
    Modesw{i} = Modesw{i} ./ Modesw{i}(:, end);
    Modesw_exact{i} = Modesw_exact{i} ./ Modesw_exact{i}(:, end);
    
    % Standard L2 Norm approach:
    diff = Modesw_exact{i} - Modesw{i};
    abs_error_L2 = sqrt(trapz(x, diff.^2, 2));
    exact_L2 = sqrt(trapz(x, Modesw_exact{i}.^2, 2));
    
    % Relative Global Error
    error_L2w{i} = abs_error_L2 ./ exact_L2;
end

error_L2w_modes = cell(1, 10);
for modes = 1:10
    for i = 1:imax
        if length(error_L2w{i}) >= modes
            error_L2w_modes{modes}(1, i) = error_L2w{i}(modes, 1);
        end
    end
end

% Plot Convergence Analysis - Bending
x1 = n_elem0:n_elemend;
figure(1); clf; hold on;
color_map = lines(5); 
styles = {'-','-','-','-','-'};
legend_labels = cell(1, 5);

for modes = 1:5
    plot(x1, error_L2w_modes{modes}, 'Color', color_map(modes,:), 'LineStyle', styles{modes});
    legend_labels{modes} = sprintf('Mode $%d$', modes);
end

% 1% Error Reference Line
plot(x1, 0.01 .* ones(size(x1)), 'k--', 'HandleVisibility', 'off', 'LineWidth', 1.0); 

xlabel('Number of elements $n_{elem}$');
ylabel('$L_2$ Error norm $\|e\|_2$');
ylim([0 0.05]);
legend(legend_labels);

% --- Torsion Loop ---
Modest = cell(n_elemend - n_elem0, 1);
Modest_exact = cell(n_elemend - n_elem0, 1);
i = 0;

for n_elem = n_elem0:n_elemend
    i = i + 1;
    x = 0:1/n_elem:1;
    x_av = 1/2 * (x(1:end-1) + x(2:end));
    GJ_v = GJ * x_av.^0;
    Jp_v = Jp * x_av.^0;
    
    No = n_elem;
    if No > 10
        No = 10;
    end
    
    [ft, Modest{i}] = FEM_beam_torsion(x, GJ_v, Jp_v, L, cdofst, No);
    [ft_exact, Modest_exact{i}] = ExactBeamTorsionModes(type, x, No);
    
    Modest{i} = Modest{i} ./ Modest{i}(:, end);
    Modest_exact{i} = Modest_exact{i} ./ (Modest_exact{i}(:, end) + eps);
    
    % Standard L2 Norm approach:
    diff_t = Modest_exact{i} - Modest{i};
    abs_error_L2_t = sqrt(trapz(x, diff_t.^2, 2));
    exact_L2_t = sqrt(trapz(x, Modest_exact{i}.^2, 2));
    
    % Relative Global Error
    error_L2t{i} = abs_error_L2_t ./ exact_L2_t;
end

error_L2_modest = cell(1, 10);
for modes = 1:10
    for i = 1:imax
        if length(error_L2t{i}) >= modes
            error_L2_modest{modes}(1, i) = error_L2t{i}(modes, 1);
        end
    end
end

% Plot Convergence Analysis - Torsion
figure(2); clf; hold on;
for modes = 1:5
    plot(x1, error_L2_modest{modes}, 'Color', color_map(modes,:), 'LineStyle', styles{modes});
end

% 1% Error Reference Line
plot(x1, 0.01 .* ones(size(x1)), 'k--', 'HandleVisibility', 'off', 'LineWidth', 1.0); 

xlabel('Number of elements $n_{elem}$');
ylabel('$L_2$ Error norm $\|e\|_2$'); 
ylim([0 0.05]);
legend(legend_labels);


%% =======================================================================
% TASK 2: Coupled Frequencies Convergence (Richardson Extrapolation)
% =======================================================================

% Nw and Nt needed for frequency precision; First calculation of 'exact' solution
w0 = 6;
t0 = 6;
f = Goland_coupled_using_modal(w0, t0);
f11 = f(1:6);
f = Goland_coupled_using_modal(2*w0, t0);
f12 = f(1:6);
f_ex_w = f12 + (f12 - f11) / 3; % Extrapolate Bending (Fixed Torsion)

f = Goland_coupled_using_modal(w0, 2*t0);
f21 = f(1:6);
f = Goland_coupled_using_modal(2*w0, 2*t0);
f22 = f(1:6);
f_ex_w2 = f22 + (f22 - f21) / 3; % Extrapolate Bending (Fixed Double Torsion)

f_ex = f_ex_w2 + (f_ex_w2 - f_ex_w) / 3; % Extrapolate Torsion (using results above)

% Calculate error with increasing Nw and Ntheta
n_modb = 10;
n_modt = 10;
n_modb_v = 1:n_modb;
n_modt_v = 1:n_modt;

f_2 = zeros(length(n_modb_v), length(n_modt_v), 6);
error_f_2 = f_2;

for i = 1:n_modb
    for j = 1:n_modt
        if i + j >= 6
            aux = Goland_coupled_using_modal(n_modb_v(i), n_modt_v(j));
            f_2(i, j, 1:6) = aux(1:6);
            error_f_2(i, j, 1:6) = abs(squeeze(f_2(i, j, 1:6)) - f_ex) ./ f_ex * 100;
        else
            f_2(i, j, :) = NaN;
            error_f_2(i, j, :) = NaN;
        end
    end
end

% --- Surface Plotting Loop ---
for k = 1:6
    figure(k+2); clf;
    
    % Create surface with contour
    s = surfc(n_modb_v, n_modt_v, 100 * error_f_2(:, :, k));
    s(2).LineWidth = 1.5; % Thicken contour lines
    
    [NB, NT] = meshgrid(n_modb_v, n_modt_v); 
    z_plane = 1; % 1% Error threshold plane
    Z = z_plane * ones(size(NB));
    
    % Draw translucent threshold plane
    hold on; 
    h_plane = surf(NB, NT, Z); 
    h_plane.FaceColor = 'red';      
    h_plane.EdgeColor = 'none';     
    h_plane.FaceAlpha = 0.3;        
    
    xlabel('$N_{\omega}$'); 
    ylabel('$N_{\theta}$');
    zlabel('Relative Error [\%]');
    
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex', 'ZScale', 'log');
    colormap(parula);
    view(135, 30);
    box on; 
end


%% =======================================================================
% TASK 3: Flutter Speed Modal Convergence
% =======================================================================

num_modes = 10;
U_flutter3 = zeros(1, num_modes);

for i = 2:num_modes
    U_flutter3(i) = Goland_flutter_w_Modes(i);
end

% Analytical reference value (Goland Wing standard ~137 m/s)
U_exact = 137.0; 
margin_sup = U_exact * 1.01;
margin_inf = U_exact * 0.99;

x_axis = 2:num_modes;
y_data = U_flutter3(2:num_modes);

% --- Plot Convergence of Flutter ---
figure(9); clf; 
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
hold on; box on; grid on;

% Draw +/- 1% bounds
h_bounds = yline([margin_inf, margin_sup], '--r', 'LineWidth', 2, 'Alpha', 0.5);
set(h_bounds, 'HandleVisibility', 'off'); 

% Analytical Solution
plot([2 num_modes], [U_exact U_exact], 'k--', 'LineWidth', 2, 'DisplayName', 'Analytical Solution');

% Numerical Solution
plot(x_axis, y_data, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2, ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 6, 'DisplayName', 'FEM Estimate');

xlabel('Number of Modes $N_{modes}$', 'Interpreter', 'latex');
ylabel('Flutter Speed $U_f$ [m/s]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
xlim([2 num_modes]);
ylim([min(y_data)*0.98, max(y_data)*1.02]);

% Console Output: Minimal modes required
error_rel = abs(y_data - U_exact) ./ U_exact;
valid_indices = find(error_rel < 0.01);

if ~isempty(valid_indices)
    min_modes = x_axis(valid_indices(1));
    fprintf('------------------------------------------------\n');
    fprintf('RESULT: Minimum number of modes for error < 1%% is: %d\n', min_modes);
    fprintf('Obtained speed: %.4f m/s (Error: %.4f%%)\n', ...
            y_data(valid_indices(1)), error_rel(valid_indices(1)) * 100);
    fprintf('------------------------------------------------\n');
    plot(min_modes, y_data(valid_indices(1)), 'rs', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility','off');
else
    fprintf('With %d modes, error is still not below 1%%.\n', num_modes);
end


%% =======================================================================
% TASK 4: Sensitivity Analysis to Physical Parameters
% =======================================================================

variation_range = 0.98:0.005:1.02;
param_vector = [EJ, GJ, xt_n, Cla];
U_flutter4 = zeros(length(variation_range), 4);

for i = 1:length(variation_range)
    for j = 1:4
        modifier = ones(1, 4);
        modifier(j) = variation_range(i);
        mod_param_vector = param_vector .* modifier;
        U_flutter4(i, j) = Goland_flutter(mod_param_vector);
    end
end

% --- Plot Sensitivity Analysis ---
figure(15); clf; hold on;
var_percent = (variation_range - 1) * 100;

param_legends = { ...
    'Bending Stiffness $EI$', ...
    'Torsion Stiffness $GJ$', ...
    'Static Unbalance $x_{cg}$', ...
    'Lift Slope $C_{L\alpha}$' ...
};

plot_colors = lines(4); 
for j = 1:4
    plot(var_percent, U_flutter4(:, j), 'Color', plot_colors(j,:), 'DisplayName', param_legends{j}); 
end

xline(0, 'k--', 'HandleVisibility', 'off'); % Nominal design reference
xlabel('Parameter Variation [\%]');
ylabel('Flutter Speed $U_f$ [m/s]');
xlim([-2 2]); 
legend('Location', 'north');


%% =======================================================================
% TASK 5: Effect of Tip Mass Balancing
% =======================================================================

% Key cases definitions
st_EA = 0;                % Case A: Elastic Axis (33% chord)
st_5  = (0.05 - 0.33);    % Case B: 5% Chord
st_50 = (0.50 - 0.33);    % Case C: 50% Chord

fprintf('Calculating specific tip mass cases...\n');
U_a = Goland_TipWeight(st_EA * chord);
U_b = Goland_TipWeight(st_5 * chord);
U_c = Goland_TipWeight(st_50 * chord);

% Target speed (+15% of Case A)
U_target = U_a * 1.15;

% Sweep from Leading Edge to Trailing Edge
st_v = linspace(-0.33, 0.67, 100);  
U_sweep = zeros(size(st_v));

fprintf('Calculating tip mass sweep...\n');
for k = 1:length(st_v)
    U_sweep(k) = Goland_TipWeight(st_v(k) * chord);
end

x_percent = (st_v + 0.33) * 100;

% --- Plot Tip Mass Influence ---
figure(16); clf; hold on;
plot(x_percent, U_sweep, 'DisplayName', 'Flutter Speed Trend'); 

plot((st_EA+0.33)*100, U_a, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Case A (Elastic Axis)');
plot((st_5 +0.33)*100, U_b, 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'Case B (5\% Chord)');
plot((st_50+0.33)*100, U_c, 'rd', 'MarkerFaceColor', 'r', 'DisplayName', 'Case C (50\% Chord)');

yline(U_target, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Target +15\\%% ($%.1f$ m/s)', U_target));

valid_indices = find(U_sweep >= U_target);
if ~isempty(valid_indices)
    area_x = x_percent(valid_indices);
    area_y = U_sweep(valid_indices);
    fill([area_x fliplr(area_x)], [area_y ones(size(area_y))*U_target], ...
         'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    [max_U, idx_max] = max(U_sweep);
    plot(x_percent(idx_max), max_U, 'p', 'MarkerSize', 12, ...
         'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Max Flutter Speed');
         
    is_possible_txt = 'YES';
else
    is_possible_txt = 'NO';
end

xlabel('Tip Weight CG Position $x_{cg}$ [\% Chord]', 'Interpreter', 'latex');
ylabel('Flutter Speed $U_f$ [m/s]', 'Interpreter', 'latex');
xlim([0 100]); 
legend('Location', 'northeast', 'NumColumns', 1); 
grid on; box on;

% Console Output
fprintf('\n================ RESULTS TASK 5 ================\n');
fprintf('A) CG at Elastic Axis (33%%):   %.2f m/s\n', U_a);
fprintf('B) CG at 5%% of the chord:      %.2f m/s\n', U_b);
fprintf('C) CG at 50%% of the chord:     %.2f m/s\n', U_c);
fprintf('-------------------------------------------------------\n');
fprintf('Target (+15%% of A):            %.2f m/s\n', U_target);
fprintf('Is it achievable?:             %s\n', is_possible_txt);
if ~isempty(valid_indices)
    fprintf('Valid range approx: between %.1f%% and %.1f%% of the chord\n', ...
        min(x_percent(valid_indices)), max(x_percent(valid_indices)));
end
fprintf('=======================================================\n');


%% =======================================================================
% TASK 6: Compressibility Effects (Mach Number)
% =======================================================================

% Prandtl-Glauert Correction: Beta = sqrt(1 - M^2)
Mach_v = 0:0.05:0.7;
U_flutter6 = zeros(size(Mach_v));
fprintf('\nCalculating Mach number sweep...\n');

for i = 1:length(Mach_v)
    U_flutter6(i) = Goland_flutter_Mach(EJ, GJ, xt_n, Mach_v(i));
end

% --- Plot Compressibility Effects ---
figure(17); clf; hold on;
plot(Mach_v, U_flutter6, '-o', 'Color', [0.8500 0.3250 0.0980], ...
    'MarkerSize', 5, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Compressible Limit'); 

yline(U_flutter6(1), 'k--', 'LineWidth', 1, 'DisplayName', 'Incompressible Baseline ($M=0$)');

xlabel('Mach Number $M_{\infty}$');
ylabel('Flutter Speed $U_f$ [m/s]');
xlim([0 0.7]);
ylim([min(U_flutter6)*0.98, max(U_flutter6)*1.02]); 
legend('Location', 'southwest');

%% =======================================================================
% TASK 7: Aileron Effects (Flutter, Reversal, and Divergence)
% =======================================================================