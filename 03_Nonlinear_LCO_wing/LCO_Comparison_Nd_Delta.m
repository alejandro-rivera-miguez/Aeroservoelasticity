%% =======================================================================
% SCRIPT: LCO_Comparison_Nd_Delta
% DESCRIPTION: Sweeps through different numbers of modes (Nd) and 
%              free-play gaps (delta) to compare the resulting Limit 
%              Cycle Oscillation branches on the Goland wing.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================

clear; clc; close all;

addpath("Auxiliary Functions\")

% --- Input Parameters ---
Nd_vec = [2, 6];               % Number of coupled modes               
delta_vec = [0.1, 0.5, 1, 2];  % Free-play gaps in degrees
colores = lines(length(delta_vec));

% Initialize global comparison figures
figU_norm_comp = figure('Name', 'Normalized Velocity Comparison', 'Color', 'w');
figW_norm_comp = figure('Name', 'Normalized Frequency Comparison', 'Color', 'w');
figU_abs_comp  = figure('Name', 'Absolute Velocity Comparison', 'Color', 'w');
figW_abs_comp  = figure('Name', 'Absolute Frequency Comparison', 'Color', 'w');

% --- Main Loop Generation ---
for n = 1:length(Nd_vec)
    Nd = Nd_vec(n);
    
    % Baseline Flutter solution for scaling
    [UF, wF, qF] = NoKB_FlapPKGoland(100, Nd, 200, 1000, 0);
    
    % Define line style based on structural basis size
    if Nd == 2
        ls = '-';
    else
        ls = '--';
    end
    
    % Initialize individual figures per Nd
    figU_norm = figure('Name', ['Normalized Velocity Nd = ', num2str(Nd)], 'Color', 'w');
    figW_norm = figure('Name', ['Normalized Frequency Nd = ', num2str(Nd)], 'Color', 'w');
    figU_abs  = figure('Name', ['Absolute Velocity Nd = ', num2str(Nd)], 'Color', 'w');
    figW_abs  = figure('Name', ['Absolute Frequency Nd = ', num2str(Nd)], 'Color', 'w');
    
    for d = 1:length(delta_vec)
        delta = delta_vec(d);
        c = colores(d,:);
        
        % Calculate LCO Branch
        [Ulc, wlc, A_v_rad] = Goland_LCO_identification_PK_function(delta, Nd);
        A_v = A_v_rad * 180 / pi; % Convert to degrees
        
        figs = [figU_norm, figW_norm, figU_abs, figW_abs, ...
                figU_norm_comp, figW_norm_comp, figU_abs_comp, figW_abs_comp];
        
        % Data for X and Y axes
        dataX = {Ulc/UF, wlc/wF, Ulc, wlc, Ulc/UF, wlc/wF, Ulc, wlc};
        dataY = {A_v/delta, A_v/delta, A_v, A_v, A_v/delta, A_v/delta, A_v, A_v};
        start_pt_X = {1, 1, UF, wF, 1, 1, UF, wF};
        
        % Legend entry 
        leg_entry = ['$\delta = ', num2str(delta), '^\circ$ ($N_d=', num2str(Nd), '$)'];
        
        % Loop through all 8 target figures
        for i = 1:8
            figure(figs(i)); hold on;
            % Flutter point connection line (dashed reference)
            plot([start_pt_X{i}, dataX{i}(1)], [0, dataY{i}(1)], 'LineStyle', ls, 'LineWidth', 1.5, 'Color', c, 'HandleVisibility', 'off');
            % Actual LCO Curve
            plot(dataX{i}, dataY{i}, 'LineStyle', ls, 'LineWidth', 1.5, 'Color', c, 'DisplayName', leg_entry);
        end
    end
    
    % Format individual figures
    format_figs([figU_norm, figW_norm, figU_abs, figW_abs], Nd);
end

% Format aggregate comparison figures
format_figs([figU_norm_comp, figW_norm_comp, figU_abs_comp, figW_abs_comp], '2 \textrm{ vs } 6');

%% =======================================================================
% LOCAL HELPER FUNCTION: format_figs
% DESCRIPTION: Applies standard LaTeX formatting to arrays of figures.
% =======================================================================
function format_figs(handles, Nd_val)
    x_labels = {'$U/U_F$', '$\omega/\omega_F$', '$U$ [m/s]', '$\omega$ [rad/s]'};
    y_labels = {'$A/\delta$ [$^\circ$]', '$A/\delta$ [$^\circ$]', '$A$ [$^\circ$]', '$A$ [$^\circ$]'};
    % titles   = {'Normalized Velocity', 'Normalized Frequency', 'Absolute Velocity', 'Absolute Frequency'};
    
    if isnumeric(Nd_val)
        txt_nd = num2str(Nd_val);
    else
        txt_nd = Nd_val; 
    end
    
    for i = 1:4
        figure(handles(i));
        grid on; grid minor;
        xlabel(x_labels{i}, 'Interpreter', 'latex', 'FontSize', 16);
        ylabel(y_labels{i}, 'Interpreter', 'latex', 'FontSize', 16);
        legend('show', 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    end
end