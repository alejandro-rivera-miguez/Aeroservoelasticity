%% =======================================================================
% FUNCTION: Goland_flutter_w_Modes
% DESCRIPTION: Computes the flutter speed for the Goland Wing based on a 
%              specified number of retained coupled modes.
%              Includes plotting functionality for the V-g diagram, 
%              frequencies, and damping.
%
% AUTHORS: Giuseppe Quaranta <giuseppe.quaranta@polimi.it>
%          Lucrezia Mapelli <lucrezia.mapelli@polimi.it>
%          Aeroservoelasticity of Fixed and Rotary wing Aircraft (2025)
%
% INPUT:
%   No    - Number of coupled modes to retain
% OUTPUT:
%   Uflut - Computed flutter speed [m/s]
%% =======================================================================
function Uflut = Goland_flutter_w_Modes(No)
    
    % Number of elements to be used
    n_elem = 100;
    
    % Uniform node grids 
    x = 0:1/n_elem:1;
    x_av = 1/2 * (x(1:end-1) + x(2:end));
    
    %% Step 1: Compute the bending modes
    % Properties of the beam
    EJ = 9.77e6 * x_av.^0;
    m  = 35.75 * x_av.^0;
    L  = 6.096;
    
    % Properties for torsion
    GJ = 9.876e5 * x_av.^0;
    Jp = 8.65 * x_av.^0;
    cdofs = [1, 2];
    
    % Number of modes to be computed 
    % (Modes are scaled to unit mass)
    n_modb = 20;
    [fb, Modesw, Modeswp] = FEM_beam_bending(x, EJ, m, L, cdofs, n_modb);
    
    %% Step 2: Compute the torsion modes
    n_modt = 20;
    cdofs_t = 1;
    [ft, Modest] = FEM_beam_torsion(x, GJ, Jp, L, cdofs_t, n_modt);
    
    %% Step 3: Assembly of the cross-mass terms
    % Modes are known at node positions, so mass and CG position 
    % must be evaluated there.
    chord = 1.829;
    xt_n = 0.1 * chord * x.^0;
    m_n  = 35.75 * x.^0;
    
    Mbt = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            Mbt(i,j) = trapz(x, Modesw(i, :) .* m_n .* xt_n .* Modest(j, :)) * L;
        end
    end
    
    %% Step 4: Assembly of the global matrices
    % DOFs: 'n' amplitudes of bending modes + 'n' amplitudes of torsional modes 
    n_mod = n_modb + n_modt;
    K = zeros(n_mod, n_mod);
    M = eye(n_mod, n_mod);
    
    K(1:n_modb, 1:n_modb) = diag(fb.^2);
    K(n_modb+1:n_mod, n_modb+1:n_mod) = diag(ft.^2);
    
    M(1:n_modb, n_modb+1:end) = -Mbt;
    M(n_modb+1:end, 1:n_modb) = -Mbt';
    
    %% Step 5: Computation of coupled eigenvalues and eigenvectors
    [V, E] = eig(-K, M);
    [Eo, I] = sort(diag(E), 'descend');
    f = imag(sqrt(Eo));
    V = V(:, I);
    
    Modes = zeros(No, (n_elem+1)*3);
    
    % Normalize modes with unit modal mass
    Modm = V' * M * V;
    V = V * diag(1 ./ sqrt(diag(Modm)));
    
    Modes(:, 1:3:end) = (V(1:n_modb, 1:No)') * Modesw(1:n_modb, :);
    Modes(:, 2:3:end) = (V(1:n_modb, 1:No)') * Modeswp(1:n_modb, :);
    Modes(:, 3:3:end) = (V(n_modb+1:n_mod, 1:No)') * Modest(1:n_modb, :);
    
    %% Step 6: Aeroelastic Model Setup
    % Aerodynamic matrices construction 
    b = 0.5 * chord;
    eTh = (0.33 - 0.50) * 2;
    rho = 1.225;
    CLa = 2 * pi;
    
    % Integral computation 
    I_wh = zeros(No, No);
    I_tt = zeros(No, No);
    I_wt = zeros(No, No);
    
    for i = 1:No
        for j = 1:No
            I_wh(i,j) = -trapz(x, Modes(i, 1:3:end) .* Modes(j, 1:3:end)) * L;
            I_tt(i,j) =  trapz(x, Modes(i, 3:3:end) .* Modes(j, 3:3:end)) * L;
            I_wt(i,j) =  trapz(x, Modes(i, 1:3:end) .* Modes(j, 3:3:end)) * L;
        end
    end
    I_th = -I_wt';
    
    % MATRIX A_acc (Acceleration)
    Mae = zeros(2, 2);
    Mae(1,1) =  pi;
    Mae(1,2) = -pi * b * eTh;
    Mae(2,1) =  pi * b * eTh;
    Mae(2,2) = -pi * b^2 * (1/8 + eTh^2);
    Mae = Mae * b^2;
    
    % MATRIX A_vel (Velocity)
    Cae = zeros(2, 2);
    Cae(1,2) =  pi;
    Cae(2,2) = -pi * b * (1/2 - eTh);
    Cae = Cae * b^2;
    
    % MATRIX A_pos (Position)
    Kae = zeros(2, 2);
    Kae = Kae * b^2;
    
    % MATRICES Cw, Cw_cap, Bw
    Cw = [0, 1];
    Cw_cap = [1, b * (1/2 - eTh)];
    Bw = b * CLa * [1; b * (eTh + 1/2)];
    
    % Theodorsen function C(k) = H1(k) / (H1(k) + jH0(k))
    Theodorsen = @(k) besselh(1,2,k) ./ (besselh(1,2,k) + 1i*besselh(0,2,k));
    
    % Base matrices for the generalized eigenvalue problem
    A1 = [eye(No), zeros(No);
          zeros(No), eye(No)];
          
    A2 = [zeros(No), eye(No);
         -diag(f(1:No).^2), zeros(No)];
         
    % Range of velocity
    U_range = 0:2:200;
    
    % Preallocation
    eigenvectors = zeros(2*No, 2*No, length(U_range));
    eigenvalues  = zeros(2*No, length(U_range));
    
    % Structural eigenvalues (initialization)
    [eigenvectors(:,:,1), eig_diag] = eig(A2, A1);
    eigenvalues(:,1) = diag(eig_diag);
    
    toll = 1e-6;
    iter_max = 200;
    
    %% Step 7: p-k Iterative Solver
    for i = 2:length(U_range)
        U = U_range(i);
        for j = 1:2*No
            
            n_iter = 0;
            eigenval_0 = 0;
            eigenval = eigenvalues(j, i-1);
            
            while (n_iter < iter_max) && (abs(eigenval - eigenval_0) > toll)
                omega_guess = imag(eigenval);        
                k = max(abs(omega_guess) * b / U, 1e-6);
                eigenval_0 = eigenval;
                Ck = Theodorsen(k);
                
                % Aerodynamic matrices
                K_aer = U^2 * rho * (Kae + Bw * real(Ck) * Cw - k/b * Bw * imag(Ck) * Cw_cap - k^2/b^2 * Mae);
                C_aer = U * rho * (Cae + b/k * Bw * imag(Ck) * Cw + Bw * real(Ck) * Cw_cap);
                
                K_aer_full = I_wh*K_aer(1,1) + I_wt*K_aer(1,2) + I_th*K_aer(2,1) + I_tt*K_aer(2,2);
                C_aer_full = I_wh*C_aer(1,1) + I_wt*C_aer(1,2) + I_th*C_aer(2,1) + I_tt*C_aer(2,2);
                
                % New dynamic matrix
                A2_new = [zeros(No), eye(No);
                         -diag(f(1:No).^2) + K_aer_full, C_aer_full];
         
                % Critical eigenvalue tracking
                [eigenvectors(:,j,i), eig_diag] = eigs(A2_new, 1, eigenval_0);
                eigenvalues(j,i) = diag(eig_diag);
                eigenval = eigenvalues(j,i);
                n_iter = n_iter + 1;
            end
        end
    end
    
    %% Step 8: Find Exact Flutter Speed (Interpolation)
    flutternotfound = true;
    idx_U = NaN;
    idx_mode = NaN;
    
    for j = 1:length(U_range)
        for i = 1:size(eigenvalues, 1) 
            if real(eigenvalues(i, j)) >= 1e-6 && flutternotfound
                flutternotfound = false;
                idx_U = j;
                idx_mode = i;
            end
        end
    end
    
    if ~isnan(idx_U) && ~isnan(idx_mode)
        % Linear interpolation
        m_slope = (real(eigenvalues(idx_mode, idx_U-1)) - real(eigenvalues(idx_mode, idx_U))) / (U_range(idx_U-1) - U_range(idx_U));
        n_intercept = real(eigenvalues(idx_mode, idx_U-1)) - m_slope * U_range(idx_U-1);
        Uflut = -n_intercept / m_slope;
        fprintf('\n>>> FLUTTER SPEED (%d Modes): %.2f m/s <<<\n', No, Uflut);
    else
        fprintf('\n>>> No Flutter detected in range <<<\n');
        Uflut = NaN;
    end
    
    %% Step 9: Plotting (V-g, Frequency, and Sigma)
    
    % Data Preparation
    omega = abs(imag(eigenvalues));
    sigma = real(eigenvalues);
    omega(omega < 1e-6) = 1e-6; % Avoid division by zero
    
    g_damping = 2 * sigma ./ omega;
    freq_hz = omega / (2*pi);
    n_modes_total = size(eigenvalues, 1);
    colors = turbo(n_modes_total);
    
    % --- 1. V-g Diagram ---
    figure(100); clf;
    hold on; grid on; box on;
    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
    
    for m = 1:n_modes_total
        if mean(freq_hz(m,:)) > 0.1 % Filter out rigid/static modes
            plot(U_range, g_damping(m,:), 'LineWidth', 2.0, 'Color', colors(m,:));
        end
    end
    
    if ~isnan(Uflut)
        plot(Uflut, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
    end
    
    ylabel('Damping $g$ [-]', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, max(U_range)]);
    ylim([-0.2, 0.1]);
    
    % --- 2. Frequency vs. Speed Plot ---
    figure(101); clf;
    hold on; grid on; box on;
    for m = 1:n_modes_total
        if mean(freq_hz(m,:)) > 0.1
            plot(U_range, freq_hz(m,:), 'LineWidth', 2.0, 'Color', colors(m,:));
        end
    end
    xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex');
    ylabel('Frequency [Hz]', 'Interpreter', 'latex');
    xlim([0, max(U_range)]);
    ylim([0, 100]);
    
    % --- 3. Sigma vs. Speed Plot ---
    figure(102); clf;
    hold on; grid on; box on;
    for m = 1:n_modes_total
        if mean(freq_hz(m,:)) > 0.1
            plot(U_range, -sigma(m,:) ./ abs(eigenvalues(m,:)), 'LineWidth', 2.0, 'Color', colors(m,:));
        end
    end
    xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex');
    ylabel('$\sigma$', 'Interpreter', 'latex');
    xlim([0, max(U_range)]);
    ylim([-0.3, 0.3]);
    
end