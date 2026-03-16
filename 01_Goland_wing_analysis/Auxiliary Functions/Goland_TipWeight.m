%% =======================================================================
% FUNCTION: Goland_TipWeight
% DESCRIPTION: Computes the flutter and divergence speeds for the Goland 
%              Wing when a concentrated mass (tip weight) is added at the 
%              wing tip. 
%
% AUTHORS: Giuseppe Quaranta <giuseppe.quaranta@polimi.it>
%          Lucrezia Mapelli <lucrezia.mapelli@polimi.it>
%          Aeroservoelasticity of Fixed and Rotary wing Aircraft (2025)
%          Modifications by Alejandro Rivera Míguez
%
% INPUT:
%   st_meters - Distance from Elastic Axis (EA) to the CG of the tip mass [m].
%               Positive: Mass is BEHIND the EA.
%               Negative: Mass is IN FRONT OF the EA.
%
% OUTPUT:
%   Uflut     - Computed flutter speed [m/s]
%% =======================================================================
function Uflut = Goland_TipWeight(st_meters)
    
    %% 1. Wing and Tip Mass Parameters
    mt = 80;               % Tip mass [kg]
    Jpt_cg = 15;           % Inertia of the tip mass at its CG [kg m^2]
    
    % Parallel Axis Theorem: Inertia seen from the Elastic Axis (EA)
    Jpt_EA = Jpt_cg + mt * st_meters^2; 
    
    EJ = 9.77e6;  m_wing = 35.75;
    GJ = 9.876e5; Jp_wing = 8.65;
    L = 6.096;    chord = 1.829;
    
    xt_dist = 0.1 * chord; % Distributed static unbalance of the clean wing
    
    %% 2. Obtain Clean Wing Modes (FEM)
    n_elem = 100; % Sufficient elements for accurate integration
    x = linspace(0, 1, n_elem+1);
    
    EJ_v = EJ * ones(1, n_elem); m_v  = m_wing * ones(1, n_elem);
    GJ_v = GJ * ones(1, n_elem); Jp_v = Jp_wing * ones(1, n_elem);
    
    n_modb = 10; % Number of base bending modes
    n_modt = 10; % Number of base torsion modes
    
    [fb, Modesw] = FEM_beam_bending(x, EJ_v, m_v, L, [1,2], n_modb);
    [ft, Modest] = FEM_beam_torsion(x, GJ_v, Jp_v, L, 1, n_modt);
    
    %% 3. Construct Generalized Matrices WITH Tip Mass
    % Extract modal values at the TIP (last node)
    phi_h_tip = Modesw(:, end);    % [n_modb x 1]
    phi_a_tip = Modest(:, end);    % [n_modt x 1]
    
    % --- STIFFNESS MATRIX (K_gen) ---
    % Tip mass adds no stiffness, only clean wing stiffness is used
    K_gen = diag([fb.^2; ft.^2]);
    
    % --- MASS MATRIX (M_gen) ---
    % 1. Clean Wing Diagonal (Identity due to normalization)
    M_gen = eye(n_modb + n_modt);
    
    % 2. Clean Wing Distributed Coupling
    Mbt_dist = zeros(n_modb, n_modt);
    m_n = 35.75; 
    for i = 1:n_modb
        for j = 1:n_modt
            % Negative sign to match standard formulation (mass behind EA)
            Mbt_dist(i,j) = -trapz(x, Modesw(i,:) .* m_n .* xt_dist .* Modest(j,:)) * L;
        end
    end
    
    % Insert off-diagonal coupling blocks
    M_gen(1:n_modb, n_modb+1:end) = Mbt_dist; 
    M_gen(n_modb+1:end, 1:n_modb) = Mbt_dist';
    
    % 3. ADD TIP MASS EFFECTS
    % A) Bending-Bending (Add mass mt)
    M_gen(1:n_modb, 1:n_modb) = M_gen(1:n_modb, 1:n_modb) + mt * (phi_h_tip * phi_h_tip');
    
    % B) Torsion-Torsion (Add Inertia Jpt_EA)
    idx_t = n_modb+1 : n_modb+n_modt;
    M_gen(idx_t, idx_t) = M_gen(idx_t, idx_t) + Jpt_EA * (phi_a_tip * phi_a_tip');
    
    % C) Bending-Torsion Coupling at Tip (mt * st)
    M_coupling_tip = -mt * st_meters * (phi_h_tip * phi_a_tip');
    M_gen(1:n_modb, idx_t) = M_gen(1:n_modb, idx_t) + M_coupling_tip;
    M_gen(idx_t, 1:n_modb) = M_gen(idx_t, 1:n_modb) + M_coupling_tip';
    
    %% 4. Solve the New Structural Problem (Loaded Wing)
    [V_new, D_new] = eig(K_gen, M_gen);
    [freqs_new, idx_sort] = sort(sqrt(diag(D_new))); 
    V_new = V_new(:, idx_sort);
    
    % Normalize new modes with respect to the new mass matrix (CRITICAL for aeroelasticity)
    m_modal_new = diag(V_new' * M_gen * V_new);
    V_new = V_new ./ sqrt(m_modal_new');
    
    %% 5. Reconstruct Physical Modal Shapes
    No = 6; % Retain first 6 coupled modes
    Modes_aero_h = zeros(No, length(x));
    Modes_aero_a = zeros(No, length(x));
    
    for i = 1:No
        contrib_h = V_new(1:n_modb, i); 
        Modes_aero_h(i,:) = (contrib_h' * Modesw);
        
        contrib_a = V_new(n_modb+1:end, i);
        Modes_aero_a(i,:) = (contrib_a' * Modest);
    end
    f_pk = freqs_new(1:No);
    
    %% 6. Generalized Aerodynamic Matrices
    I_wh = zeros(No, No); 
    I_tt = zeros(No, No); 
    I_wt = zeros(No, No);
    
    for i = 1:No
        for j = 1:No
            I_wh(i,j) = -trapz(x, Modes_aero_h(i,:) .* Modes_aero_h(j,:)) * L; 
            I_tt(i,j) =  trapz(x, Modes_aero_a(i,:) .* Modes_aero_a(j,:)) * L;
            I_wt(i,j) =  trapz(x, Modes_aero_h(i,:) .* Modes_aero_a(j,:)) * L;
        end
    end
    I_th = -I_wt';
    
    %% 7. p-k Solver (Theodorsen)
    b = 0.5 * chord;
    eTh = -0.34; % (0.33 - 0.5)*2
    rho = 1.225;
    CLa = 2 * pi;
    
    % Aerodynamic Strip Matrices
    Mae = b^2 * [ pi,         -pi*b*eTh;
                  pi*b*eTh,   -pi*b^2*(1/8 + eTh^2) ];
              
    Cae = b^2 * [ 0,           pi;
                  0,          -pi*b*(1/2 - eTh) ];
              
    Kae = zeros(2, 2); 
    
    Cw = [0, 1];
    Cw_cap = [1, b * (1/2 - eTh)];
    Bw = b * CLa * [1; b * (eTh + 1/2)];
    
    Theodorsen = @(k) besselh(1,2,k) ./ (besselh(1,2,k) + 1i*besselh(0,2,k));
    
    % Range extended to capture stable/high-speed cases
    U_range = 1:2:300; 
    eigenvalues  = zeros(2*No, length(U_range));
    
    % Initialization (Quasi-steady / structural solution)
    A2_init = [zeros(No), eye(No); -diag(f_pk.^2), zeros(No)];
    [~, eig_diag] = eig(A2_init); 
    eigenvalues(:,1) = diag(eig_diag); 
    
    for i = 2:length(U_range)
        U = U_range(i);
        for j = 1:2*No
            eigenval = eigenvalues(j, i-1);
            
            for iter = 1:50
                k = max(abs(imag(eigenval)) * b / U, 1e-4);
                Ck = Theodorsen(k);
                
                K_aer_2D = U^2 * rho * (Kae + Bw * real(Ck) * Cw - k/b * Bw * imag(Ck) * Cw_cap - k^2/b^2 * Mae);
                C_aer_2D = U * rho * (Cae + b/k * Bw * imag(Ck) * Cw + Bw * real(Ck) * Cw_cap);
                
                K_aer_gen = I_wh*K_aer_2D(1,1) + I_wt*K_aer_2D(1,2) + I_th*K_aer_2D(2,1) + I_tt*K_aer_2D(2,2);
                C_aer_gen = I_wh*C_aer_2D(1,1) + I_wt*C_aer_2D(1,2) + I_th*C_aer_2D(2,1) + I_tt*C_aer_2D(2,2);
                
                A_sys = [zeros(No), eye(No);
                        -diag(f_pk.^2) + K_aer_gen, C_aer_gen];
                
                eigenval_old = eigenval;
                try
                    val_new = eigs(A_sys, 1, eigenval_old);
                catch
                    [vv, dd] = eig(A_sys);
                    dd = diag(dd);
                    [~, idx] = min(abs(dd - eigenval_old));
                    val_new = dd(idx);
                end
                
                eigenval = val_new;
                if abs(eigenval - eigenval_old) < 1e-4
                    break;
                end
            end
            eigenvalues(j,i) = eigenval;
        end
    end
    
    %% 8. Detect Flutter Speed (Zero Damping Crossing)
    flutternotfound = true;
    Uflut = NaN;
    
    for j = 2:length(U_range)
        for m = 1:2*No
            sigma = real(eigenvalues(m,j));
            omega = abs(imag(eigenvalues(m,j)));
            
            % If damping becomes positive and it is an oscillatory mode (>0.1 rad/s)
            if sigma > 0 && omega > 0.1
                if flutternotfound
                    % Linear Interpolation
                    U1 = U_range(j-1); U2 = U_range(j);
                    s1 = real(eigenvalues(m,j-1)); s2 = real(eigenvalues(m,j));
                    
                    Uflut = U1 - s1 * (U2 - U1) / (s2 - s1);
                    flutternotfound = false;
                end
            end
        end
        if ~flutternotfound, break; end
    end
    
    %% 9. Divergence Speed Calculation
    % Extract the torsion block from the generalized stiffness matrix
    n_mod = n_modt; 
    idx_torsion = (n_mod + 1) : (2 * n_mod); 
    K_tors = K_gen(idx_torsion, idx_torsion);
    
    dist_AC_EA = b * (eTh - (-0.5)); % Lever arm (AC to EA)
    term_div = chord * CLa * dist_AC_EA;
    
    Q_alpha = zeros(n_mod, n_mod);
    for i = 1:n_mod
        for j = 1:n_mod
            Q_alpha(i,j) = trapz(x, Modest(i,:) .* Modest(j,:)) * term_div * L;
        end
    end
    
    % Solve Generalized Eigenvalue Problem (Torsion only)
    [~, D_div] = eig(K_tors, Q_alpha);
    q_vals = diag(D_div);
    
    % Find lowest positive real value
    q_divs = q_vals(q_vals > 0 & imag(q_vals) == 0);
    
    if ~isempty(q_divs)
        q_crit = min(q_divs);
        U_divergence = sqrt(2 * q_crit / rho);
        % fprintf('\n>>> DIVERGENCE SPEED: %.2f m/s <<<\n', U_divergence);
    else
        % fprintf('\n>>> No divergence found. <<<\n');
        U_divergence = NaN;
    end
    
    %% 10. Optional V-g / V-f Plotting (Uncomment to use)
    % figure(100); clf;
    % set(gcf, 'Color', 'w', 'Position', [100 100 1000 600]);
    % 
    % subplot(2,1,1); hold on; grid on; box on;
    % yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
    % 
    % g_mat = zeros(size(eigenvalues));
    % f_mat = zeros(size(eigenvalues));
    % 
    % for m = 1:size(eigenvalues,1)
    %     omega = abs(imag(eigenvalues(m,:)));
    %     g_mat(m,:) = 2 * real(eigenvalues(m,:)) ./ (omega + eps);
    %     f_mat(m,:) = omega / (2*pi); 
    % 
    %     if mean(f_mat(m,:)) > 0.5
    %         plot(U_range, g_mat(m,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Mode %d', m));
    %     end
    % end
    % 
    % title('\textbf{V-g Diagram (Damping)}', 'Interpreter', 'latex', 'FontSize', 14);
    % ylabel('Damping $g$ [-]', 'Interpreter', 'latex');
    % xlim([0, max(U_range)]);
    % ylim([-0.3, 0.3]); 
    % legend('show', 'Location', 'eastoutside');
    % 
    % subplot(2,1,2); hold on; grid on; box on;
    % for m = 1:size(eigenvalues,1)
    %     if mean(f_mat(m,:)) > 0.5
    %         plot(U_range, f_mat(m,:), 'LineWidth', 1.5);
    %     end
    % end
    % 
    % title('\textbf{V-f Diagram (Frequency)}', 'Interpreter', 'latex', 'FontSize', 14);
    % xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex');
    % ylabel('Frequency [Hz]', 'Interpreter', 'latex');
    % xlim([0, max(U_range)]);
end