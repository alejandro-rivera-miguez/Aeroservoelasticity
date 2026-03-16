%% =======================================================================
% FUNCTION: PumaBladeFEM_5
% DESCRIPTION: Computes the aeroelastic stability of the Puma helicopter
%              blade in hover using an EXACT UNSTEADY formulation.
%              Unlike previous tasks, the Theodorsen function C(k) is 
%              evaluated at every spanwise station and kept inside the 
%              aerodynamic integrals. Uses a full basis (6 modes).
%
% AUTHORS: Alejandro Rivera Míguez
%          (Based on course material by G. Quaranta)
%
% INPUTS:
%   link_stiff - Pitch link stiffness [N/m or Nm/rad]
%   cg_off     - CG offset in % of chord
% OUTPUTS:
%   stability  - Boolean indicating if the system is stable
%   damping_box- Array containing damping values for the coupled modes
%% =======================================================================
function [stability, damping_box] = PumaBladeFEM_5(link_stiff, cg_off)

    %% 1. DATA AND INITIALIZATION
    Data = puma_data(); 
    Omega = Data.omega;
    L = Data.blade_radius;
    c = Data.blade_chord;
    
    % Spatial Mesh Generation
    n_elem_f = 6; 
    n_elem_b = 80;
    x_grid = (0:1/n_elem_f:1)*(Data.pitch_bearing - Data.flap_hinge_axis) + Data.flap_hinge_axis;
    xb = 0.610 / (Data.blade_radius - Data.pitch_bearing);
    x_grid = [x_grid, (xb:(1-xb)/n_elem_b:1)*(Data.blade_radius - Data.pitch_bearing) + Data.pitch_bearing];
    x_norm = x_grid / L;
    x_col = x_grid(:);
    
    %% 2. STRUCTURAL RESOLUTION
    % Parameter Interpolation over the blade span
    [EJ_el, ~] = w_interpolate(Data.blade_flap_chord_stiffness(:,1), Data.blade_flap_chord_stiffness(:,2), x_grid); 
    [m_el, ~]  = w_interpolate(Data.blade_mass(:,1), Data.blade_mass(:,2), x_grid); 
    [GJ_el, ~] = w_interpolate(Data.blade_torsional_stiffness(:,1), Data.blade_torsional_stiffness(:,2), x_grid);
    [J_el, ~]  = w_interpolate(Data.blade_torsional_inertia(:,1), Data.blade_torsional_inertia(:,2), x_grid);
    [xi_el, ~] = w_interpolate(Data.blade_cg_offset(:,1), Data.blade_cg_offset(:,2), x_grid); 
    
    % Shift CG offset based on user input
    offset_meters = (cg_off - 25)/100 * c;
    xi_el = ones(size(xi_el)) * offset_meters; 
    
    % Root Cutout Modification (Ritz Matching)
    GJ_el(1 : 2*n_elem_f) = 0; 
    J_el(1 : 2*n_elem_f)  = 0;
    
    % --- FEM MODES COMPUTATION ---
    n_modb = 10;
    [fb, Modesw, Modeswp] = FEM_beam_bending(x_norm, EJ_el(1:2:end), m_el(1:2:end), L, 1, [], n_modb);
    
    node_bearing = n_elem_f + 1;
    cdofs_t = 1:(node_bearing - 1); 
    gdofs_t = [node_bearing, link_stiff]; 
    n_modt = 10;
    [ft, Modest] = FEM_beam_torsion(x_norm, GJ_el(1:2:end), J_el(1:2:end), L, cdofs_t, gdofs_t, n_modt);
    
    psi    = Modesw;
    psi_p  = Modeswp;
    phi    = Modest;     
    
    % --- COUPLED BENDING-TORSION PROBLEM ASSEMBLY ---
    n_mod = n_modb + n_modt;
    K = zeros(n_mod, n_mod);
    Kg = zeros(n_mod, n_mod);
    M = eye(n_mod, n_mod);
    K(1:n_modb, 1:n_modb) = diag(fb.^2);
    K(n_modb+1:n_mod, n_modb+1:n_mod) = diag(ft.^2);
    
    % Coupled Mass Matrix
    Mc = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            Mc(i,j) = trapz(x_grid, psi(i,:) .* m_el([1,2:2:end])' .* xi_el([1,2:2:end])' .* phi(j,:));
        end
    end
    M(1:n_modb, n_modb+1:n_mod) = M(1:n_modb, n_modb+1:n_mod) - Mc;
    M(n_modb+1:n_mod, 1:n_modb) = M(n_modb+1:n_mod, 1:n_modb) - Mc';
    
    % Centrifugal Bending Stiffness (Geometric Stiffness)
    N = zeros(size(x_grid));
    N(end-1) = m_el(end)/2 * (x_grid(end)^2 - x_grid(end-1)^2);
    for i = length(N)-2:-1:1
        N(i) = N(i+1) + m_el(i*2)/2 * (x_grid(i+1)^2 - x_grid(i)^2);
    end
    for i = 1:n_modb
        for j = 1:n_modb
            Kg(i,j) = trapz(x_grid, psi_p(i, :) .* N .* psi_p(j,:));
        end
    end
    
    % Centrifugal Torsional Stiffness
    for i = 1:n_modt
        for j = 1:n_modt
            Kg(n_modb+i, n_modb+j) = trapz(x_grid, phi(i, :) .* J_el([1, 2:2:end])' .* phi(j,:));
        end
    end
    
    % Centrifugal Coupling (Bending - Torsion)
    Kgc = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            Kgc(i,j) = trapz(x_grid, Modeswp(i, :) .* x_grid .* m_el([1,2:2:end])' .* xi_el([1,2:2:end])' .* Modest(j,:));
        end
    end
    Kg(1:n_modb, n_modb+1:n_mod) = Kg(1:n_modb, n_modb+1:n_mod) - Kgc;
    Kg(n_modb+1:n_mod, 1:n_modb) = Kg(n_modb+1:n_mod, 1:n_modb) - Kgc';
    
    % Total Stiffness Matrix
    Kt = K + Omega^2 * Kg;
    
    % Solve initial Eigenvalue Problem
    % "Raw" eigenvectors (MATLAB normalizes them to Euclidean norm 1, not to mass)
    [V_raw, E_val] = eig(-Kt, M); 
    [Eo, I_sort] = sort(diag(E_val), 'descend'); 
    
    % Select the modes of interest (Full basis: 6 modes)
    idx = [1, 2, 3, 4, 5, 6];
    f = imag(sqrt(Eo(idx)));
    No = length(idx);
    
    stability = false;
    
    if f > 1e-6
        V_sel = V_raw(:, I_sort(idx)); 
        
        %% --- CRITICAL STEP: MASS NORMALIZATION ---
        % We enforce V_sel(:,i)' * M * V_sel(:,i) = 1
        for i = 1:No
            vec = V_sel(:, i);
            
            % 1. Calculate the Modal Mass for this specific mode
            modal_mass_i = vec.' * M * vec; 
            
            % 2. Calculate scale factor
            scale_factor = sqrt(modal_mass_i);
            
            % 3. Normalize the eigenvector
            V_sel(:, i) = vec / scale_factor;
        end
        
        %% --- PROJECTION TO PHYSICAL MODES ---
        % Project normalized vectors onto bending and torsion shapes
        
        % Bending Part
        q_b = V_sel(1:n_modb, :); 
        Modesgw  = q_b.' * psi;
        Modesgwp = q_b.' * psi_p;
        
        % Torsion Part
        q_t = V_sel(n_modb+1:end, :);
        Modesgt = q_t.' * phi;
        
        %% --- MASS ORTHOGONALITY VERIFICATION ---
        % Verify if the generalized mass matrix is actually the identity
        M_gen_check = V_sel.' * M * V_sel; 
        
        % fprintf('\n=== MASS ORTHOGONALITY VERIFICATION ===\n');
        % disp('Generalized Mass Matrix (Should be Identity):');
        % disp(M_gen_check);
        % error_norm = norm(M_gen_check - eye(No));
        % if error_norm < 1e-10
        %     fprintf('SUCCESS: Matrix is Unitary. Error: %e\n', error_norm);
        % else
        %     fprintf('WARNING: Matrix is NOT Unitary. Error: %e\n', error_norm);
        % end
        
        % --- FINAL MATRIX DEFINITION FOR AEROELASTICITY ---
        M = eye(No);                 
        K = diag(f.^2);              
        
        psi_c  = Modesgw';
        psip_c = Modesgwp';
        phi_c  = Modesgt';
        
        %% 3. AERODYNAMICS
        ep = -0.5;
        b = c/2;
        rho = 1.225;
        Cla = 2*pi;
        
        %% Aerodynamic integral shape functions
        mask = x_grid >= Data.blade_cutout;
        x_aero = x_col(mask);
        
        psi_a   = psi_c(mask,:); 
        phi_a   = phi_c(mask,:);
        psi_p_a = psip_c(mask,:);
        Uvec = Omega * x_aero;
        
        % Preallocate matrices for numerical integration
        WW_K = zeros(No); WW_M = zeros(No); WT_K = zeros(No); WT_M = zeros(No);
        TW_K = zeros(No); TW_M = zeros(No); TT_K = zeros(No); TT_M = zeros(No);
        WW_C = zeros(No); WT_C_Theta = zeros(No); WT_C_Omega = zeros(No);
        TW_C = zeros(No); TT_C_Theta = zeros(No); TT_C_Omega = zeros(No);
        WW_Ck = zeros(No); WT_Ck_Theta = zeros(No); WT_Ck_Omega = zeros(No);
        TW_Ck = zeros(No); TT_Ck_Theta = zeros(No); TT_Ck_Omega = zeros(No);
        
        %% 4. AEROELASTIC SOLVER DEFINITIONS (p-k METHOD)
        tolerancia = 10e-6;
        Amatrix = zeros(2*No,  2*No);
        Amatrix (1:No, (No+1):end) = eye(No);
        tol = 100;
        counter1 = 0;
        QuasiSteady = false;
        
        %% 5. ITERATIVE RESOLUTION
        [eigenvectors, E_vac] = eig(K, M);
        eigenvalues = sqrt(diag(E_vac)); 
        [eigenvalues, idx_sort] = sort(eigenvalues,'descend');
        eigenvectors = eigenvectors(:, idx_sort);
        
        k_box = zeros(1, No);
        damping_box = zeros(1, No);
        freq_box = zeros(1, No);
        
        for modes = 1:No
            while tol > tolerancia
                % EXACT UNSTEADY: Reduced frequency is a vector along the span
                k = 4/3 * eigenvalues(modes) * b / Omega / L ./ x_aero;
                
                if QuasiSteady
                    Ck = 1;
                else
                    % Ck evaluated at every spanwise station
                    Ck = besselh(1, 2, k) ./ (besselh(1, 2, k) + 1i*besselh(0, 2, k));
                end
    
                Ma = rho*b^2 * [pi, -pi*b*ep; ep*pi*b, -pi*b^2*(1/8+ep^2)];
                Ca = rho*b^2 * [0, pi; 0, -pi*b*(1/2-ep)];
                Ka = rho*b^2 * [0, 0; 0, 0];
                Bw = rho*b*Cla * [1; b*(ep+1/2)];
                Cw = [0, 1];
                Chatw = [1, b*(1/2-ep)];
                
                Kaero  = Bw * Cw;
                Caero1 = Ca;
                Caero2 = Bw * Chatw;
                Maero  = Ma;
                
                % Spanwise integration including C(k(x)) inside the integral
                for i = 1:No
                    for j = 1:No
                        WW_K(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j).*Uvec.^2.*Ck);
                        WW_M(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j));
                        WT_K(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j).*Uvec.^2.*Ck);
                        WT_M(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j));
                        TW_K(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j).*Uvec.^2.*Ck);
                        TW_M(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j));
                        TT_K(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j).*Uvec.^2.*Ck); 
                        TT_M(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j)); 
                        
                        WW_C(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j).*Uvec);
                        WW_Ck(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j).*Uvec.*Ck);
                        WT_C_Theta(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j).*Uvec);
                        WT_Ck_Theta(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j).*Uvec.*Ck);
                        WT_C_Omega(i, j) = trapz(x_aero, Omega*Uvec.*psi_a(:, i).*psi_p_a(:, j)); 
                        WT_Ck_Omega(i, j) = trapz(x_aero, Omega*Uvec.*psi_a(:, i).*psi_p_a(:, j).*Ck); 
                        
                        TW_C(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j).*Uvec); 
                        TW_Ck(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j).*Uvec.*Ck); 
                        TT_C_Theta(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j).*Uvec);
                        TT_Ck_Theta(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j).*Uvec.*Ck);
                        TT_C_Omega(i, j) = trapz(x_aero, Omega*Uvec.*phi_a(:, i).*psi_p_a(:, j)); 
                        TT_Ck_Omega(i, j) = trapz(x_aero, Omega*Uvec.*phi_a(:, i).*psi_p_a(:, j).*Ck); 
                    end
                end
                
                % Generalized Aero Matrices Assembly
                KK = -Kaero(1,1)*WW_K + Caero1(1, 2)*WT_C_Omega + Caero2(1, 2)*WT_Ck_Omega + Kaero(1,2)*WT_K -Kaero(2,1)*TW_K + Caero1(2,2)*TT_C_Omega + Caero2(2,2)*TT_Ck_Omega+ Kaero(2,2)*TT_K;
                CC = -Caero1(1,1)*WW_C - Caero2(1,1)*WW_Ck + Caero1(1,2)*WT_C_Theta + Caero2(1, 2)*WT_Ck_Theta -Caero1(2,1)*TW_C - Caero2(2, 1)*TW_Ck + Caero1(2,2)*TT_C_Theta + Caero2(2,2)*TT_Ck_Theta;
                MM = -Maero(1,1)*WW_M+ Maero(1,2)*WT_M -Maero(2,1)*TW_M + Maero(2,2)*TT_M;
        
                TM2 = M - MM;
                TM1 = -CC;
                TM0 = K - KK;
                
                Amatrix((No+1):end, 1:(No)) = -TM2\TM0;
                Amatrix((No+1):end, (No+1):end) = -TM2\TM1;
                
                [AVEC, z_eig] = eig(Amatrix);
        
                if counter1 == 0
                    eigenvectors_track = AVEC((No+1):end, :);
                    counter1 = 1;
                else 
                    eigenvectors_track = AVEC;
                end
                
               % Modal Tracking of the eigenvector via maximum complex scalar product
               ScalarProducts = zeros(1, length(eigenvectors_track));
               for tracking = 1:length(eigenvectors_track)
                   ScalarProducts(tracking) = ps_complex(eigenvectors_track(:, tracking), eigenvectors(:, modes));
               end
        
               [~, Index2] = max(ScalarProducts); 
        
               z_val = diag(z_eig);
               z_val = z_val(Index2); % Eigenvalue of the mode we are studying
        
               % Check for convergence of omega between iterations
               diff1 = abs(eigenvalues(modes) - sqrt(imag(z_val)^2 + real(z_val)^2)); 
        
               eigenvalues(modes) = sqrt(imag(z_val)^2 + real(z_val)^2);
               k = 4/3 * eigenvalues(modes) * b / Omega / L;
               damping = real(z_val) / eigenvalues(modes);
               tol = diff1;
               
               % Update eigenvectors for next iteration of convergence
               eigenvectors(1:(No), modes) = AVEC(1:(No), Index2);
               eigenvectors((No+1):(2*No), modes) = AVEC((No+1):end, Index2);
            end
            
            tol = 100;
            k_box(modes) = k;
            damping_box(modes) = damping;
            freq_box(modes) = eigenvalues(modes);
            counter1 = 0;
    end
    
    stability = true;
    for kkk = 1:No
        if damping_box(kkk) > 1e-3 && stability 
            stability = false;
        end
    end
    
    %% 6. RESULTS OUTPUT
    % fprintf('\n=========================================\n');
    % fprintf('     TASK 5 RESULTS (EXACT UNSTEADY)\n');
    % fprintf('=========================================\n');
    % fprintf('Mode 1 (Flap):    Freq = %.2f rad/s | Zeta = %.4f\n', freq_box(2), damping_box(2));
    % fprintf('Mode 2 (Torsion): Freq = %.2f rad/s | Zeta = %.4f\n', freq_box(1), damping_box(1));
    % fprintf('-----------------------------------------\n');
else
    damping_box = zeros(1, No);
end

end

%% =======================================================================
% HELPER FUNCTION: ps_complex
% DESCRIPTION: Computes a normalized complex scalar product robustly.
% =======================================================================
function Res = ps_complex(X, Y)
    S1 = 0; S2 = 0; S3 = 0; S4 = 0;
    for jjj = 1:length(X)
        S1 = S1 + real(X(jjj))*real(Y(jjj)) + imag(X(jjj))*imag(Y(jjj));
        S2 = S2 + real(X(jjj))*imag(Y(jjj)) - imag(X(jjj))*real(Y(jjj));
        S3 = norm(X(jjj))^2 + S3;
        S4 = norm(Y(jjj))^2 + S4;
    end
    Res = sqrt(S1^2 + S2^2) / sqrt(S3*S4);
end