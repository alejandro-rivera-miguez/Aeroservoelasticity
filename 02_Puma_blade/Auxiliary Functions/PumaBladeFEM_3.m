%% =======================================================================
% FUNCTION: PumaBladeFEM_3_V5
% DESCRIPTION: Computes the aeroelastic stability of the Puma helicopter
%              blade in hover with variable pitch link stiffness and CG
%              offset, using a full basis (6 modes).
%
% AUTHORS: Alejandro Rivera Míguez
%          (Based on course material by G. Quaranta)
% =======================================================================
function [stability, damping_box] = PumaBladeFEM_3(link_stiff, cg_off)

    %% 1. DATA AND INITIALIZATION
    Data = puma_data(); 
    Omega = Data.omega;
    L = Data.blade_radius;
    c = Data.blade_chord;
    
    % Mesh
    n_elem_f = 6; 
    n_elem_b = 80;
    x_grid = (0:1/n_elem_f:1)*(Data.pitch_bearing - Data.flap_hinge_axis) + Data.flap_hinge_axis;
    xb = 0.610/(Data.blade_radius - Data.pitch_bearing);
    x_grid = [x_grid, (xb:(1-xb)/n_elem_b:1)*(Data.blade_radius - Data.pitch_bearing) + Data.pitch_bearing];
    x_norm = x_grid / L;
    x_col = x_grid(:);
    
    %% 2. STRUCTURAL RESOLUTION
    % Interpolation
    [EJ_el, ~] = w_interpolate(Data.blade_flap_chord_stiffness(:,1), Data.blade_flap_chord_stiffness(:,2), x_grid); 
    [m_el, ~]  = w_interpolate(Data.blade_mass(:,1), Data.blade_mass(:,2), x_grid); 
    [GJ_el, ~] = w_interpolate(Data.blade_torsional_stiffness(:,1), Data.blade_torsional_stiffness(:,2), x_grid);
    [J_el, ~]  = w_interpolate(Data.blade_torsional_inertia(:,1), Data.blade_torsional_inertia(:,2), x_grid);
    [xi_el, ~] = w_interpolate(Data.blade_cg_offset(:,1), Data.blade_cg_offset(:,2), x_grid); 
    
    offset_meters = (cg_off - 25)/100 * c;
    xi_el = ones(size(xi_el)) * offset_meters; 
    
    % Root Cutout (Ritz Match)
    GJ_el(1 : 2*n_elem_f) = 0; 
    J_el(1 : 2*n_elem_f)  = 0;
    
    % --- FEM MODES ---
    n_modb = 10;
    [fb, Modesw, Modeswp] = FEM_beam_bending(x_norm, EJ_el(1:2:end), m_el(1:2:end), L, 1, [], n_modb);
    node_bearing = n_elem_f + 1;
    cdofs_t = 1:(node_bearing - 1); 
    gdofs_t = [node_bearing, link_stiff]; 
    n_modt = 10;
    [ft, Modest] = FEM_beam_torsion(x_norm, GJ_el(1:2:end), J_el(1:2:end), L, cdofs_t, gdofs_t, n_modt);
    
    psi   = Modesw;
    psi_p = Modeswp;
    phi   = Modest;     
    
    % Coupled problem
    n_mod = n_modb + n_modt;
    K = zeros(n_mod, n_mod);
    Kg = zeros(n_mod, n_mod);
    M = eye(n_mod, n_mod);
    K(1:n_modb,1:n_modb) = diag(fb.^2);
    K(n_modb+1:n_mod,n_modb+1:n_mod) = diag(ft.^2);
    
    % Coupled mass matrix
    Mc = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            Mc(i,j) = trapz(x_grid, psi(i,:) .* m_el([1,2:2:end])' .* xi_el([1,2:2:end])' .* phi(j,:));
        end
    end
    M(1:n_modb,n_modb+1:n_mod) = M(1:n_modb,n_modb+1:n_mod) - Mc;
    M(n_modb+1:n_mod,1:n_modb) = M(n_modb+1:n_mod,1:n_modb) - Mc';
    
    % Centrifugal bending stiffness
    N = zeros(size(x_grid));
    N(end-1) = m_el(end)/2*(x_grid(end)^2 - x_grid(end-1)^2);
    for i = length(N)-2:-1:1
        N(i) = N(i+1) + m_el(i*2)/2*(x_grid(i+1)^2 - x_grid(i)^2);
    end
    for i = 1:n_modb
        for j = 1:n_modb
            Kg(i,j) = trapz(x_grid, psi_p(i, :) .* N .* psi_p(j,:));
        end
    end
    
    % Centrifugal torsional stiffness
    for i = 1:n_modt
        for j = 1:n_modt
            Kg(n_modb+i,n_modb+j) = trapz(x_grid, phi(i, :) .* J_el([1, 2:2:end])' .* phi(j,:));
        end
    end
    
    Kgc = zeros(n_modb, n_modt);
    
    % Centrifugal coupling bending and torsion
    for i = 1:n_modb
        for j = 1:n_modt
            Kgc(i,j) = trapz(x_grid, Modeswp(i, :) .* x_grid .* m_el([1,2:2:end])' .* xi_el([1,2:2:end])' .* Modest(j,:));
        end
    end
    Kg(1:n_modb,n_modb+1:n_mod) = Kg(1:n_modb,n_modb+1:n_mod) - Kgc;
    Kg(n_modb+1:n_mod,1:n_modb) = Kg(n_modb+1:n_mod,1:n_modb) - Kgc';
    Kt = K + Omega^2*Kg;
    
    [V_raw, E] = eig(-Kt, M); % "Raw" eigenvectors
    [Eo, I_sort] = sort(diag(E), 'descend'); 
    
    idx = [1,2];
    f = imag(sqrt(Eo(idx)));
    No = length(idx);
    
    % Select only the modes of interest
    stability = false;
    if f > 1e-6
        V_sel = V_raw(:, I_sort(idx)); 
        
        %% --- CRITICAL STEP: MASS NORMALIZATION ---
        for i = 1:No
            vec = V_sel(:, i);
            modal_mass_i = vec.' * M * vec; 
            scale_factor = sqrt(modal_mass_i);
            V_sel(:, i) = vec / scale_factor;
        end
        
        %% --- PROJECTION TO PHYSICAL MODES ---
        q_b = V_sel(1:n_modb, :); 
        Modesgw = q_b.' * psi;
        Modesgwp = q_b.' * psi_p;
        
        q_t = V_sel(n_modb+1:end, :);
        Modesgt = q_t.' * phi;
        
        % --- FINAL MATRIX DEFINITION FOR AEROELASTICITY ---
        M = eye(No);                 
        K = diag(f.^2);              
        
        psi_c = Modesgw';
        psip_c = Modesgwp';
        phi_c = Modesgt';
        
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
        
        for i=1:No
            for j=1:No
                WW_K(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j).*Uvec.^2);
                WW_M(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j));
                WT_K(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j).*Uvec.^2);
                WT_M(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j));
                TW_K(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j).*Uvec.^2);
                TW_M(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j));
                TT_K(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j).*Uvec.^2); 
                TT_M(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j)); 
                WW_C(i, j) = trapz(x_aero, psi_a(:, i).*psi_a(:, j).*Uvec);
                WT_C_Theta(i, j) = trapz(x_aero, psi_a(:, i).*phi_a(:, j).*Uvec);
                WT_C_Omega(i, j) =  trapz(x_aero, Omega*Uvec.*psi_a(:, i).*psi_p_a(:, j)); 
                TW_C(i, j) = trapz(x_aero, phi_a(:, i).*psi_a(:, j).*Uvec); 
                TT_C_Theta(i, j) = trapz(x_aero, phi_a(:, i).*phi_a(:, j).*Uvec);
                TT_C_Omega(i, j) = trapz(x_aero, Omega*Uvec.*phi_a(:, i).*psi_p_a(:, j)); 
            end
        end
        
        %% Resolution Definitions
        tolerancia = 10e-6;
        Amatrix = zeros(2*No,  2*No);
        Amatrix (1:No, (No+1):end) = eye(No);
        tol = 100;
        counter1 = 0;
        QuasiSteady = false;
        
        %% Resolution
        [Autovectores, E_vac] = eig(K, M);
        Autovalores = sqrt(diag(E_vac)); 
        [Autovalores, idx_sort] = sort(Autovalores,'descend');
        Autovectores = Autovectores(:, idx_sort);
        
        for modes=1:No
            k = 4/3*Autovalores(modes)*b/Omega/L;
            while tol > tolerancia
                if QuasiSteady
                    Ck = 1;
                else
                    Ck = besselh(1, 2, k)/(besselh(1, 2, k) + 1i*besselh(0, 2, k));
                end
    
                Ma = rho*b^2*[pi,-pi*b*ep;
                              ep*pi*b,-pi*b^2*(1/8+ep^2)];
                
                Ca = rho*b^2*[0,pi;
                              0,-pi*b*(1/2-ep)];
                              
                Ka = rho*b^2*[0,0;
                              0,0];
                
                Bw = rho*b*Cla*[1;b*(ep+1/2)];
                Cw = [0, 1];
                Chatw = [1,b*(1/2-ep)];
                
                Kaero = (Ka + Bw*Ck*Cw);
                Caero = (Ca + Bw*Ck*Chatw);
                Maero = Ma;
            
                KK = -Kaero(1,1)*WW_K + Caero(1, 2)*WT_C_Omega + Kaero(1,2)*WT_K -Kaero(2,1)*TW_K + Caero(2,2)*TT_C_Omega+ Kaero(2,2)*TT_K;
                CC = -Caero(1,1)*WW_C + Caero(1,2)*WT_C_Theta -Caero(2,1)*TW_C + Caero(2,2)*TT_C_Theta;
                MM = -Maero(1,1)*WW_M + Maero(1,2)*WT_M -Maero(2,1)*TW_M + Maero(2,2)*TT_M;
    
                TM2 = M - MM;
                TM1 = -CC;
                TM0 = K - KK;
                
                Amatrix((No+1):end, 1:(No)) = -TM2\TM0;
                Amatrix((No+1):end, (No+1):end) = -TM2\TM1;
                
                [AVEC, z] = eig(Amatrix);
    
                if counter1 == 0
                    Autovectores_track = AVEC((No+1):end, :);
                    counter1 = 1;
                else 
                    Autovectores_track = AVEC;
                end
                
                % Tracking of the eigenvector
                ScalarProducts = zeros(1, length(Autovectores_track));
                for tracking=1:length(Autovectores_track)
                   ScalarProducts(tracking) = ps_complex(Autovectores_track(:, tracking), Autovectores(:, modes));
                end
    
                [aux1, Index2] = max(ScalarProducts); 
    
                z = diag(z);
                z = z(Index2); % Eigenvalue of the mode we are studying
    
                diff1 = abs(Autovalores(modes) - sqrt(imag(z)^2 + real(z)^2)); 
    
                Autovalores(modes) = sqrt(imag(z)^2 + real(z)^2);
                k = 4/3*Autovalores(modes)*b/Omega/L;
                
                % ORIGINAL EXACT DAMPING EQUATION PRESERVED
                damping = real(z)/Autovalores(modes);
                
                tol = diff1;
                
                % Redefinition of eigenvectors
                Autovectores(1:(No), modes) = AVEC(1:(No), Index2);
                Autovectores((No+1):(2*No), modes) = AVEC((No+1):end, Index2);
            end
            tol = 100;
            k_box(modes) = k;
            damping_box(modes) = damping;
            freq_box(modes) = Autovalores(modes);
            counter1 = 0;
        end
        
        stability = true;
        for kkk = 1:No
            if damping_box(kkk) > 1e-3 && stability 
                stability = false;
            end
        end
        
    else
        damping_box = zeros(1,No);
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