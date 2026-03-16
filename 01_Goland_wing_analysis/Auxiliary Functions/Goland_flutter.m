%% =======================================================================
% FUNCTION: Goland_flutter
% DESCRIPTION: Computes the flutter speed for the Goland Wing given a 
%              vector of modified physical parameters.
%
% AUTHORS: Giuseppe Quaranta <giuseppe.quaranta@polimi.it>
%          Lucrezia Mapelli <lucrezia.mapelli@polimi.it>
%          Aeroservoelasticity of Fixed and Rotary wing Aircraft (2025)
%
% INPUT:
%   vec_aux - Vector containing [EJ, GJ, xt, Cla]
% OUTPUT:
%   Uflut   - Computed flutter speed [m/s]
%% =======================================================================
function Uflut = Goland_flutter(vec_aux)
    
    EJ_val  = vec_aux(1);
    GJ_val  = vec_aux(2);
    xt_val  = vec_aux(3);
    Cla_val = vec_aux(4);
    
    No = 6; % Number of coupled modes to retain
    n_elem = 120; % Number of elements for discretization
    
    % Uniform node grid
    x = 0:1/n_elem:1;
    x_av = 1/2 * (x(1:end-1) + x(2:end));
    
    %% Step 1: Compute Bending Modes
    EJ = EJ_val * x_av.^0;
    m  = 35.75 * x_av.^0;
    L  = 6.096;
    cdofs = [1, 2];
    n_modb = 20; % Number of bending modes to compute
    
    [fb, Modesw, Modeswp] = FEM_beam_bending(x, EJ, m, L, cdofs, n_modb);
    
    %% Step 2: Compute Torsion Modes
    GJ = GJ_val * x_av.^0;
    Jp = 8.65 * x_av.^0;
    cdofs_t = 1;
    n_modt = 20; % Number of torsion modes to compute
    
    [ft, Modest] = FEM_beam_torsion(x, GJ, Jp, L, cdofs_t, n_modt);
    
    %% Step 3: Assembly of Cross-Mass Terms (Inertial Coupling)
    chord = 1.829;
    xt_n = xt_val * x.^0;
    m_n  = 35.75 * x.^0;
    
    Mbt = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            Mbt(i,j) = trapz(x, Modesw(i, :) .* m_n .* xt_n .* Modest(j,:)) * L;
        end
    end
    
    %% Step 4: Assembly of Global Structural Matrices
    n_mod = n_modb + n_modt;
    K = zeros(n_mod, n_mod);
    M = eye(n_mod, n_mod);
    
    K(1:n_modb, 1:n_modb) = diag(fb.^2);
    K(n_modb+1:n_mod, n_modb+1:n_mod) = diag(ft.^2);
    
    M(1:n_modb, n_modb+1:end) = -Mbt;
    M(n_modb+1:end, 1:n_modb) = -Mbt';
    
    %% Step 5: Coupled Eigenvalues and Eigenvectors
    [V, E] = eig(-K, M);
    [Eo, I] = sort(diag(E), 'descend');
    f = imag(sqrt(Eo));
    V = V(:, I);
    
    Modes = zeros(No, (n_elem+1)*3);
    Modm = V' * M * V; % Modal mass
    V = V * diag(1 ./ sqrt(diag(Modm))); % Normalize
    
    Modes(:, 1:3:end) = (V(1:n_modb, 1:No)') * Modesw(1:n_modb, :);
    Modes(:, 2:3:end) = (V(1:n_modb, 1:No)') * Modeswp(1:n_modb, :);
    Modes(:, 3:3:end) = (V(n_modb+1:n_mod, 1:No)') * Modest(1:n_modb, :);
    
    %% Step 6: Aeroelastic Model Setup
    b = 0.5 * chord;
    eTh = (0.33 - 0.50) * 2;
    rho = 1.225;
    CLa = Cla_val;
    
    % Aerodynamic Integrals
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
    
    % Strip Theory Matrices (Acceleration, Velocity, Position)
    Mae = b^2 * [ pi,         -pi*b*eTh;
                  pi*b*eTh,   -pi*b^2*(1/8 + eTh^2) ];
              
    Cae = b^2 * [ 0,          pi;
                  0,         -pi*b*(1/2 - eTh) ];
              
    Kae = zeros(2, 2);
    
    Cw = [0, 1];
    Cw_cap = [1, b*(1/2 - eTh)];
    Bw = b * CLa * [1; b*(eTh + 1/2)];
    
    Theodorsen = @(k) besselh(1,2,k) ./ (besselh(1,2,k) + 1i*besselh(0,2,k));
    
    % State-Space Base Matrices
    A1 = [eye(No), zeros(No); zeros(No), eye(No)];
    A2 = [zeros(No), eye(No); -diag(f(1:No).^2), zeros(No)];
    
    %% Step 7: p-k Method Solver
    U_range = 0:2:300;
    eigenvectors = zeros(2*No, 2*No, length(U_range));
    eigenvalues  = zeros(2*No, length(U_range));
    
    [eigenvectors(:,:,1), eig_diag] = eig(A2, A1);
    eigenvalues(:,1) = diag(eig_diag);
    
    toll = 1e-6;
    iter_max = 200;
    
    for i = 2:length(U_range)
        U = U_range(i);
        for j = 1:2*No
            n_iter = 0;
            eigenval_0 = 0;
            eigenval = eigenvalues(j, i-1);
            
            while (n_iter < iter_max) && (abs(eigenval - eigenval_0) > toll)
                omega_guess = imag(eigenval);        
                k = max(abs(omega_guess)*b/U, 1e-6);
                eigenval_0 = eigenval;
                Ck = Theodorsen(k);
                
                K_aer = U^2 * rho * (Kae + Bw * real(Ck) * Cw - k/b * Bw * imag(Ck) * Cw_cap - k^2/b^2 * Mae);
                C_aer = U * rho * (Cae + b/k * Bw * imag(Ck) * Cw + Bw * real(Ck) * Cw_cap);
                
                K_aer_full = I_wh*K_aer(1,1) + I_wt*K_aer(1,2) + I_th*K_aer(2,1) + I_tt*K_aer(2,2);
                C_aer_full = I_wh*C_aer(1,1) + I_wt*C_aer(1,2) + I_th*C_aer(2,1) + I_tt*C_aer(2,2);
                
                A2_new = [zeros(No), eye(No); -diag(f(1:No).^2) + K_aer_full, C_aer_full];
                
                [eigenvectors(:,j,i), eig_diag] = eigs(A2_new, 1, eigenval_0);
                eigenval = diag(eig_diag);
                eigenvalues(j,i) = eigenval;
                n_iter = n_iter + 1;
            end
        end
    end
    
    %% Step 8: Find Flutter Speed (Zero Damping Crossing)
    flutternotfound = true;
    for j = 1:length(U_range)
        for i = 1:size(eigenvalues, 1) 
            if real(eigenvalues(i, j)) >= 1e-6 && flutternotfound
                flutternotfound = false;
                idx_U = j;
                idx_mode = i;
            end
        end
    end
    
    % Linear interpolation to find precise crossing speed
    m_slope = (real(eigenvalues(idx_mode, idx_U-1)) - real(eigenvalues(idx_mode, idx_U))) / (U_range(idx_U-1) - U_range(idx_U));
    n_intercept = real(eigenvalues(idx_mode, idx_U-1)) - m_slope * U_range(idx_U-1);
    Uflut = -n_intercept / m_slope;
end