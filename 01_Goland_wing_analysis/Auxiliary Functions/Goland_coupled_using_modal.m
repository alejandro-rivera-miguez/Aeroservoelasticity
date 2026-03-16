%% =======================================================================
% FUNCTION: Goland_coupled_using_modal
% DESCRIPTION: Computes coupled natural frequencies of the Goland wing 
%              using a Rayleigh-Ritz modal synthesis approach with 
%              uncoupled bending and torsion modes.
%
% AUTHOR: Alejandro Rivera Míguez (Adapted for Aeroservoelasticity course)
% =======================================================================
function f = Goland_coupled_using_modal(n_modb, n_modt)
    % 1. Fixed grid for integration (high fidelity)
    n_elem = 50; 
    x = linspace(0, 1, n_elem + 1);
    
    %% Step 1: Compute Uncoupled Basis Functions
    % Bending 
    EJ = 9.77e6 * ones(size(x)); 
    m  = 35.75 * ones(size(x));
    L  = 6.096;
    cdofs_b = [1, 2];
    [fb_unc, Modesw] = FEM_beam_bending(x, EJ(1:end-1), m(1:end-1), L, cdofs_b, n_modb);
    
    % Torsion 
    GJ = 9.876e5 * ones(size(x));
    Jp = 8.65  * ones(size(x));
    cdofs_t = 1;
    [ft_unc, Modest] = FEM_beam_torsion(x, GJ(1:end-1), Jp(1:end-1), L, cdofs_t, n_modt);
    
    %% Step 2: Assemble Generalized Matrices
    % Mass Coupling Term (Static Unbalance)
    chord = 1.829;
    xt    = 0.1 * chord; % Constant unbalance
    m_val = 35.75;       % Constant mass
    
    % Mbt(i,j) = Integral (Phi_h * m * x_alpha * Phi_alpha)
    Mbt = zeros(n_modb, n_modt);
    for i = 1:n_modb
        for j = 1:n_modt
            integrand = Modesw(i, :) .* m_val .* xt .* Modest(j, :);
            Mbt(i,j) = trapz(x, integrand) * L; 
        end
    end
    
    %% Step 3: Global System Construction
    n_total = n_modb + n_modt;
    
    % Stiffness Matrix (K) - Diagonal because uncoupled modes are orthogonal
    K_gen = zeros(n_total, n_total);
    K_gen(1:n_modb, 1:n_modb) = diag(fb_unc.^2);
    K_gen(n_modb+1:end, n_modb+1:end) = diag(ft_unc.^2);
    
    % Mass Matrix (M) - Block diagonal is Identity (modes are mass-normalized)
    M_gen = eye(n_total); 
    M_gen(1:n_modb, n_modb+1:end) = Mbt;
    M_gen(n_modb+1:end, 1:n_modb) = Mbt';
    
    %% Step 4: Solve Coupled Eigenvalues
    [~, D] = eig(K_gen, M_gen);
    
    % Sort ASCENDING (Low freq -> High freq)
    [eigenvalues, ~] = sort(diag(D), 'ascend');
    
    % Frequencies in rad/s (ensure real output)
    f = real(sqrt(eigenvalues)); 
end