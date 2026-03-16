%% =======================================================================
% FUNCTION: NoKB_FEM_FlapCoupledBendingTorsion
% DESCRIPTION: FEM model for the Goland wing WITHOUT actuator stiffness 
%              (KBeta = 0). Solves the coupled bending-torsion-flap problem.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================
function [fr, qr, FEMdata] = NoKB_FEM_FlapCoupledBendingTorsion(Ne, Nd)
    % --- Physical Parameters ---
    % Wing
    L = 6.096;        % [m]
    c = 1.829;        % [m]
    shat = 0.1;       % [-] Distance between elastic and mass axes
    EI = 9.77 * 1e6;  % [N*m2]
    m = 35.72;        % [kg/m]
    GJ = 987600;      % [N m^2]
    Io = 7.452;       % [kg m] 
    s = shat * c;
    Itheta = Io + m * s^2;
    xe = 0.33 * c;
    
    % Flap
    Mf = 8.92;            % [kg]
    If = 0.09;            % [kg m^2]
    cf = 0.45;            % [m]
    yL = 3.35;            % [m]
    yU = yL + 1.8;        % [m]
    yf = (yL + yU) / 2;   % [m]
    xf = c - cf / 2;      % [m]
    Ih = If + (cf / 2)^2 * Mf;  % [kg m^2]
    KBeta = 0;            % [Nm/rad] Actuator stiffness
    
    % --- Mesh Definitions ---
    NnodesB = Ne + 1;
    NdofsB = 2 * NnodesB;
    NnodesT = 2 * Ne + 1;
    NdofsT = NnodesT;
    Ndofs = NdofsB + NdofsT + 1; % Overall DOFs including Beta
    deltay = L / Ne;
    
    % Define connectivity matrix (7 DOFs per element: 4 bending, 3 torsion)
    Connect = nan(Ne, 7);
    for i = 1:Ne
        Connect(i, 1:4) = [1+(i-1)*2, 2+(i-1)*2, 3+(i-1)*2, 4+(i-1)*2];
        Connect(i, 5:7) = NdofsB + [1+2*(i-1), 2+2*(i-1), 3+2*(i-1)];
    end
    
    % --- Element Matrices ---
    % Bending problem
    Kww_local = EI/deltay^3 * [12, 6*deltay, -12, 6*deltay; 6*deltay, 4*deltay^2, -6*deltay, 2*deltay^2; ...
                               -12, -6*deltay, 12, -6*deltay; 6*deltay, 2*deltay^2, -6*deltay, 4*deltay^2];
    Mww_local = m*deltay/420 * [156, 22*deltay, 54, -13*deltay; 22*deltay, 4*deltay^2, 13*deltay, -3*deltay^2; ...
                                54, 13*deltay, 156, -22*deltay; -13*deltay, -3*deltay^2, -22*deltay, 4*deltay^2];
    % Torsion problem
    Ktt_local = GJ/deltay/3 * [7, -8, 1; -8, 16, -8; 1, -8, 7];
    Mtt_local = Itheta*deltay/30 * [4, 2, -1; 2, 16, 2; -1, 2, 4];
    
    % Mass coupling
    Mwt_local = -m*deltay*s/60 * [11, 20, -1; deltay, 4*deltay, 0; -1, 20, 11; 0, -4*deltay, -deltay];
    
    % Coupled Bending-Torsion Element
    K_local = [Kww_local, zeros(4,3); zeros(3,4), Ktt_local];
    M_local = [Mww_local, Mwt_local; Mwt_local', Mtt_local]; 
    
    % --- Global Assembly ---
    KK = zeros(Ndofs, Ndofs);
    MM = zeros(Ndofs, Ndofs);
    for i = 1:Ne
        KK(Connect(i,:), Connect(i,:)) = KK(Connect(i,:), Connect(i,:)) + K_local;
        MM(Connect(i,:), Connect(i,:)) = MM(Connect(i,:), Connect(i,:)) + M_local;
    end
    
    % --- Flap Contributions ---
    EID = floor(yf/deltay) + 1; 
    NIDs = Connect(EID, :);    
    
    % Shape functions at flap CG
    xi = (yf - (EID - 1)*deltay) / deltay; 
    psiw1 = 1 - 3*xi^2 + 2*xi^3;
    psiw2 = (xi - 2*xi^2 + xi^3)*deltay;
    psiw3 = 3*xi^2 - 2*xi^3;
    psiw4 = (-xi^2 + xi^3)*deltay;
    psiw = [psiw1, psiw2, psiw3, psiw4]';
    
    psit1 = 1 - 3*xi + 2*xi^2;
    psit2 = 4*xi - 4*xi^2;
    psit3 = 2*xi^2 - xi;
    psit = [psit1, psit2, psit3]';
    
    % Add Flap mass coupling
    MM(NIDs(1:4), NIDs(1:4)) = MM(NIDs(1:4), NIDs(1:4)) + psiw * Mf * psiw';
    MM(NIDs(5:7), NIDs(5:7)) = MM(NIDs(5:7), NIDs(5:7)) + psit * (If + (xf - xe)^2 * Mf) * psit';
    MM(NIDs(1:4), NIDs(5:7)) = MM(NIDs(1:4), NIDs(5:7)) - psiw * Mf * (xf - xe) * psit';
    MM(NIDs(5:7), NIDs(1:4)) = MM(NIDs(5:7), NIDs(1:4)) - psit * Mf * (xf - xe) * psiw';
    
    % Coupling with flap deflection (beta)
    MM(NIDs(1:4), end) = MM(NIDs(1:4), end) - psiw * Mf * cf / 2;
    MM(NIDs(5:7), end) = MM(NIDs(5:7), end) + psit * (If + Mf * cf / 2 * (xf - xe));
    
    % Beta equation direct terms
    MM(end, end) = Ih;
    KK(end, end) = KBeta;
    
    MM(end, NIDs(1:4)) = MM(end, NIDs(1:4)) - Mf * cf / 2 * psiw';
    MM(end, NIDs(5:7)) = MM(end, NIDs(5:7)) + (If + Mf * cf / 2 * (xf - xe)) * psit';
    
    % --- Apply Boundary Conditions ---
    % Eliminate torsion of the root (first torsional dof)
    MMr = MM([1:NdofsB, NdofsB+2:end], [1:NdofsB, NdofsB+2:end]);
    KKr = KK([1:NdofsB, NdofsB+2:end], [1:NdofsB, NdofsB+2:end]);
    
    % Eliminate vertical displacement and bending angle of the root
    MMr = MMr(3:end, 3:end);
    KKr = KKr(3:end, 3:end);
    
    % --- Solve Eigenproblem ---
    [Vr, D] = eig(-KKr, MMr);
    [Frequencies, sorting] = sort(imag(sqrt(diag(D))));
    Vr = Vr(:, sorting); % Bona fide MATLAB magic to order the matrix by columns
    
    % Reintroduce eliminated DOFs
    V = [zeros(2, Ndofs - 3); Vr];
    V = [V(1:NdofsB, :); zeros(1, Ndofs - 3); V(NdofsB+1:end, :)];
    
    fr = Frequencies(1:Nd);
    qr = nan(Ndofs, Nd);
    
    % Return eigenvectors normalized by mass matrix
    for i = 1:Nd
        % WARNING: Original authors noted a potential normalization issue here!
        qr(:, i) = V(:, i) / sqrt(V(:, i)' * MM * V(:, i)); 
    end
    
    FEMdata.deltay = deltay;
    FEMdata.Connect = Connect; 
end