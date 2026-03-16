%% =======================================================================
% MAIN SCRIPT: Workshop 1 - Task 7 (Aileron Effects)
% DESCRIPTION: Comprehensive aeroelastic analysis of the Goland Wing 
%              equipped with a trailing-edge aileron.
%              Uses EXACT SYMBOLIC integration for Rayleigh-Ritz matrices.
%
% AUTHORS: Alejandro Rivera Míguez
%          (Based on course material by G. Quaranta & L. Mapelli)
%% =======================================================================

clc
clear all
y = sym('y');
addpath("Auxiliary Functions\")

%% 1. PARAMETERS (GOLAND + AILERON)
orden = 6;
L = 6.096; chord = 1.829; b = chord/2;

% Structural Properties
EJ = 9.77e6; GJ = 9.876e5;
m = 35.75; Jp = 8.65; xt_n = 0.1*chord;

% Flap Geometry
x0 = 3.35; xf = 5.15; % (Approx 1.8m span)
L_flap_span = xf - x0;
cf = 0.45;

% Flap Mass and Stiffness
mf_dist = 8.92 / 1.8;
K_beta = 6.48e3;

% Section Geometry
x_EA = 0.33 * chord;
x_hinge = chord - cf;
x_cg_flap = x_hinge + cf/2;

% Distances
d_cg_hinge = x_cg_flap - x_hinge;
d_cg_ea    = x_cg_flap - x_EA;
d_hinge_ea = x_hinge - x_EA;

% Inertias and Static Unbalances
S_beta = mf_dist * d_cg_hinge;
I_beta = (0.09/1.8) + mf_dist * d_cg_hinge^2;
S_flap_on_wing = mf_dist * d_cg_ea;
I_flap_on_wing = (0.09/1.8) + mf_dist * d_cg_ea^2;

% Aerodynamics
Cla = 2*pi; rho = 1.225;
ep = -0.34;
cp = (x_hinge - b)/b; 

%% 2. RAYLEIGH-RITZ SHAPE FUNCTIONS
for n = 1:orden
    phi(n) = sin(n*pi/(2*L)*y);
    psi(n) = 1-cos(n*pi/(2*L)*y);
end

%% 3. STRUCTURAL MATRICES (M and K)
% Derivatives of shape functions
psi_pp = diff(psi, y, 2); 
phi_p = diff(phi, y, 1);

% Stiffness Matrices (Symbolic Integration)
Kww = double(int(psi_pp'*EJ*psi_pp, y, 0, L));
Ktt = double(int(phi_p'*GJ*phi_p, y, 0, L));
Kbb = K_beta;

% Mass Matrices (Symbolic Integration)
Mww = double(int(psi'*m*psi, y, 0, L)) + double(int(psi'*mf_dist*psi, y, x0, xf));
Mtt = double(int(phi'*Jp*phi, y, 0, L)) + double(int(phi'*I_flap_on_wing*phi, y, x0, xf));
Mbb = I_beta * L_flap_span;

% Coupling Mass Terms
Mwt = -double(int(psi'*m*xt_n*phi, y, 0, L)) - double(int(psi'*S_flap_on_wing*phi, y, x0, xf));
Mtw = Mwt';
Mwb = -double(int(psi'*S_beta, y, x0, xf)); 
Mbw = Mwb';

term_coupl = 0.5 * mf_dist * d_cg_ea * d_cg_hinge + (0.09/1.8); 
Mtb = double(int(phi'*term_coupl, y, x0, xf)); 
Mbt = Mtb';

% Global Structural Matrices Assembly
M = [Mww, Mwt, Mwb; 
     Mtw, Mtt, Mtb; 
     Mbw, Mbt, Mbb];
     
K = [Kww,          zeros(orden), zeros(orden,1); 
     zeros(orden), Ktt,          zeros(orden,1); 
     zeros(orden,1)', zeros(orden,1)', Kbb ];

%% 4. STRUCTURAL EIGENPROBLEM
% Solve for natural frequencies and mode shapes
[Autovectores, Autovalores] = eig(M\K);
aux = Autovectores;
[Autovalores, sortng] = sort(sqrt(diag(Autovalores)), 'ascend');

% Sort eigenvectors to match frequency order
for j=1:length(Autovalores)
    Autovectores(:, j) = aux(:, sortng(j));
end
aux = Autovectores;
aux2 = Autovalores;

%% 5. THEODORSEN PARAMETERS
T = TheodorsenCoefficients(cp,ep);

%% 6. AERODYNAMIC INTEGRAL SHAPE FUNCTIONS
% Pre-computing integrals for aeroelastic generalized forces
WW = double(int(psi'*psi, y, 0, L));
TT = double(int(phi'*phi, y, 0, L));
WT = double(int(psi'*phi, y, 0, L));
TW = double(int(phi'*psi, y, 0, L));

% Integrals specific to the flap span
WB = double(int(psi', y, x0, xf));
BW = double(int(psi, y, x0, xf));
TB = double(int(phi', y, x0, xf));
BT = double(int(phi, y, x0, xf));
BB = double(int(1, y, x0, xf));

%% 7. SOLVER DEFINITIONS
tolerancia = 10e-4;
Uvec = linspace(1, 300, 100);

% Initialize State-Space Matrix
Amatrix = zeros(4*length(phi)+2,  4*length(phi)+2);
Amatrix (1:(2*length(phi)+1), (2*length(phi)+2):end) = eye(2*length(psi)+1);
tol = 100;
counter1 = 0;

%% 8. AEROELASTIC LOOP (p-k METHOD)
for modes=1:2*length(psi)
    for jj=1:length(Uvec)
        U = Uvec(jj);
        k = Autovalores(modes)*b/U;
        
        while tol>tolerancia
            % Unsteady Aerodynamics calculation
            Ck = besselh(1, 2, k)/(besselh(1, 2, k) + 1i*besselh(0, 2, k));
            
            % Strip theory aero matrices (Mass, Damping, Stiffness)
            Ma = rho*b^2*[pi,-pi*b*ep,-T(1)*b;
                        ep*pi*b,-pi*b^2*(1/8+ep^2),(T(7)+(0.5079-ep)*T(1))*b^2;
                        T(1)*b,-2*T(13)*b^2,1/pi*T(3)*b^2];
            
            Ca = rho*b^2*[0,pi,-T(4);
                            0,-pi*b*(1/2-ep),-b*(T(1)-T(8)-(0.5079-ep)*T(4)+T(11)/2);
                            0, -b*(-2*T(9)-T(1)+T(4)*(ep-1/2)),b/(2*pi)*T(4)*T(11)];
                            
            Ka = rho*b^2*[0,0,0;
                          0,0,-(T(4)+T(10));
                          0,0,-(1/pi)*(T(5)-T(4)*T(10))];
            
            Bw = rho*b*Cla*[1;b*(ep+1/2);-b*T(12)/(2*pi)];
            Cw = [0, 1, T(10)/pi];
            Chatw = [1,b*(1/2-ep),b*T(11)/(2*pi)];
            
            % Combined 2D aero matrices
            Kaero = (Ka + Bw*Ck*Cw);
            Caero = (Ca+Bw*Ck*Chatw);
            Maero = Ma;
        
            % Generalized aero matrices (Projected onto modes)
            KK = U^2*[-Kaero(1,1)*WW, Kaero(1,2)*WT, Kaero(1,3)*WB; -Kaero(2,1)*TW, Kaero(2,2)*TT, Kaero(2,3)*TB; -Kaero(3, 1)*BW, Kaero(3,2)*BT, Kaero(3,3)*BB];
            CC = U*[-Caero(1,1)*WW, Caero(1,2)*WT, Caero(1,3)*WB; -Caero(2,1)*TW, Caero(2,2)*TT, Caero(2,3)*TB; -Caero(3, 1)*BW, Caero(3,2)*BT, Caero(3,3)*BB];
            MM = [-Maero(1,1)*WW, Maero(1,2)*WT, Maero(1,3)*WB; -Maero(2,1)*TW, Maero(2,2)*TT, Maero(2,3)*TB; -Maero(3, 1)*BW, Maero(3,2)*BT, Maero(3,3)*BB];
            
            % System State-Space Components
            TM2 = M -MM;
            TM1 = -CC;
            TM0 = K-KK;
            
            % Populate State-Space Matrix
            Amatrix((2*length(psi)+2):end, 1:(2*length(psi)+1)) = -TM2\TM0;
            Amatrix((2*length(psi)+2):end, (2*length(psi)+2):end) = -TM2\TM1;
            
            [AVEC, z] = eig(Amatrix);
            
            % Initialization of tracking vector
            if counter1 == 0
                Autovectores_track = AVEC((2*length(phi)+2):end, :);
                counter1 = 1;
           else 
               Autovectores_track = AVEC;
           end
           
           % Modal Tracking using complex scalar product
           for tracking=1:length(Autovectores)
               ScalarProducts(tracking) = ps_complex(Autovectores_track(:, tracking), Autovectores(:, modes));
           end
           
           % Select matched mode
           [aux1, Index2] = max(ScalarProducts); 
           z = diag(z);
           z = z(Index2); % Eigenvalue of the mode we are studying
           
           % Check for convergence
           diff1= abs(Autovalores(modes) - sqrt(imag(z)^2 + real(z)^2)); 
           Autovalores(modes) = sqrt(imag(z)^2 + real(z)^2);
           k = Autovalores(modes)*b/U;
           damping = real(z)/Autovalores(modes);
           tol = diff1;
           
           % Redefinition of eigenvectors for next iteration of convergence
           Autovectores(1:(2*length(phi)+1), modes) = AVEC(1:(2*length(phi)+1), Index2);
           Autovectores((2*length(phi)+2):(4*length(phi)+2), modes) = AVEC((2*length(phi)+2):end, Index2);
    
        end
        tol = 100;
        k_box (jj, modes) = k;
        damping_box (jj, modes) = damping;
    end
    counter1 = 0;
end

%% 9. FLUTTER SPEED DETECTION
flutternotfound = true;
for jj=1:length(Uvec)
    for kkk =1:2*length(phi)
        % Detect crossing of stability boundary (damping > 0)
        if damping_box(jj, kkk)>0 && flutternotfound
            jjj = jj;
            kkkk = kkk;
            flutternotfound = false;
        end
    end
end

if flutternotfound 
    U_flutter = NaN;
else
    % Linear interpolation for precise flutter speed
    m = (real(damping_box(jjj-1, kkkk)) - real(damping_box(jjj, kkkk)))/(Uvec(jjj-1) - Uvec(jjj));
    n = real(damping_box(jjj-1, kkkk)) - m*Uvec(jjj-1);
    U_flutter = -n/m
end

%% 10. V-g DIAGRAM PLOTTING
figure(101); clf;
U_plot = Uvec;  
g_data = 2 * damping_box';  
hold on; grid on; box on;

% Plot damping for all tracked modes
for m = 1:12
    plot(U_plot, g_data(m, :), 'LineWidth', 1.5);
end

% Zero stability line
yline(0, 'k--', 'LineWidth', 1.5);

% Mark flutter speed
if exist('U_flutter', 'var')
    plot([U_flutter, U_flutter], ylim, 'r--', 'LineWidth', 1.5);
    plot(U_flutter, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
end

ylabel('Damping $g$ [-]', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
xlim([0, max(U_plot)]);
ylim([-0.15, 0.05]);
fprintf('Flutter speed: %.2f m/s\n', U_flutter);

%% 11. CONTROL REVERSAL ANALYSIS
fprintf('\n=== CONTROL REVERSAL ANALYSIS ===\n');

UVec_CR = linspace(0, 200, 400);
rolling_moment = zeros(length(UVec_CR), 1);

% Discretization for numerical integration
n_points = 100;
y_num = linspace(0, L, n_points)';
psi_num = zeros(n_points, orden);
phi_num = zeros(n_points, orden);

for n = 1:orden
    phi_num(:, n) = sin(n*pi/(2*L)*y_num);
    psi_num(:, n) = 1 - cos(n*pi/(2*L)*y_num);
end

% Theodorsen coefficients for k=0 (static, C(k)=1)
Ck = 1;
T = TheodorsenCoefficients(cp, ep);

% Aerodynamic matrices for static case (k=0)
Ka_static = rho*b^2*[0, 0, 0;
                     0, 0, -(T(4)+T(10));
                     0, 0, -(1/pi)*(T(5)-T(4)*T(10))];
                
Bw_static = rho*b*Cla*[1; b*(ep+1/2); -b*T(12)/(2*pi)];
Cw_static = [0, 1, T(10)/pi];

% Generalized aerodynamic stiffness matrix
Kaero_static = (Ka_static + Bw_static*Ck*Cw_static);
KK_static = @(U) U^2 * [-Kaero_static(1,1)*WW, Kaero_static(1,2)*WT, Kaero_static(1,3)*WB; 
                        -Kaero_static(2,1)*TW, Kaero_static(2,2)*TT, Kaero_static(2,3)*TB; 
                        -Kaero_static(3,1)*BW, Kaero_static(3,2)*BT, Kaero_static(3,3)*BB];

% Compute rolling moment vs velocity
for i = 1:length(UVec_CR)
    U = UVec_CR(i);
    
    % Aeroelastic stiffness matrix
    K_aeroelastic = K - KK_static(U);
    
    % Force vector from unit aileron deflection
    delta = 1;  % Unit deflection [rad]
    F_beta = zeros(size(K, 1), 1);
    F_beta(end) = K_beta * delta;
    
    % Static solution
    q_static = K_aeroelastic \ F_beta;
    
    % Extract modal coefficients
    q_theta = q_static(orden+1:2*orden);  % Torsion coefficients
    q_beta = q_static(end);               % Aileron deflection
    
    % Calculate physical torsion distribution
    theta_dist = phi_num * q_theta;  
    in_aileron = (y_num >= x0) & (y_num <= xf);
    
    % Effective angle of attack
    alpha_eff = theta_dist;
    alpha_eff(in_aileron) = alpha_eff(in_aileron) + (T(10)/pi) * q_beta;
    
    % Lift distribution (strip theory)
    L_dist = 0.5 * rho * U^2 * chord * Cla * alpha_eff;
    
    % Rolling moment: integral of lift × moment arm
    rolling_moment(i) = trapz(y_num, L_dist .* y_num);
end

% Find reversal speed (where rolling moment crosses zero)
reversal_found = false;
U_reversal = NaN;
for i = 2:length(UVec_CR)
    if rolling_moment(i-1) * rolling_moment(i) <= 0 && rolling_moment(i-1) ~= 0
        % Linear interpolation
        U1 = UVec_CR(i-1);
        U2 = UVec_CR(i);
        M1 = rolling_moment(i-1);
        M2 = rolling_moment(i);
        
        U_reversal = U1 - M1 * (U2 - U1) / (M2 - M1);
        reversal_found = true;
        break;
    end
end

figure(24); clf;
hold on; grid on; box on;
plot(UVec_CR, rolling_moment, 'b-', 'LineWidth', 2);
xlabel('Airspeed $U_\infty$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$M_{\mathrm{roll}}$ [N$\cdot$m]', 'Interpreter', 'latex', 'FontSize', 14);
yline(0, 'k--', 'LineWidth', 1.5);

if reversal_found
    plot(U_reversal, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
    xline(U_reversal, 'r--', 'LineWidth', 1.5);
    fprintf('Control reversal speed: %.1f m/s\n', U_reversal);
else
    fprintf('No control reversal found in the velocity range.\n');
end
xlim([0, 250]);
ylim([min(rolling_moment)*1.1, max(rolling_moment)*1.1]);

%% 12. DIVERGENCE ANALYSIS
fprintf('\n=== DIVERGENCE ANALYSIS ===\n');
U_div_range = linspace(0, 400, 200);
det_values = zeros(size(U_div_range));
min_eig_values = zeros(size(U_div_range));

% For divergence, we use static aerodynamics (k=0, C(k)=1)
Ck_div = 1;
Kaero_static = (Ka_static + Bw_static*Ck_div*Cw_static);
KK_static_div = @(U) U^2 * [-Kaero_static(1,1)*WW, Kaero_static(1,2)*WT, Kaero_static(1,3)*WB; 
                            -Kaero_static(2,1)*TW, Kaero_static(2,2)*TT, Kaero_static(2,3)*TB; 
                            -Kaero_static(3,1)*BW, Kaero_static(3,2)*BT, Kaero_static(3,3)*BB];

divergence_found = false;
U_divergence = NaN;

for i = 1:length(U_div_range)
    U = U_div_range(i);
    K_aeroelastic = K - KK_static_div(U);
    
    % Eigenvalues of the stiffness matrix
    eig_vals = eig(K_aeroelastic);
    min_eig_values(i) = min(real(eig_vals));
    det_values(i) = det(K_aeroelastic);
    
    % Check for divergence (minimum eigenvalue crosses zero)
    if i > 1 && min_eig_values(i-1) > 0 && min_eig_values(i) <= 0
        U1 = U_div_range(i-1);
        U2 = U_div_range(i);
        eig1 = min_eig_values(i-1);
        eig2 = min_eig_values(i);
        
        U_divergence = U1 - eig1 * (U2 - U1) / (eig2 - eig1);
        divergence_found = true;
    end
end

figure(26); clf;
set(gcf, 'Position', [100, 100, 700, 500]);
subplot(2,1,1); hold on; grid on; box on;
plot(U_div_range, min_eig_values, 'b-', 'LineWidth', 2);
xlabel('$U$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Min Eigenvalue of $K_{aero}$', 'Interpreter', 'latex', 'FontSize', 12);
yline(0, 'k--', 'LineWidth', 1.5);

if divergence_found
    plot(U_divergence, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xline(U_divergence, 'r--', 'LineWidth', 1.5);
    text(U_divergence + 5, max(min_eig_values)*0.1, ...
         sprintf('$U_{div} = %.1f$ m/s', U_divergence), ...
         'Interpreter', 'latex', 'FontSize', 12, ...
         'BackgroundColor', 'white');
end

subplot(2,1,2); hold on; grid on; box on;
plot(U_div_range, det_values, 'r-', 'LineWidth', 2);
xlabel('$U$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('det$(K_{aero})$', 'Interpreter', 'latex', 'FontSize', 12);
yline(0, 'k--', 'LineWidth', 1.5);
if divergence_found
    xline(U_divergence, 'r--', 'LineWidth', 1.5);
end

%% 13. RESULTS SUMMARY
fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('Flutter speed: %.1f m/s\n', U_flutter);
if reversal_found
    fprintf('Control reversal: %.1f m/s\n', U_reversal);
end
if divergence_found
    fprintf('Divergence speed: %.1f m/s\n', U_divergence);
else
    fprintf('No divergence found below %.0f m/s\n', max(U_div_range));
end

% Determine which phenomenon limits the flight envelope
fprintf('\n=== FLIGHT ENVELOPE LIMITATION ===\n');
if reversal_found && divergence_found
    speeds = [U_flutter, U_reversal, U_divergence];
    [min_speed, idx] = min(speeds);
    
    if idx == 1
        fprintf('FLUTTER limits the flight envelope at %.1f m/s\n', U_flutter);
    elseif idx == 2
        fprintf('CONTROL REVERSAL limits the flight envelope at %.1f m/s\n', U_reversal);
    else
        fprintf('DIVERGENCE limits the flight envelope at %.1f m/s\n', U_divergence);
    end
    
elseif reversal_found
    if U_reversal < U_flutter
        fprintf('CONTROL REVERSAL limits the flight envelope at %.1f m/s\n', U_reversal);
    else
        fprintf('FLUTTER limits the flight envelope at %.1f m/s\n', U_flutter);
    end
else
    fprintf('FLUTTER limits the flight envelope at %.1f m/s\n', U_flutter);
end

%% 14. COMPARATIVE PLOT
figure(27); clf;
set(gcf, 'Position', [100, 100, 700, 500]);
if reversal_found && divergence_found
    speeds = [U_flutter, U_reversal, U_divergence];
    labels = {'Flutter', 'Control Reversal', 'Divergence'};
    bar_colors = [0.8, 0.2, 0.2; 0.2, 0.4, 0.8; 0.2, 0.8, 0.4];
    
    h_bar = bar(speeds);
    set(h_bar, 'FaceColor', 'flat');
    for i = 1:3
        h_bar.CData(i,:) = bar_colors(i,:);
    end
    set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'latex');
    ylabel('Speed [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
    
    for i = 1:3
        text(i, speeds(i) + 5, sprintf('%.1f', speeds(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    [min_speed, idx] = min(speeds);
    text(idx, min_speed - 10, 'LIMITING', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white');
     
elseif reversal_found
    speeds = [U_flutter, U_reversal];
    labels = {'Flutter', 'Control Reversal'};
    bar_colors = [0.8, 0.2, 0.2; 0.2, 0.4, 0.8];
    
    h_bar = bar(speeds);
    set(h_bar, 'FaceColor', 'flat');
    for i = 1:2
        h_bar.CData(i,:) = bar_colors(i,:);
    end
    set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'latex');
    ylabel('Speed [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
    
    for i = 1:2
        text(i, speeds(i) + 5, sprintf('%.1f', speeds(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    [min_speed, idx] = min(speeds);
    text(idx, min_speed - 10, 'LIMITING', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white');
else
    speeds = U_flutter;
    labels = {'Flutter'};
    bar_colors = [0.8, 0.2, 0.2];
    
    h_bar = bar(speeds);
    set(h_bar, 'FaceColor', bar_colors);
    set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'latex');
    ylabel('Speed [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
    text(1, speeds + 5, sprintf('%.1f', speeds), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
grid on; box on;
ylim([0, max([U_flutter, U_reversal, U_divergence]) * 1.2]);


%% =======================================================================
% HELPER FUNCTION: ps_complex
% DESCRIPTION: Computes a normalized complex scalar product robustly.
% =======================================================================
function Res = ps_complex(X, Y)
    S1 = 0;
    S2 = 0;
    S3 = 0;
    S4 = 0;
    for jjj=1:length(X)
        S1 = S1 + real(X(jjj))*real(Y(jjj)) + imag(X(jjj))*imag(Y(jjj));
        S2 = S2 + real(X(jjj))*imag(Y(jjj)) - imag(X(jjj))*real(Y(jjj));
        S3 = norm(X(jjj))^2 + S3;
        S4 = norm(Y(jjj))^2 + S4;
    end
    Res = sqrt(S1^2 + S2^2)/sqrt(S3*S4);
end