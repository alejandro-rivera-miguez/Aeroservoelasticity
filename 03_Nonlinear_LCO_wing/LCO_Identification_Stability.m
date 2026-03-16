%% =======================================================================
% SCRIPT: LCO_Identification_Stability
% DESCRIPTION: Performs LCO identification for the Goland wing with 
%              free-play using a continuation algorithm. It ALSO performs 
%              a local dynamic stability check on the identified LCOs 
%              using the exact derivative of the Theodorsen function.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================

clc; clear all; close all;

addpath("Auxiliary Functions\")

% --- Solver Parameters ---
tol = 1e-8;
maxiter = 100;
Nd = 2;
Ne = 10;

% --- Free-play and Amplitude Range ---
Kbeta = 6.48 * 1e3;     % [Nm/rad]
delta = deg2rad(1);
Amax = 1.45 * delta;
A_v = linspace(delta, Amax, 100);

% --- Load flutter results as initial conditions ---
[UF, wF, qF] = NoKB_FlapPKGoland(Ne, Nd, 200, 1000, 0);

% --- Physical Parameters ---
L = 6.096;        % [m]
c = 1.829;        % [m]
cf = 0.45;        % [m]

b = c / 2;              % [m]
ehat = (0.33*c - b)/b;  % [-]
chat = (c - cf - b)/b;  % [-]
rho = 1.225;            % [kg/m^3]
Clalpha = 2 * pi;       % [1/rad]
yL = 3.35;
yU = 3.35 + 1.8;

% --- Theodorsen Constants ---
mu = acos(chat);
T1 = -1/3 * sqrt(1 - chat^2) * (2 + chat^2) + chat * mu;
T3 = -(1/8 + chat^2) * mu^2 + 1/4 * chat * sqrt(1 - chat^2) * mu * (7 + 2*chat^2) - 1/8 * (1 - chat^2) * (5*chat^2 + 4);
T4 = -mu + chat * sqrt(1 - chat^2);
T5 = -(1 - chat^2) - mu^2 + 2 * chat * sqrt(1 - chat^2) * mu;
T7 = -(1/8 + chat^2) * mu + 1/8 * chat * sqrt(1 - chat^2) * (7 + 2*chat^2);
T8 = -1/3 * sqrt(1 - chat^2) * (2*chat^2 + 1) + chat * mu;
T9 = 1/2 * (1/3 * sqrt(1 - chat^2)^3 + ehat * T4);
T10 = sqrt(1 - chat^2) + mu;
T11 = mu * (1 - 2*chat) + sqrt(1 - chat^2) * (2 - chat);
T12 = sqrt(1 - chat^2) * (2 + chat) - mu * (2*chat + 1);
T13 = -1/2 * (T7 + (chat - ehat) * T1);

Ca = rho*b^2 * [0, pi, -T4; 0, -pi*b*(1/2 - ehat), -b*(T1 - T8 - (chat - ehat)*T4 + T11/2); 0, -b*(-2*T9 - T1 + T4*(ehat - 1/2)), b*T4*T11/2/pi];
Ka = rho*b^2 * [0, 0, 0; 0, 0, -(T4 + T10); 0, 0, -(T5 - T4*T10)/pi];
Ma = rho*b^2 * [pi, -pi*b*ehat, -b*T1; pi*b*ehat, -pi*b^2*(1/8 + ehat^2), b^2*(T7 + (chat - ehat)*T1); b*T1, -2*b^2*T13, b^2*T3/pi];
Bw = rho*b*Clalpha * [1; b*(ehat + 1/2); -T12*b/2/pi];
Cw = [0, 1, T10/pi];
Cwhat = [1, b*(1/2 - ehat), b*T11/2/pi];

% --- FEM Subroutine ---
[wstruct, qstruct, FEMdata] = NoKB_FEM_FlapCoupledBendingTorsion(Ne, Nd);
M = eye(Nd, Nd);
K = diag(wstruct.^2);
Cnl = qstruct(end, :);
Knl = Cnl' * Cnl;

modes_cell = cell(1, Nd);
for i = 1:Nd
    modes_cell{i} = @(y) [EvaluateCoupledBendingDisplacement(y, qstruct(:,i), FEMdata.deltay, FEMdata.Connect);...
                          EvaluateCoupledTorsionDisplacement(y, qstruct(:,i), FEMdata.deltay, FEMdata.Connect);
                          Betafun(y, qstruct(:,i))];
end

% --- Aerodynamic Integrals ---
Kstar1 = nan(Nd, Nd); Kstar2 = nan(Nd, Nd); Kstar3 = nan(Nd, Nd); Kstar4 = nan(Nd, Nd);
Cstar1 = nan(Nd, Nd); Cstar2 = nan(Nd, Nd); Cstar3 = nan(Nd, Nd);

for i = 1:Nd
    for j = 1:Nd
        Kstar1(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ma), 0, L, 'Waypoints', [yL, yU]);
        Kstar2(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cw), 0, L, 'Waypoints', [yL, yU]);
        Kstar3(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cwhat), 0, L, 'Waypoints', [yL, yU]);
        Kstar4(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ka), 0, L, 'Waypoints', [yL, yU]);
        
        Cstar1(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ca), 0, L, 'Waypoints', [yL, yU]);
        Cstar2(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cw), 0, L, 'Waypoints', [yL, yU]);
        Cstar3(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cwhat), 0, L, 'Waypoints', [yL, yU]);
    end
end

% --- Initialize storage ---
Ulc  = nan(size(A_v));
wlc  = nan(size(A_v));
qrlc = nan(Nd, length(A_v));
qilc = nan(Nd, length(A_v));
stability = nan(size(A_v));

U = UF;
w = wF;
qr = real(delta * qF);
qi = imag(delta * qF);

% --- First Iterant ---
NA = 0;
dNA = 0;

k = w * b / U;
Ck = besselh(1, 2, k) / (besselh(1, 2, k) + 1i * besselh(0, 2, k));
Ktheo = -(k/b)^2 * Kstar1 + real(Ck) * Kstar2 - k/b * imag(Ck) * Kstar3 + Kstar4;
Ctheo = Cstar1 + b/k * imag(Ck) * Cstar2 + real(Ck) * Cstar3;

error_v = [ (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qr + w*U*Ctheo*qi;
            (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qi - w*U*Ctheo*qr;
            Cnl*qr - A_v(1);
            Cnl*qi ];

iter = 0;
while norm(error_v) > tol && iter <= maxiter
    DD = [ (-2*U*Ktheo*qr + w*Ctheo*qi), (-2*M*w*qr + U*Ctheo*qi), (-M*w^2 + K + Knl*NA - U^2*Ktheo), w*U*Ctheo;
           (-2*U*Ktheo*qi - w*Ctheo*qr), (-2*M*w*qi - U*Ctheo*qr), -w*U*Ctheo, (-M*w^2 + K + Knl*NA - U^2*Ktheo);
           0, 0, Cnl, zeros(1, Nd);
           0, 0, zeros(1, Nd), Cnl ];
    dd = DD \ (-error_v);
    
    U  = U  + dd(1);
    w  = w  + dd(2);
    qr = qr + dd(3 : 3+Nd-1);
    qi = qi + dd(3+Nd : 3+2*Nd-1);
    
    k = w * b / U;
    Ck = besselh(1, 2, k) / (besselh(1, 2, k) + 1i * besselh(0, 2, k));
    Ktheo = -(k/b)^2 * Kstar1 + real(Ck) * Kstar2 - k/b * imag(Ck) * Kstar3 + Kstar4;
    Ctheo = Cstar1 + b/k * imag(Ck) * Cstar2 + real(Ck) * Cstar3;
    
    error_v = [ (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qr + w*U*Ctheo*qi;
                (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qi - w*U*Ctheo*qr;
                Cnl*qr - A_v(1);
                Cnl*qi ];
    iter = iter + 1;
end

Ulc(1) = U;
wlc(1) = w;
qrlc(:,1) = qr;
qilc(:,1) = qi;

% --- LCO STABILITY CHECK (Bessel Derivatives) ---
dH1dk = 1/2 * (besselh(0, 2, k) - besselh(2, 2, k));
dH0dk = 1/2 * (besselh(-1, 2, k) - besselh(1, 2, k));
dCkdk = 1i * (dH1dk * besselh(0, 2, k) - besselh(1, 2, k) * dH0dk) / ((besselh(1, 2, k) + 1i * besselh(0, 2, k))^2);

dCtheodk = -b/(k^2) * imag(Ck) * Cstar2 + b/k * imag(dCkdk) * Cstar2 + real(dCkdk) * Cstar3;
dKtheodk = -2*k/b * Kstar1 + real(dCkdk) * Kstar2 - 1/b * imag(Ck) * Kstar3 - k/b * imag(dCkdk) * Kstar3;

Astab = [ (-2*M*w - U*b*dKtheodk)*qr + (U*Ctheo + b*w*dCtheodk)*qi, U*Ctheo*qr - 2*M*w*qi, -M*w^2 + K + Knl*NA - U^2*Ktheo, U*w*Ctheo;
          (-U*Ctheo - b*w*dCtheodk)*qr + (-2*M*w - U*b*dKtheodk)*qi, 2*M*w*qr - U*Ctheo*qi, -U*w*Ctheo, -M*w^2 + K + Knl*NA - U^2*Ktheo;
          0, 0, Cnl, zeros(1, Nd);
          0, 0, zeros(1, Nd), Cnl ];
bstab = [-Knl*dNA*qr; -Knl*dNA*qi; 1; 0];
vstab = Astab \ bstab;

if vstab(2) < 0
    stability(1) = 1;
else
    stability(1) = 0;
end

% --- Continuation Loop ---
for i = 2:length(A_v)
    NA = Free_Play_SIDF(A_v(i), delta, Kbeta);
    dNA = Free_Play_dSIDFdA(A_v(i), delta, Kbeta);
    
    U = Ulc(i-1);
    w = wlc(i-1);
    qr = qrlc(:, i-1);
    qi = qilc(:, i-1);
    
    k = w * b / U;
    Ck = besselh(1, 2, k) / (besselh(1, 2, k) + 1i * besselh(0, 2, k));
    Ktheo = -(k/b)^2 * Kstar1 + real(Ck) * Kstar2 - k/b * imag(Ck) * Kstar3 + Kstar4;
    Ctheo = Cstar1 + b/k * imag(Ck) * Cstar2 + real(Ck) * Cstar3;
    
    DDA = [ (-2*U*Ktheo*qr + w*Ctheo*qi), (-2*M*w*qr + U*Ctheo*qi), (-M*w^2 + K + Knl*NA - U^2*Ktheo), w*U*Ctheo;
            (-2*U*Ktheo*qi - w*Ctheo*qr), (-2*M*w*qi - U*Ctheo*qr), -w*U*Ctheo, (-M*w^2 + K + Knl*NA - U^2*Ktheo);
            0, 0, Cnl, zeros(1, Nd);
            0, 0, zeros(1, Nd), Cnl ];
    DDb = [-Knl*qr*dNA; -Knl*qi*dNA; 1; 0];
    DDv = DDA \ DDb;
    
    U  = U  + DDv(1)*(A_v(i) - A_v(i-1));
    w  = w  + DDv(2)*(A_v(i) - A_v(i-1));
    qr = qr + DDv(3 : 3+Nd-1)*(A_v(i) - A_v(i-1));
    qi = qi + DDv(3+Nd : 3+2*Nd-1)*(A_v(i) - A_v(i-1));
    
    k = w * b / U;
    Ck = besselh(1, 2, k) / (besselh(1, 2, k) + 1i * besselh(0, 2, k));
    Ktheo = -(k/b)^2 * Kstar1 + real(Ck) * Kstar2 - k/b * imag(Ck) * Kstar3 + Kstar4;
    Ctheo = Cstar1 + b/k * imag(Ck) * Cstar2 + real(Ck) * Cstar3;
    
    error_v = [ (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qr + w*U*Ctheo*qi;
                (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qi - w*U*Ctheo*qr;
                Cnl*qr - A_v(i);
                Cnl*qi ];
                
    iter = 0;
    while norm(error_v) > tol && iter <= maxiter
        DD = [ (-2*U*Ktheo*qr + w*Ctheo*qi), (-2*M*w*qr + U*Ctheo*qi), (-M*w^2 + K + Knl*NA - U^2*Ktheo), w*U*Ctheo;
               (-2*U*Ktheo*qi - w*Ctheo*qr), (-2*M*w*qi - U*Ctheo*qr), -w*U*Ctheo, (-M*w^2 + K + Knl*NA - U^2*Ktheo);
               0, 0, Cnl, zeros(1, Nd);
               0, 0, zeros(1, Nd), Cnl ];
        dd = DD \ (-error_v);
        
        U  = U  + dd(1);
        w  = w  + dd(2);
        qr = qr + dd(3 : 3+Nd-1);
        qi = qi + dd(3+Nd : 3+2*Nd-1);        
        
        k = w * b / U;
        Ck = besselh(1, 2, k) / (besselh(1, 2, k) + 1i * besselh(0, 2, k));
        Ktheo = -(k/b)^2 * Kstar1 + real(Ck) * Kstar2 - k/b * imag(Ck) * Kstar3 + Kstar4;
        Ctheo = Cstar1 + b/k * imag(Ck) * Cstar2 + real(Ck) * Cstar3;
        
        error_v = [ (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qr + w*U*Ctheo*qi;
                    (-M*w^2 + K + Knl*NA - U^2*Ktheo)*qi - w*U*Ctheo*qr;
                    Cnl*qr - A_v(i);
                    Cnl*qi ];
        iter = iter + 1;
    end
    
    Ulc(i)  = U;
    wlc(i)  = w;
    qrlc(:,i) = qr;
    qilc(:,i) = qi;
    
    % --- Stability check at converged step ---
    dH1dk = 1/2 * (besselh(0, 2, k) - besselh(2, 2, k));
    dH0dk = 1/2 * (besselh(-1, 2, k) - besselh(1, 2, k));
    dCkdk = 1i * (dH1dk * besselh(0, 2, k) - besselh(1, 2, k) * dH0dk) / ((besselh(1, 2, k) + 1i * besselh(0, 2, k))^2);
    
    dCtheodk = -b/(k^2) * imag(Ck) * Cstar2 + b/k * imag(dCkdk) * Cstar2 + real(dCkdk) * Cstar3;
    dKtheodk = -2*k/b * Kstar1 + real(dCkdk) * Kstar2 - 1/b * imag(Ck) * Kstar3 - k/b * imag(dCkdk) * Kstar3;
    
    Astab = [ (-2*M*w - U*b*dKtheodk)*qr + (U*Ctheo + b*w*dCtheodk)*qi, U*Ctheo*qr - 2*M*w*qi, -M*w^2 + K + Knl*NA - U^2*Ktheo, U*w*Ctheo;
              (-U*Ctheo - b*w*dCtheodk)*qr + (-2*M*w - U*b*dKtheodk)*qi, 2*M*w*qr - U*Ctheo*qi, -U*w*Ctheo, -M*w^2 + K + Knl*NA - U^2*Ktheo;
              0, 0, Cnl, zeros(1, Nd);
              0, 0, zeros(1, Nd), Cnl ];
    bstab = [-Knl*dNA*qr; -Knl*dNA*qi; 1; 0];
    vstab = Astab \ bstab;
    
    if vstab(2) < 0
        stability(i) = 1;
    else
        stability(i) = 0;
    end
end

%% --- Plot Stability Results ---
figure('Name', 'LCO Stability Scatter', 'Color', 'w');
hold on;
for i = 1:length(A_v)
    if stability(i)
        scatter(Ulc(i)/UF, A_v(i)/delta, 'g*', 'LineWidth', 1.5)
    else
        scatter(Ulc(i)/UF, A_v(i)/delta, 'r*', 'LineWidth', 1.5)
    end
end
xlabel('$\frac{U}{U_F}$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\frac{A}{\delta}$', 'Interpreter', 'latex', 'FontSize', 16)
grid on; grid minor;

% Final continuous stability plot
istab = find(stability, 1);
figure('Name', 'LCO Branch Stability', 'Color', 'w');
hold on;
% Unstable branch: Dashed line
plot([Ulc(1), Ulc(1:istab-1)]/UF, [0, A_v(1:istab-1)]/delta, 'k--', 'LineWidth', 2)
% Stable branch: Solid line
plot(Ulc(istab:end)/UF, A_v(istab:end)/delta, 'k', 'LineWidth', 2)

xlabel('$U/U_F$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$A/\delta$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
legend('Unstable', 'Stable', 'FontSize', 12, 'Location', 'best', 'Interpreter', 'latex')
grid on; grid minor;