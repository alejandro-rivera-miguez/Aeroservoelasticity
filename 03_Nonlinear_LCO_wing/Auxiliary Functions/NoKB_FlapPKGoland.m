%% =======================================================================
% FUNCTION: NoKB_FlapPKGoland
% DESCRIPTION: Adapts the PK formulation to the case with null actuator 
%              stiffness (KBeta = 0). Provides initial conditions for LCO.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
% =======================================================================
function [UF, wF, qF] = NoKB_FlapPKGoland(Ne, Nd, Umax, Nu, Display)
    L=6.096; c=1.829; cf=0.45;
    b=c/2; ehat=(0.33*c-b)/b; chat=(c-cf-b)/b;
    rho=1.225; Clalpha=2*pi; yL=3.35; yU=3.35+1.8;
    
    mu=acos(chat);
    T1=-1/3*sqrt(1-chat^2)*(2+chat^2)+chat*mu;
    T3=-(1/8+chat^2)*mu^2+1/4*chat*sqrt(1-chat^2)*mu*(7+2*chat^2)-1/8*(1-chat^2)*(5*chat^2+4);
    T4=-mu+chat*sqrt(1-chat^2);
    T5=-(1-chat^2)-mu^2+2*chat*sqrt(1-chat^2)*mu;
    T7=-(1/8+chat^2)*mu+1/8*chat*sqrt(1-chat^2)*(7+2*chat^2);
    T8=-1/3*sqrt(1-chat^2)*(2*chat^2+1)+chat*mu;
    T9=1/2*(1/3*sqrt(1-chat^2)^3+ehat*T4);
    T10=sqrt(1-chat^2)+mu;
    T11=mu*(1-2*chat)+sqrt(1-chat^2)*(2-chat);
    T12=sqrt(1-chat^2)*(2+chat)-mu*(2*chat+1);
    T13=-1/2*(T7+(chat-ehat)*T1);
    
    [fr, qr, FEMdata] = NoKB_FEM_FlapCoupledBendingTorsion(Ne, Nd);
    
    modes_cell = cell(1, Nd);
    for i = 1:Nd
        modes_cell{i} = @(y)[EvaluateCoupledBendingDisplacement(y, qr(:,i), FEMdata.deltay, FEMdata.Connect);...
                             EvaluateCoupledTorsionDisplacement(y, qr(:,i), FEMdata.deltay, FEMdata.Connect); Betafun(y, qr(:,i))];
    end
    
    MM = eye(Nd, Nd);
    KK = diag(fr.^2);
    
    Ca = rho*b^2*[0,pi,-T4; 0,-pi*b*(1/2-ehat),-b*(T1-T8-(chat-ehat)*T4+T11/2); 0,-b*(-2*T9-T1+T4*(ehat-1/2)),b*T4*T11/2/pi];
    Ka = rho*b^2*[0,0,0; 0,0,-(T4+T10); 0,0,-(T5-T4*T10)/pi];
    Ma = rho*b^2*[pi,-pi*b*ehat,-b*T1; pi*b*ehat,-pi*b^2*(1/8+ehat^2),b^2*(T7+(chat-ehat)*T1); b*T1,-2*b^2*T13,b^2*T3/pi];
    Bw = rho*b*Clalpha*[1;b*(ehat+1/2);-T12*b/2/pi];
    Cw = [0,1,T10/pi];
    Cwhat = [1,b*(1/2-ehat),b*T11/2/pi];
    
    Kstar1 = nan(Nd,Nd); Kstar2 = nan(Nd,Nd); Kstar3 = nan(Nd,Nd); Kstar4 = nan(Nd,Nd);
    Cstar1 = nan(Nd,Nd); Cstar2 = nan(Nd,Nd); Cstar3 = nan(Nd,Nd);
    
    for i = 1:Nd
        for j = 1:Nd
            Kstar1(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ma), 0, L, 'Waypoints', [yL,yU]);
            Kstar2(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cw), 0, L, 'Waypoints', [yL,yU]);
            Kstar3(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cwhat), 0, L, 'Waypoints', [yL,yU]);
            Kstar4(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ka), 0, L, 'Waypoints', [yL,yU]);
            
            Cstar1(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Ca), 0, L, 'Waypoints', [yL,yU]);
            Cstar2(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cw), 0, L, 'Waypoints', [yL,yU]);
            Cstar3(i,j) = integral(@(y)AeroIntegrand(y, modes_cell{i}, modes_cell{j}, Bw*Cwhat), 0, L, 'Waypoints', [yL,yU]);
        end
    end
    
    UU = linspace(0, Umax, Nu); 
    Eigenvalues = nan(Nd, length(UU));
    qq = nan(Nd, Nd, length(UU));
    Eigenvalues(:,1) = fr * 1i;
    qq(:,:,1) = eye(Nd, Nd); 
    
    tol = 1e-6;   
    maxiter = 20; 
    UF = nan; 
    
    for i = 2:length(UU)
        U = UU(i);
        Um1 = UU(i-1);
        for j = 1:Nd
            lambda = Eigenvalues(j, i-1);
            q = qq(:, j, i-1);
            if Um1 == 0 
                Ctheo = zeros(Nd, Nd);
                Ktheo = zeros(Nd, Nd);
            else
                k = imag(lambda)*b/Um1;
                if k ~= 0
                    Ck = besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
                    Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                    Ctheo = Cstar1 + b/k*imag(Ck)*Cstar2 + real(Ck)*Cstar3;
                else
                    Ck = 1;
                    Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                    Ctheo = Cstar1 + real(Ck)*Cstar3;
                end
            end
            
            ddu = [2*MM*lambda*q - Um1*Ctheo*q, MM*lambda^2 - Um1*Ctheo*lambda + (KK - Um1^2*Ktheo); 0, 2*q'] \ [Ctheo*lambda*q - (KK - 2*Um1*Ktheo)*q; 0];
            lambda = lambda + ddu(1)*(U - Um1);
            q = q + ddu(2:end)*(U - Um1);
            
            k = imag(lambda)*b/U;
            if k ~= 0
                Ck = besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
                Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                Ctheo = Cstar1 + b/k*imag(Ck)*Cstar2 + real(Ck)*Cstar3;
            else
                Ck = 1;
                Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                Ctheo = Cstar1 + real(Ck)*Cstar3;
            end
            error = norm([(MM*lambda^2 - U*Ctheo*lambda + (KK - U^2*Ktheo))*q; 1 - q'*q]);
            
            iter = 0;
            while error > tol && iter <= maxiter
                A = [2*lambda*MM*q - U*Ctheo*q, MM*lambda^2 - U*Ctheo*lambda + (KK - U^2*Ktheo); 0, 2*q'];
                v = [-MM*lambda^2*q + U*Ctheo*lambda*q - (KK - U^2*Ktheo)*q; 1 - q'*q];
                delta = A\v;
                
                lambda = lambda + delta(1);
                q = q + delta(2:end);
                
                k = imag(lambda)*b/U;
                if k ~= 0
                    Ck = besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
                    Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                    Ctheo = Cstar1 + b/k*imag(Ck)*Cstar2 + real(Ck)*Cstar3;
                else
                    Ck = 1;
                    Ktheo = -(k/b)^2*Kstar1 + real(Ck)*Kstar2 - k/b*imag(Ck)*Kstar3 + Kstar4;
                    Ctheo = Cstar1 + real(Ck)*Cstar3;
                end            
                iter = iter + 1;
                error = norm([(MM*lambda^2 - U*Ctheo*lambda + (KK - U^2*Ktheo))*q; 1 - q'*q]);
            end
            Eigenvalues(j,i) = lambda;
            qq(:,j,i) = q;
        end
        if max(real(Eigenvalues(:,i))) > 0 
            [~, I] = max(real(Eigenvalues(:,i)));
            gammaim1 = real(Eigenvalues(I, i-1));
            gammai = real(Eigenvalues(I, i));
            UF = Um1 - gammaim1/(gammai - gammaim1)*(U - Um1);
            break
        end
    end
    
    iF = i;
    if Display
        Frequencies = imag(Eigenvalues)/2/pi; 
        figure 
        for i = 1:Nd
            plot(UU, Frequencies(i,:), 'DisplayName', sprintf('Mode %d', i))
            hold on
        end
        grid minor; legend; title('Modal Frequencies vs Airspeed');
        xlabel('$U \; [m/s]$', 'Interpreter', 'latex'); ylabel('$f \; [Hz]$', 'Interpreter', 'latex');
        
        figure
        for i = 1:Nd
            plot(UU, real(Eigenvalues(i,:)), 'DisplayName', sprintf('Mode %d', i))
            hold on
        end
        grid minor; legend; title('Modal real part vs Airspeed');
        xlabel('$U \; [m/s]$', 'Interpreter', 'latex'); ylabel('$\Re(\lambda) \; [1/s]$', 'Interpreter', 'latex');
    end
    
    % Interpolate qF and wF for return
    qF = qq(:, 2, iF) - gammaim1/(gammai - gammaim1)*(qq(:, 2, iF) - qq(:, 2, iF-1));
    Cnl = qr(end, :);
    Aq = [Cnl*real(qF), -Cnl*imag(qF); Cnl*imag(qF), Cnl*real(qF)];
    ab = Aq \ [1; 0];
    qF = (ab(1) + 1i*ab(2)) * qF;
    wF = imag(Eigenvalues(2, iF)) - gammaim1/(gammai - gammaim1)*(imag(Eigenvalues(2, iF)) - imag(Eigenvalues(2, iF-1)));
end