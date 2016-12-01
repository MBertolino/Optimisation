clear all; close all; clc;

data_input

% Quadprog to minimise std given expected return
rho = 0.1; % Seeked return
Aineq = []; % Olikheter, 0 för oss
bineq = []; % Olikheter, 0 för oss
Aeq = [r;ones(size(r))]; % likhets constraint
beq = [rho; 1]; % likhets constraints
ub = ones(size(r))'; % Upper bound blir 1 1 1 1 1 - Högst vikt i en aktie i %
lb = zeros(size(r))'; % Lower bound blir 0 eller -1, där 0 är ingen aktie, -1 är shortsell av aktie i %
%lb = -ub;
f = zeros(size(r)); % Nothing

x0 = []; % Startvärde? - Interior point convex algoritm slumpar startvärde.

options = optimoptions('quadprog','Algorithm','interior-point-convex');
options = optimoptions(options,'Display','iter','TolCon', 1e-9,'TolFun',1e-10);
% Display iter - visar itterationen medan den räknar ut f(x)
% TolCon - Tolerans Condition
% TolFun - Tolerans Funktion
% quadprog - 
% Algorithm - 


rho = 0.01:0.01:0.2; % Different rate of return conditions

for ii = 1:length(rho)
    beq = [rho(ii); 1];
    [x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
    
    sigma(ii) = fval;
    
end

figure(1)
plot(sigma.^2, rho)


alpha = 0.05:0.05:1;
rho = 0.01:0.01:0.2;

for ii = 1:length(rho)
    for kk = 1:length(alpha)
        beq = [rho(ii); 1];
        aH = alpha(kk)*H;
        f = (1-alpha(kk))*r;
        [x,fval,exitflag] = quadprog(aH,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
        
        sigma_ra(ii,kk) = fval;
    end
end

surf(alpha, rho, sigma_ra.^2)



