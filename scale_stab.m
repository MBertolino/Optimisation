clear all; close all; clc;
load('Sharad.mat')

%%
% Initial guesses 
Ai = 1e-26;
ni = 2.1;
p0 = [Ai ni];

Ascale = 1e-26;

% Optimization options
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton','OptimalityTolerance', 1e-9, 'Display', 'iter', 'MaxFunctionEvaluations', 1024);


% try differenet Ascales
A_vec = logspace(-23,-27, 128);

for(ii = 1:128)
   Ascale = A_vec(ii);
   p0(1) = Ai / Ascale;
   

   % Miniminze H_cost3 to find both A and n simultaniously
    H_cost3 = @(p) norm((-a.*(p(2)+2)./(2*p(1)*Ascale).*(rho*g).^(-p(2)).*abs(dhdx).^(1-p(2))./dhdx).^(1/(p(2)+2)) - H_obs)^2;
    
    [p,f3val,exitflag3,output3] = fminunc(H_cost3,p0,options);
    
    Z(ii) = f3val;
end


semilogx(A_vec, Z)
xlabel('Ascale')
ylabel('Fval')
