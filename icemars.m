clear all; close all; clc;
load('Sharad.mat')

% Initial guesses 
Ai = 1e-25;
ni = 1.1;
p0 = [Ai ni];

% Optimization options
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton');

% Objective function
H_cost = @(n) norm((-a.*(n+2)./(2*Ai).*(rho*g).^(-n).*abs(dhdx).^(1-n)./dhdx).^(1/(n+2)) - H_obs)^2;

% Minimize H_cost to find n
[n,Hval,exitflag1,output1] = fminunc(H_cost,ni,options);

% Minimize H_cost2 to find A
H_cost2 = @(A) norm((-a.*(n+2)./(2*A).*(rho*g).^(-n).*abs(dhdx).^(1-n)./dhdx).^(1/(n+2)) - H_obs)^2;
[A,Aval,exitflag2,output2] = fminunc(H_cost2,Ai,options);

% Miniminze H_cost3 to find both A and n simultaniously
H_cost3 = @(p) norm((-a.*(p(2)+2)./(2*p(1)).*(rho*g).^(-p(2)).*abs(dhdx).^(1-p(2))./dhdx).^(1/(p(2)+2)) - H_obs)^2;
[p,fval,exitflag3,output3] = fminunc(H_cost3,p0,options);


n
A
p

A_plot = linspace(1.089e-25,1.1e-25);
n_plot = linspace(2.2,2.6);
for i = 1:length(n_plot)
    for j = 1:length(A_plot)
        H_plot(i,j) = H_cost3([A_plot(i) n_plot(j)]);
    end
%     H_plot(i) = H_cost2(A_plot(i));
end

% plot(1:100, H_plot)
surf(A_plot,n_plot,H_plot)