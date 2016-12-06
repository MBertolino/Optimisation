clear all; close all; clc;

data_input

% Quadprog to minimise std given expected return
rho = 0.1;
Aineq = [];
bineq = [];
Aeq = [r;ones(size(r))];
beq = [rho; 1];
ub = ones(size(r))';
lb = zeros(size(r))';
%lb = -ub;
f = zeros(size(r));

x0 = [];

options = optimoptions('quadprog','Algorithm','interior-point-convex');
options = optimoptions(options,'Display','iter','TolCon', 1e-9,'TolFun',1e-10);


rho = 0.01:0.01:0.2;

for ii = 1:length(rho)
    beq = [rho(ii); 1];
    [x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
    
    sigma(ii) = fval;
    
end

plot(sigma, rho)
hold on
plot(diag(H), r, 'o')
hold off