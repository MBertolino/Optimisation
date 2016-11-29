clear all; close all; clc;

% Load data
data_input

% Quadprog to minimise std given expected return
rho = 0.1;
Aineq = [];
bineq = [];
Aeq = [r;ones(size(r))];
beq = [rho; 1];
ub = ones(size(r))'*0.5;
lb = zeros(size(r))';
lb = -ub;
f = zeros(size(r));
x0 = [];

options = optimoptions('quadprog','Algorithm','interior-point-convex');
options = optimoptions(options,'Display','iter','TolCon', 1e-9,'TolFun',1e-10);

% Efficient frontier, stepping in expected return
rho = 0.01:0.002:0.2;

sigma = zeros(size(r));
for ii = 1:length(rho)
    beq = [rho(ii); 1];
    [x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
    
    sigma(ii) = fval;
    w(ii,:) = x;
end

figure()
plot(sigma, rho) 
hold on;
plot(diag(H),r,'o')
ylabel('Expected return')
xlabel('Volatility')

figure()
area(sigma,w)
xlabel('Volatility')

