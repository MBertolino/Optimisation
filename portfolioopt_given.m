clear all; close all; clc;

data_input

% Quadprog to minimise std given expected return
rho = 0.1;
Aineq = [];
bineq = [];
% Aeq = blkdiag(r,ones(size(r)))';
Aeq = [r; ones(size(r))];
beq = [rho; 1];
ub = ones(size(r));
lb = zeros(size(r));
x0 = [];
f = zeros(size(r));

options = optimset('quadprog')
options = optimset(options,'Display','iter')

[x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);