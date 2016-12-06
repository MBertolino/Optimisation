clear all; close all; clc;

% Load data
data_input

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> MBertolino/master
=======
>>>>>>> MBertolino/master
%% 2. KKT system. Risk avert
A11 = H;
A12 = ones(length(r),1);
A21 = A12';
A22 = 0;

A = [A11 A12; A21 A22];
b = [zeros(size(r))'; 1];

% The first 5 are the weights
w2 = A\b;
w2 = w2(1:length(r));

%% 3. Quadprog to minimise std given expected return
% b) rho_b = 0.1
rho_b = 0.1;
Aeq_b = [r; ones(size(r))];
beq_b = [rho_b; 1];

% c) rho_c = 0.2
rho_c = 0.2;
Aeq_c = [r; ones(size(r))];
beq_c = [rho_c; 1];

ub = ones(size(r))';
lb = zeros(size(r))';
%lb = -ub;
f = zeros(size(r));
x0 = [];
>>>>>>> MBertolino/master

options = optimoptions('quadprog','Algorithm','interior-point-convex');
options = optimoptions(options,'Display','iter','TolCon', 1e-9,'TolFun',1e-10);
% Display iter - visar itterationen medan den räknar ut f(x)
% TolCon - Tolerans Condition
% TolFun - Tolerans Funktion
% quadprog - 
% Algorithm - 

% b)
w3b = quadprog(H,f,[],[],Aeq_b,beq_b,lb,ub,x0,options);
sigma2_b = w3b'*H*w3b;

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
rho = 0.01:0.01:0.2; % Different rate of return conditions

for ii = 1:length(rho)
    beq = [rho(ii); 1];
    [x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
    
    sigma(ii) = fval;
    
end

figure(1)
plot(sigma.^2, rho)
=======
% c)
w3c = quadprog(H,f,[],[],Aeq_c,beq_c,lb,ub,x0,options);
sigma2_c = w3c'*H*w3c;
>>>>>>> MBertolino/master
=======
% c)
w3c = quadprog(H,f,[],[],Aeq_c,beq_c,lb,ub,x0,options);
sigma2_c = w3c'*H*w3c;
>>>>>>> MBertolino/master
=======
% c)
w3c = quadprog(H,f,[],[],Aeq_c,beq_c,lb,ub,x0,options);
sigma2_c = w3c'*H*w3c;
>>>>>>> MBertolino/master

%% 4. Efficient frontier, stepping in alpha
Aeq = ones(size(r));
beq = 1;

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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


=======
alpha = 0:0.05:1;
w4 = zeros(length(alpha),length(r));
sigma_ef = zeros(length(alpha),1);
ret_ef = sigma_ef;
for ai = 1:length(alpha)
    Q = 2*alpha(ai)*H;
    f = (alpha(ai)-1)*r';
    [w4(ai,:),fval,exitflag] = quadprog(Q,f,[],[],Aeq,beq,lb,ub,[],options);
    sigma_ef(ai) = sqrt(w4(ai,:)*H*w4(ai,:)');
    ret_ef(ai) = r*w4(ai,:)';
end

=======
alpha = 0:0.05:1;
w4 = zeros(length(alpha),length(r));
sigma_ef = zeros(length(alpha),1);
ret_ef = sigma_ef;
for ai = 1:length(alpha)
    Q = 2*alpha(ai)*H;
    f = (alpha(ai)-1)*r';
    [w4(ai,:),fval,exitflag] = quadprog(Q,f,[],[],Aeq,beq,lb,ub,[],options);
    sigma_ef(ai) = sqrt(w4(ai,:)*H*w4(ai,:)');
    ret_ef(ai) = r*w4(ai,:)';
end

>>>>>>> MBertolino/master
=======
alpha = 0:0.05:1;
w4 = zeros(length(alpha),length(r));
sigma_ef = zeros(length(alpha),1);
ret_ef = sigma_ef;
for ai = 1:length(alpha)
    Q = 2*alpha(ai)*H;
    f = (alpha(ai)-1)*r';
    [w4(ai,:),fval,exitflag] = quadprog(Q,f,[],[],Aeq,beq,lb,ub,[],options);
    sigma_ef(ai) = sqrt(w4(ai,:)*H*w4(ai,:)');
    ret_ef(ai) = r*w4(ai,:)';
end

>>>>>>> MBertolino/master
% Plot efficient frontier
figure()
subplot(2,1,1)
plot(stdev,r,'o')
hold on;
plot(sigma_ef,ret_ef)
xlabel('Volatility')
ylabel('Expected return')
legend('Individual assets', 'Efficient frontier')
xlim([0.05 max(stdev)*1.1])
ylim([0 max(r)*1.2])
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> MBertolino/master
=======
>>>>>>> MBertolino/master
=======
>>>>>>> MBertolino/master

subplot(2,1,2)
area(sigma_ef,w4)
xlim([sigma_ef(end) sigma_ef(1)])
ylim([0 1])
xlabel('Volatility')
ylabel('Shares in asset')
legend('1','2','3','4','5')