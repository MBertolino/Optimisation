%% Analysis
clear all, close all
load('Gdata2')
% timeframe for rate of return
TF = 30;

figure(1)
hold on
for ii = 1:nMarkets
    data.(markets{ii}).index = data.(markets{ii}).Close/data.(markets{ii}).Close(1);
    plot(data.(markets{ii}).index)
    
    
    R(:, ii) = (data.(markets{ii}).Close(TF:end-TF)-data.(markets{ii}).Close(1:end-2*TF+1))./data.(markets{ii}).Close(1:end-2*TF+1);
    %R(:, ii) = data.(markets{ii}).index(1:end);
end

r = mean(R)
H = cov(R)
stdev = sqrt(diag(H));


% hold on
% for m = 1:100
%     v = 1;
%     for kk = 2:173
%         v(kk) = v(kk-1)*random('Normal',1+r(1),stdev(1),1,1)^(1/TF);
%     end
% 
%     plot(v)
% end


%% 4. Efficient frontier, stepping in alpha
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

options = optimoptions('quadprog','Algorithm','interior-point-convex');
options = optimoptions(options,'Display','iter','TolCon', 1e-9,'TolFun',1e-10);

Aeq = ones(size(r));
beq = 1;

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

%% Plot efficient frontier
figure(2)
subplot(2,1,1)
plot(stdev,r,'o')
hold on;
plot(sigma_ef,ret_ef)
xlabel('Volatility')
ylabel('Expected return')
legend('Individual assets', 'Efficient frontier')
xlim([0 max(stdev)*1.1])
ylim([min(r) max(r)*1.2])

subplot(2,1,2)
area(sigma_ef,w4)
xlim([sigma_ef(end) sigma_ef(1)])
ylim([0 1])
xlabel('Volatility')
ylabel('Shares in asset')
%legend('1','2','3','4','5')


%% Validate on last part of set
for kk = 1:length(sigma_ef)
    for ii = 1:nMarkets
        GT(ii) = (data.(markets{ii}).Close(end) - data.(markets{ii}).Close(end-TF))/data.(markets{ii}).Close(end-TF); 
    end
    
    ROR(kk) = w4(kk,:)*GT';
end

figure(3)
plot(sigma_ef, ROR)
ylabel('Actual return on investment')
xlabel('Volatility')