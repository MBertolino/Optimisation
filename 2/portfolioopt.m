clear all; clear all; clc;

% Lists with tickers:
% http://www.nasdaq.com/assets/indices/major-indices.aspx
% http://bigcharts.marketwatch.com/assets/indexes.asp

% Time period
start = '2015-01-01';
stop = '2016-11-26';

% Predictors: ticker and full name
tickers = ...%[{'BOUL', 'Boule Diagnostics'}; ...
    [{'SSAB-B', 'SSAB B'}; ...
    {'ATEL-A', 'AllTele A'}; ...
    {'POOL-B', 'Poolia B'}; ...
    {'PREV-B', 'Prevas B'}];

% Load data
rawData = getGoogleDailyData(tickers(1:end/2), ...
    datenum(start), datenum(stop));

% Save only the dates (col. 1) and the closing prices (col. 2) into 'data'
data = struct;
assets = fieldnames(rawData);
nAssets = length(assets);
stocks = [];
for i = 1:nAssets
    data.(assets{i}).Date = rawData.(assets{i}).Date;
    data.(assets{i}).Close = rawData.(assets{i}).Close;
    
    % Save to matrix
    stocks = [stocks data.(assets{i}).Close(2:end) - ...
        data.(assets{i}).Close(1:end-1)];
end

r = mean(stocks);
stdev = std(stocks);
Q = cov(stocks);

%-------------------------- 2 Quadprog -----------------------------
%%

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
options = optimoptions(options,'Display','iter','TolFun',1e-10);
    [x,fval,exitflag] = quadprog(Q,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);


rho = 0.01:0.01:0.2;
for ii = 1:length(rho)
    beq = [rho(ii); 1];
    [x,fval,exitflag] = quadprog(Q,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
    
    sigma(ii) = fval;
    
end

figure(1)
plot(sigma, rho)


















