clear all; clear all; clc;

% Lists with tickers:
% http://www.nasdaq.com/markets/indices/major-indices.aspx
% http://bigcharts.marketwatch.com/markets/indexes.asp

% Time period
start = '2016-01-01';
stop = today;

% Predictors: ticker and full name
tickers = ...%[{'BOUL', 'Boule Diagnostics'}; ...
    [{'HTRO', 'Hexatronic Group'}; ...
    {'ATEL-A', 'AllTele A'}; ...
    {'POOL-B', 'Poolia B'}; ...
    {'BELE', 'Beijer Electronics AB'}; ...
    {'ADDT-B', 'Addtech AB'}; ...
    {'VOLV-B', 'Volvo AB'};
    {'ABB', 'ABB Ltd'};
    {'SKF-B', 'AB SKF'};
    {'ERIC-B', 'LM Ericsson'}];

% Load data
rawData = getGoogleDailyData(tickers(1:end/2), ...
    datenum(start), datenum(stop));

% Save only the dates (col. 1) and the closing prices (col. 2) into 'data'
data = struct;
markets = fieldnames(rawData);
nMarkets = length(markets);
figure(1)
hold on
for i = 1:nMarkets
    data.(markets{i}).Date = rawData.(markets{i}).Date;
    data.(markets{i}).Close = rawData.(markets{i}).Close;
    
    plot(data.(markets{i}).Close)
end


% Adjust the dates so that all markets has the same date vector
% [data, N] = correctDates(data);
% dates = data.(markets{1}).Date;





% %%
% for ii = 1:length(rho)
%     beq = [rho(ii); 1];
%     [x,fval,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
%     
%     sigma(ii) = fval;
%     
% end