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
    {'DOM', 'Domestic Group'}];

% Load data
rawData = getGoogleDailyData(tickers(1:end/2), ...
    datenum(start), datenum(stop));

% Save only the dates (col. 1) and the closing prices (col. 2) into 'data'
data = struct;
markets = fieldnames(rawData);
nMarkets = length(markets);
for i = 1:nMarkets
    data.(markets{i}).Date = rawData.(markets{i}).Date;
    data.(markets{i}).Close = rawData.(markets{i}).Close;
end

% Adjust the dates so that all markets has the same date vector
% [data, N] = correctDates(data);
% dates = data.(markets{1}).Date;
