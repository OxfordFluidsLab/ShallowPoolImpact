% fig8a_plot.m
%
% Script that can be used to gain a basic plot for the raw data underlying
% figure 7b (fig7b.csv)
%
% Source:  https://github.com/OxfordFluidsLab/ShallowPoolImpact
% Licence: GPL-3.0 (see the Git repo).
%
% T.C. Sykes (t.c.sykes@outlook.com)
% University of Oxford (2022)

% Setup figure
clf
hold on

% Initialise a cell array to hold the legend names
lgdNames = {};


% Loop over each csv file
for ii = {dir('fig8b_*.csv').name}

    % Read in the data from the csv file
    t = readtable(ii{:});
    
    % Plot the series
    scatter(t{:,1},t{:,2},20, 'filled', 'Marker','o');

    % Get legend name
    lgdName         = extractBetween(ii{:},'fig8b_','.csv');
    lgdNames{end+1} = lgdName{:};
    
end

% Figure formatting
box on
xlim([0 0.54])                      % Limit x axis to the data
ylim([0 0.045])                     % Limit y axis to the data
xlabel('Dimensionless time')        % Label x axis
ylabel('Ejecta base thickness')     % Label y axis

% Legend
legend(lgdNames, 'Location','northwest')

% Housekeeping
clearvars ii lgdName* t