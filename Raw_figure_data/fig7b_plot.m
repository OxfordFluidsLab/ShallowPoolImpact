% fig7b_plot.m
%
% Script that can be used to gain a basic plot for the raw data underlying
% figure 7b (fig7b.csv)
%
% Source:  https://github.com/OxfordFluidsLab/ShallowPoolImpact
% Licence: GPL-3.0 (see the Git repo).
%
% T.C. Sykes (t.c.sykes@outlook.com)
% University of Oxford (2022)

% Read in the data from a csv file
t = readtable('fig7b.csv');

% Setup figure
clf
hold on

% Loop over each column in the csv file
for ii = 2:size(t,2)
    
    % Plot the series
    scatter(t{:,1},t{:,ii},20, 'filled', 'Marker','o');
    
end

% Figure formatting
box on
xlim([0 0.054])                     % Limit x axis to the data
ylim([-0.1 1.45])                   % Limit y axis to the data
xlabel('Dimensionless time')        % Label x axis
ylabel('Ejecta sheet tip velocity') % Label y axis

% Legend
legend(t.Properties.VariableNames{2:end}, 'Location','northwest')

% Housekeeping
clearvars ii t