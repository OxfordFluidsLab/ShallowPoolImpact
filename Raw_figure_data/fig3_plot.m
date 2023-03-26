% fig4_plot.m
%
% Script that can be used to gain a basic plot for the raw data underlying
% figure 4 (fig4.csv)
%
% Source:  https://github.com/OxfordFluidsLab/ShallowPoolImpact
% Licence: GPL-3.0 (see the Git repo).
%
% T.C. Sykes (t.c.sykes@outlook.com)
% University of Oxford (2022)

% Read in the data from a csv file
t = readtable('fig3.csv');

% Setup figure
clf
hold on

% Loop over each row in the csv file
for ii = 1:size(t,1)
    
    % Determine the desired marker colour from the outcome given
    switch t.outcome{ii}
        case 'lamella'
            col = 'r';
        case 'separate'
            col = 'g';
        case 'axisymmetric'
            col = 'y';
        otherwise
            continue
    end
    
    % Determine the desired marker type from whether satellites are present
    switch t.satellites(ii)
        case true
            mkr = 'o';
        case false
            mkr = 'v';
        otherwise
            continue
    end
    
    % Plot the point with the given marker and colour
    scatter(t{ii,1},t{ii,2},20,col, 'filled', 'Marker',mkr);
    
end

% Figure formatting
box on
set(gca,'xscale','log');            % Select log scale for the 
xlim([min(t.hstar),max(t.hstar)])   % Limit x axis to the data
ylim([min(t.We),max(t.We)])         % Limit y axis to the data
xlabel('hstar')                     % Label x axis
ylabel('We')                        % Label y axis
