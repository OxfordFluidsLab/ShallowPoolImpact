% ConfocalAnalysisAndPlot.m
%
% This script analyses confocal data to determine pool depth (printed to
% the terminal). It also plots the pool free surface and base position time
% series on a single graph, with the constructions (fits, change points,
% etc.) that were used to determine the reported pool depth.
%
% Confocal data should be provided as a .csv file (in the current working
% directory) named "cPool_#.csv", where # is a four digit integer
% identifying the particular experiment (expNo). Each measurement must have
% its own row. There should be 4 columns:
% 1) MATLAB datetime string
% 2) Elapsed measurement time (seconds) 
% 3) Free surface position (relative to the sensor, millimeters)
% 4) Pool base position (relative to the sensor, millimeters)
%
% IMPORTANT:
% Remember to select the appropriate material (here, ethanol) in the
% confocal UI!
%
% We generally took confocal measurements at 2 Hz.
% 
% A second .csv file is required (e2eTime.csv), giving the time between the
% last confocal measurement and event of interest (i.e. impact at t=0). The
% data was automatically evaluated by the MATLAB script used to control the
% experiment (including initiating confocal measurements), using an Arduino
% to detect when the high-speed camera was triggered, which was used to
% approximate the event time (i.e. t=0). Each measurement must have its own
% row. There should be 2 columns:
% 1) expNo
% 2) Time between the last confocal measurement and event (seconds).
%
% IMPORTANT:
% This script is distributed with the horizontal axis relabelled to time,
% in order to match the graph in the supplementary of the paper. Please see
% the lines below "% Change to plot against time". It is suggested to
% comment out the associated lines if changing the input data.
%
% Source:  github.com/OxfordFluidsLab/ShallowPoolImpact
% Licence: GPL-3.0 (see LICENCE in the root of the Git repo)
%
% T.C. Sykes (t.c.sykes@outlook.com)
% University of Oxford (2022)

% Change this number to choose different input.
expNo = 225;

% Hard-coded parameters
bThresh   = 50;    % Threshold for base positon ischange (variance)
fThresh   = 1e-3;  % Threshold for free surface position ischange (linear)
txtSize   = 8.5;
mkrSize   = 18;

% Initialise colours (https://personal.sron.nl/~pault/#sec:qualitative)
muted        = struct();
muted.green  = [17,119,51]./255;
muted.indigo = [51,34,136]./255;
muted.rose   = [204,102,119]./255;
muted.purple = [170,68,153]./255;
muted.sand   = [221,204,119]./255;

% Read in the corresponding confocal data
cDataT = readtable(strcat('.\cPool_',num2str(expNo,'%04d'),'.csv'));

% Clear a figure for plotting
clf
hold on


%% DETERMINE THE BASE POSITION ON IMPACT

% Finds abrupt changes in the variance of the base position, if any
bChangeLogical = ...
    ischange(cDataT.Var4,'variance','Threshold',bThresh);
bChangeIndices = find(bChangeLogical==true);

% Get a vector of (fixed) substrate positions
if isempty(bChangeIndices)
    bHeights = cDataT.Var4;
else
    bHeights = cDataT.Var4(bChangeIndices(end):end);
end

% Remove NaN entries; take the median as THE substrate position
bHeights = bHeights(~isnan(bHeights));
bHeight  = median(bHeights);

% Housekeeping
clearvars bHeights bThresh bChangeLogical


%% DETERMINE THE POOL FREE SURFACE POSITION ON IMPACT

% Find abrupt changes in the slope and intercept of the pool free surface
% position, assuming that the time series consists of many linear parts.
% Setting a large detection threshold reduces the number of change points
% detected (i.e. reduces the effect of noise/random variations).
[fChangeLogical,f1,f2] = ...
    ischange(cDataT.Var3,'linear','Threshold',fThresh);
fChangeIndices         = find(fChangeLogical==true);

% Restrict the table to the final linear regime of free surface positions
if isempty(fChangeIndices)
    fHeightsT = cDataT(:,1:3);
else
    fHeightsT = cDataT(fChangeIndices(end):end,1:3);
end

% Remove NaN entries from these vectors
fHeightsT = rmmissing(fHeightsT,'DataVariables',{'Var3'});

% Get the time delay between the last confocal measurement and trigger
e2eT    = readtable('.\e2eTime.csv');
e2eloc  = find(e2eT.Var1==expNo,1,'last');
e2eTime = e2eT.Var2(e2eloc);

% Construct the time vector
unShiftedTimes = fHeightsT.Var2-fHeightsT.Var2(end);
timeToTrigger  = unShiftedTimes - e2eTime;

% Fit for the final linear regime for pool positions.
lFit = polyfit(timeToTrigger,fHeightsT{:,3},1);

% Get final free surface position
fHeight = lFit(end);

% Housekeeping
clearvars e2eloc e2eT e2eTime expNo fThresh timeToTrigger


%% Plotting

% Underlying data (base)
plot(1:length(cDataT.Var4),cDataT.Var4,'k.', ...
    'MarkerSize',mkrSize);
% Change locations (base)
if ~isempty(bChangeIndices)
    xline(bChangeIndices,'--', 'LineWidth',2, 'Color',muted.green);
end
% Final position (base)
yline(bHeight,'-','LineWidth',2,'Color',muted.green);

% Linear Segments (free surface)
segline = f1.*transpose(1:length(cDataT.Var3)) + f2;
plot(1:length(cDataT.Var3),segline,'-', 'LineWidth',1.5, ...
    'Color',muted.indigo)
% Underlying Data (free surface)
plot(1:length(cDataT.Var3),cDataT.Var3,'.', ...
    'MarkerSize',mkrSize, 'Color',muted.rose)
% Final position (free surface)
yline(fHeight,'-', 'LineWidth',2,'Color', muted.purple);
% Linear fit to the final part
cRate = length(unShiftedTimes)/max(abs(unShiftedTimes));
mi    = lFit(1)/cRate;
if isempty(fChangeIndices)
    xp = 1:length(fChangeLogical);
else
    xp = fChangeIndices(end):length(fChangeLogical);
end
yp = fHeightsT{:,3}';
c  = polyfit(xp,yp - mi*xp,0);
xz = ceil( (fHeight - c)/mi );
plot([xp(1),xz],mi*[xp(1),xz]+c,'-', 'LineWidth',2, 'Color',muted.sand);

% Legend icons
% plot(92,5.37,'.','MarkerSize',markerSize, 'Color', muted.rose)
% plot(92,5.42,'k.','MarkerSize',markerSize);

% Housekeeping
clearvars bChangeIndices c cDataT cRate f1 f2 fChange* fHeightsT ...
    lFit mkrSize muted mi segline unShiftedTimes xp xz yp

% Change to plot against time
xticks([0 40 80 120 160])
ax = gca;
set(ax,'XTickLabel',num2str((-80:20:0)'))
% Explanation: To simplify our analysis, we generally plotted against time
% indices (sequential numbers corresponding to each confocal data point)
% However, for the convenience of the reader, the graph in the
% supplementary material is converted to plot against seconds. For ease of
% plotting, and retain the original characteristics of this script, we have
% just relablled the horizontal axis for the supplementary material.
% Simply comment out the 3 lines above to see our usual working graph (with
% the horizontal axis in terms of time indices).

% Figure formatting
try
    % Using figure formatting functions, if available, to match the style
    % of the graph in the supplementary material
    ax = gca;
    setupxaxis('Time, $t$ (s)', ax.XLim, txtSize);
    setupyaxis('Relative displacement (mm)', ax.YLim, txtSize);
    setfiguresize([10 6],'centimeters')
    basicfiguresetup(txtSize, 2, '');
    expandfigureaxes();
    saveasbitmap(strcat(pwd,'/','confocalExample'),'png',true);
catch
    % Label the axes
    box on
    xlabel('Time, t (s)')
    ylabel('Relative displacement (mm)')
end

% Print the pool height to the terminal
depth = round(bHeight - fHeight,2);
txt = sprintf('The pool depth should be reported as %0.2f mm (2dp)',depth);
disp(txt)

% Housekeeping
clearvars ax txtSize txt
