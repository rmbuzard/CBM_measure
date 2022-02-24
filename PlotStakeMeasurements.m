%%% PLOT STAKE MEASUREMENTS
% Plots the measurements from time-lapse cameras and CBM measurements
%
% INPUTS
% .txt file(s) of time-lapse camera measurements output by StakeMeasure
% excel baseline datasheet formatted correctly
%
% OUTPUTS
% erosion rate graph
%
% REQUIRED TOOLBOXES
% Financial Toolbox and Statistics and Machine Learning Toolbox
%
% Written by Richard Buzard
% February 22, 2022
%% Setup
clear all; close all;

% toolbox check
if license('test','financial_toolbox') == 0
    error('Error: missing Financial Toolbox. Go to Apps > Get More Apps to download and install toolbox.')
end
if license('test','statistics_toolbox') == 0
    error('Error: missing Statistics and Machine Learning Toolbox. Go to Apps > Get More Apps to download and install toolbox.')
end

addpath(fileparts(matlab.desktop.editor.getActiveFilename)) % add code folder to paths
[filename, filepath] = uigetfile('*.txt','Select Measurement Files', 'Multiselect','on');
%cd(filepath)
format long

% sort files to make life easier
if iscell(filename)                 % if there are multiple files
    fileSort = filename';           % get sort index
    filename = sortrows(fileSort);  % sort by filename
else 
    filename = {filename};          % if one file just put in cell array
end

% Get profile number
u = find(filename{1}=='_');
pnum = filename{1}(u(1)+2:u(2)-1);
%% Combine and calculate stats

% create beginning variables
ia  = length(filename);       % count number of files
ef = {};
time = {};
polyEF   = [];
polyTime = [];
efTemp   = [];
timeTemp = [];

for ii = 1:ia          % for each file of this transect
    YTemp    = [];
    XTemp    = [];
    fid      = fopen(fullfile(filepath,filename{ii}));     % open the file
    data     = textscan(fid,'%f%f');    % read into a matrix
    YTemp    = data{1,1};
    XTemp    = data{1,2};
    efTemp   = vertcat(efTemp,YTemp);   % add to array of eroding feature
    timeTemp = vertcat(timeTemp,XTemp); % add to array of time
end

[time,I] = sort(timeTemp);          % sort by date just in case
ef = efTemp(I); 

for qq = 1:length(time)             % turn time into days for linreg
    polyTime(qq,1) = time(qq)-min(time);
end

n    = length(time);                   
coeff= polyfit(polyTime,ef,1);   % Solve for linear fit
a    = coeff(1);                           % slope
b    = coeff(2);                           % intercept
x    = time';                            % stole from the internet
y    = ef';                              
Sxx  = sum((x-mean(x)).^2);
Syy  = sum((y-mean(y)).^2);
Sxy  = sum((x-mean(x)).*(y-mean(y)));
SSE  = Syy-a*Sxy;
S2yx = SSE/(n-2); 
Syx  = sqrt(SSE/(n-2));        % Standard error of estimate (meter)


for qq = 1:length(time)     % Get mx + b values
    mxb(qq) = a*polyTime(qq) + b;
end

% Make legend info
legendTransect = strcat(['PROFILE ',pnum]);
%legendPoly{kk} = strcat([num2str(a*365), ' +/- ', num2str(Syx(kk)*2), ' M/Y']);
    
baselineActive = 0;
%% Get Baseline measurements

baselineActive = 1;
base = baselineReader2();               % run baseline reader 2
for bb = 1:length(base)                 % get profile numbers indexed
    baseProf(bb) = base{bb}.prof;
end

bb = find(baseProf == str2double(pnum));    % index profile location
temp = base{bb};                            % pull struct data
efGT = temp.erosion + temp.erosion(1);  % get erosion distance
timeGT = temp.efDatenum;                % get date number

%% Correct to ground-truth values
% GT is in meters from start (0 to negative)
% ef is in m from stake
% efTL is in m from start
% First, identify a horizontal offset correction
% Find all GT times earlier than first TL time
startGT = find(timeGT<time(1));
if isempty(startGT)             % if GT starts later
    startGT = 1;                % set index to 1
end
efGTOffset = efGT(startGT(end));% GT offset
efOffset = ef(1) - efGTOffset;  % adjust time-lapse 
efTL = ef-efOffset;             % create TL offset variable
for ee = 2:length(efTL)         % for each erosion measurement
    efTL(ee) = min(efTL(ee), efTL(ee-1));   % prevent accretion
end

% for each GT point
for tt = 1:length(timeGT)
    % find TL points at/before GT
    tidx = find(timeGT(tt)>=time); 
    % adjust any TL points farther than GT to GT
    efTL(efTL(tidx)<efGT(tt)) = efGT(tt);
end   

% combine GT and TL
efCB = [efTL;efGT'];
timeCB = [time;timeGT];
[timeCB,I] = sort(timeCB);   % sort by date 
efCB = efCB(I); 

for ee = 2:length(efCB)         % for each erosion measurement
    efCB(ee) = min(efCB(ee), efCB(ee-1));   % prevent accretion
end
    
%% Plot it
%load('PTH_OTF.mat')
clf
C = {'b','r','g','y','k','m','c',[.5 .6 .7],[.8 .2 .6]};
figure(1)
% stuff to make each figure window consistent size!!!!
axes('Position',[0 0 1 1],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
set(gcf,'Units', 'inches','Position', [1 1 6.4 4])
hold on

% Correlate baseline datasheet with timelapse
if baselineActive == 1
    startGT = find(timeGT<time(1));
    if isempty(startGT)
        startGT = 1;
    end
    efGTOffset = efGT(startGT(end));
    efOffset = ef(1) - efGTOffset;
else
    efOffset = ef(1);%*0;       % offset numbers to show erosion from start
end
plot(time,3.28084*(ef-efOffset),'.','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.75 0.75 0.75],'MarkerSize',4);

if baselineActive == 1
    plot((datenum(timeGT)+0.5),3.28084*(efGT),'s',...
        'MarkerSize',4,'MarkerFaceColor',[0.2 0.2 0.2],'MarkerEdgeColor','k')
    plot(timeCB,efCB/.3048,'-k')
    legend({'Camera','Tape Measure','Shorline'})
end

grid on
grid minor
%datetick('x',1,'keepticks');
xmin = min([time(1),timeGT(1)]);
xmax = max([time(end),timeGT(end)]);
set(gca, 'Xlim', [xmin xmax]);
title(strcat(...
    [filename{1}(1:3),' P', pnum,' ', ...
    datestr(timeGT(1),1), ' to ',...
    datestr(time(end),1)]));

ax = gca;
ax.XTick = xmin:15:xmax;
%datetick('x',6,'keeplimits');
ax.Units = 'inches';
ax.Position = [0.25   0.25   5.5   3.5];
box on
ax.Layer = 'Top';
ax.XLim(2) = ax.XLim(2) + 1;    % add day to allow CBM symbol to show up

ylabel('Erosion (feet)')
ax.XMinorTick = 'on';               % add month ticks
minorLim = dateshift(datetime(datestr((ax.XLim))), 'start', 'month');
ax.XAxis.MinorTickValues = datenum(minorLim(1):calmonths(1):minorLim(2));
datetick('x','yyyy','keeplimits')
ax.YLim(2) = 1;
ax.YAxisLocation = 'right';

for tt = 1:length(ax.XAxis.MinorTickValues)
    text(ax.XAxis.MinorTickValues(tt),...       % x location
        ax.YLim(1) + (ax.YLim(1)*0.04)/(ax.Position(4)-ax.Position(2)),...   % y location
        num2str(month(datetime(datevec(ax.XAxis.MinorTickValues(tt))))),... % month
        'FontSize',6,'HorizontalAlignment','center')
end
fclose 'all';