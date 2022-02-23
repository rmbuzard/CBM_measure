function [out] = baselineReader2(~)

%%% BASELINEREADER2 reads baseline data file for stakes, outputs transect
%%% information
%
% INPUT
% Baseline datasheet for stake measurements
%
% OUTPUT
% Cell array where each cell is a structure for each stake
%
%
% Written by Richard Buzard, January 17, 2018
%% Setup
[filename, filepath ]= uigetfile('*.xlsx','Choose Baseline Datasheet');
addpath(cd)
cd(filepath)

%% Read in data
[num, txt, raw] = xlsread(filename);
stakeNames = {'A','B','C','D','E'};
numProf  = unique(num(:,1));

%% Make Transect file
for jj = 1:length(numProf)          % for each profile
    idx = (num(:,1)==numProf(jj));    % index stake row
    numStake = sum(idx);   % find how many stakes are on profile
    out{jj}.prof   = numProf(jj); % get profile number
    out{jj}.stakeX = num(idx,3);                        % get x coordinates
    out{jj}.stakeY = num(idx,4);                        % get y coordinates
    dx = out{jj}.stakeX(1)-out{jj}.stakeX(end);         % calculate offsets
    dy = out{jj}.stakeY(1)-out{jj}.stakeY(end);
    out{jj}.bearing = atan2d(dx,dy);                    % calculate bearing
    if sign(out{jj}.bearing) == -1                      % if bearing is negative
        out{jj}.bearing = out{jj}.bearing + 360;        % correct with one rotation
    end
    %bearing = asind(abs(dx)/sqrt(dx^2 + dy^2)); % calculate bearing
%     waterDir = [sign(dx),sign(dy)];                     % determine direction to water
%     if sum(waterDir) == 2                               % if NW, do nothing
%     elseif and(sum(waterDir) == 0, waterDir(1) == 1)    % if SW, add 90
%         out{jj}.bearing = out{jj}.bearing + 90;
%     elseif sum(waterDir) == -2                          % if SE, add 180
%         out{jj}.bearing = out{jj}.bearing + 1800;
%     elseif and(sum(waterDir) == 0, waterDir(2) == 1)    % if NE, add 270
%             out{jj}.bearing = out{jj}.bearing + 270;
%     end
    
    

%%% Get the measurements, corrected to distance from landward stake
    % Need to keep this datum in mind for image comparison
    kk = 7:size(num,2);                                 % get number of measurements
    %ss = 1:(numStake-1);                                % for each seaward stake
    %dx = out{jj}.stakeX(ss) - out{jj}.stakeX(end);      % get offsets
    %dy = out{jj}.stakeY(ss) - out{jj}.stakeY(end);
    %out{jj}.stakeDist(ss) = sqrt(dx.^2 + dy.^2);        % calc stake distances
    dx = out{jj}.stakeX - out{jj}.stakeX(end);          % get offsets
    dy = out{jj}.stakeY - out{jj}.stakeY(end);
    out{jj}.stakeDist = sqrt(dx.^2 + dy.^2);            % calc stake distances
    efDist = num(idx,kk)*0.3048;                        % get erosion measurements in meters
    efDist = efDist + out{jj}.stakeDist;                % calculate distance from landward stake to bluff edge
    out{jj}.efDist = nansum(efDist,1);                  % clean up measurements
    out{jj}.efDist(out{jj}.efDist == 0) = NaN;          % turn 0s to NaN
    %out{jj}.efDist = efDist(ss,:) + out{jj}.stakeDist';
    %out{jj}.efDist(isnan(out{jj}.efDist)) = [];
    out{jj}.erosion = out{jj}.efDist - out{jj}.efDist(1);   % calc erosion
    out{jj}.efDate = txt(2,kk);                             % get date
    out{jj}.efDatenum = datenum(txt(2,kk),'mm/dd/yyyy');
end

%% Plot it
clf
for ii = 1:length(out)
    %plot(out{ii}.efDatenum,out{ii}.efDist,'--d')
    w = ii/length(out);
    c = [w/2 w/2 w];

    % remove nan data value
    x = datetime(datestr(out{ii}.efDatenum(find(~isnan(out{ii}.erosion)))));
    y = out{ii}.erosion(find(~isnan(out{ii}.erosion)));
    if ~isnan(y)
    plot(x,3.28084*y,'--o','linewidth',2,'MarkerSize',4, 'color',c)
    hold on
    lglabel{ii} = num2str(ii); 
    ax = gca;
    ax.XLim = [x(1) x(end)];
    end
end

legend(lglabel,'Location','southwest');

ax = gca;
%ax.XLim = [x(1) x(end)];
ax.YAxisLocation = 'right';
ylabel('Erosion (meters)')
ylabel('Erosion (feet)')
title(filename(1:3))
ax.XMinorTick = 'on';               % add month ticks
minorLim = dateshift(ax.XLim, 'start', 'month');
ax.XAxis.MinorTickValues = minorLim(1):calmonths(1):minorLim(2);
datetick('x','mm/yyyy','keeplimits')
ax.Units = 'inches';
ax.Position = [0.25 0.25 6 3.5];
end