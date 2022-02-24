function [time,ef]=StakeMeasure(~) 
%
% This function is built to load photos of two stakes with known positions and scales
% for the identification of their distance to eroding features or shorelines.
%
%
% INPUTS
% Folder of images with acquisition dates
%
% OUTPUTS
% Graph of shoreline / eroding feature movement over time
% Table of feature's distance from stake at certain datetime
% 
% Original written by Jacquelyn R. Overbeck March 14, 2016
% Revised to make measurements on x-axis by Richard M. Buzard June 29, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP
clear all; close all; clc;
addpath(pwd)
folder_name = uigetdir(pwd,'Choose Image Folder');      % prompt user for photos   
addpath (folder_name)
photo_dir = dir(folder_name);
photo_dir(1) = [];                     % remove . and ..
photo_dir(1) = [];

% test for images vs folder(s) of images
if photo_dir(1).isdir == 0          % if contents are not directories:
   photos = photo_dir;              % index the photos as usual.
else                                % if contents are directories:
    photos = [];                    % stack them up into a big struct
    dirCount = length(photo_dir);
    for kk = 1:dirCount             
        addpath(photo_dir(kk).name)
        temp =   dir(photo_dir(kk).name);
        photos = [photos; temp(3:end)];   
    end
end

% Get date from EXIF data rather than date modified
warning('off','MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero')
for pp = 1:length(photos)
    temp = imfinfo(photos(pp).name);    % get metadata for each photo
    photos(pp).date = datestr(datetime(temp.DateTime,'InputFormat','yyyy:MM:dd HH:mm:ss'));
    photos(pp).datenum = datenum(photos(pp).date);
end

S = [photos.datenum].';             % sort photos by date
[~,S] = sort(S);
photos = photos(S);
%% Calibrate Site
%%% Calibrate the distance between stakes and get pixel-to-meter conversion
%%% factor
%%% 1. Asks user to click the base of stakes in imagery.
%%% 2. Calculates pix2m conversion factor, position of seaward stake base
%%% 3. QC check allows user to repeat for better standard deviation

buttonCal = questdlg('Does calibration data exist?',...
         '','Y','N','N');

if buttonCal == 'Y';
    [calib, calPath] = uigetfile({'*.mat'}, 'Select Transect Calibration File');
    load(strcat(calPath, calib))

else %%% Run calibration!   

    %%% Load the photo
    I     = imread(photos(1).name);     
    sz    = size(I);                     
    xline = [sz(1)*2;1];            
    yline = [sz(2)*2;1];     
    f = figure;
    dI = imreducehaze(I);
    imshow(dI)  
    hFig = gcf;
    hAx  = gca;
    set(hFig,'units','normalized','outerposition',[0 0 1 1]);
    set(hAx,'Unit','normalized','Position',[0 0 1 1]);
    hold on             
    totalImages = length(photos);

    %%% Decide if photo is good for calibration
    button = questdlg('Are stakes visible?',...
               '','Y','N','N');
    if button == 'N';
        [I, filepath]= uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
              '*.*','All Files' },'Find image with visible stakes','multiselect','off');
        cd(filepath)
    end
    clf
    imshow(I,'InitialMagnification','fit')
    hFig = gcf;
    hAx  = gca;
    set(hFig,'units','normalized','outerposition',[0 0 1 1]);
    set(hAx,'Unit','normalized','Position',[0 0 1 1]);
    hold on
    %%% Get an average distance between stakes.
    %%% Get average position of the seaward stake.
    button = 'N';
    while button == 'N';
        allx = [];
        ally = [];       
        for ii = 1:10
            set(f,'name',sprintf('Click the base of landward stake, then seaward stake. (%d of 10)', ii))
            title(sprintf('Click the base of landward stake, then seaward stake. (%d of 10)', ii))
            [tempx, tempy] = ginput(2);
            allx = [allx, tempx];
            ally = [ally, tempy];
        end

        base = [mean(allx(2,:)),mean(ally(2,:))];   
        seaward = sign(allx(2,1)-allx(1,1));
        %%% Now that we have 10 xy values for each stake, we measure the distance
        %%% between each selection, then average those distances

        alldx = [];
        for jj = 1:10
            dx_temp = abs(allx(1,jj) - allx(2,jj));
            alldx = [alldx, dx_temp]; 
        end

        aveDist = mean(abs(alldx));
        stdDist = std(alldx);

        %%% Calculate pixel-to-meter conversion factor
        prompt       = {'Enter known distance between stakes (in feet)'}; 
        realDist     = str2double(inputdlg(prompt))*0.3048; % convert to meters         
        pix2m        = aveDist/realDist;
        measuredStdm = stdDist/pix2m;

        %%% Let user check for wonky data before they go too deep!
        button = questdlg(['Distance between stakes is ',num2str(realDist),...
            ' +/- ', num2str(measuredStdm), ' meters. Is this okay?'],...
            'Quality Control','Y','N','N');
    end
    
    %%% Get true edge information, identify transect
    prompt = {'Enter known distance from seaward stake to bluff (in feet)'}; 
    insitu = str2double(inputdlg(prompt))*0.3048; % convert to meters;
    edgex  = base(1) + pix2m*insitu*seaward;    % position of bluff edge on x axis
    I      = imread(photos(1).name);
    imshow(I) 
    hold on
    plot([edgex, edgex], [yline(1),yline(2)], '--r'); % Plots bluff intersect
    set(f,'name',sprintf('Click where the line crosses the shoreline or eroding feature'))
    [~ , vertAdjust] = ginput(1);  
    
    % Get coordinates of base and bluff edge to determine transect line
    % stake coordinate is base, bluff coordinate is [edgex, vertAdjust]
    edge = [edgex, vertAdjust];
    
    %%% Make the Calibration File
    [calibFile, calPath]=uiputfile('*.mat','Save Calibration Data as .mat file',folder_name);
    %cd(calPath)
    save(strcat(calPath,calibFile),'base','pix2m','xline','edge','vertAdjust','measuredStdm')
    close(f)
end
seaward = sign(edge(1)-base(1));    % determine seaward direction

%% OPTIONAL - only use images at noon
cd(folder_name) 
button = questdlg('Only use images at Noon?',...
         'Noon Sampling','Y','N','N');
temp = []; 
if button == 'Y'   
    for ii=1:length(photos)         % identify all images at noon
        if strcmp(photos(ii,1).date(13:14),'12')
           temp = [temp; photos(ii,1)];
        end
    end
    photos = temp;                  % reduce photo index to noon images
else buttonSub = questdlg('Measure image every X Hours','Sub-Sampling',...
         '1','2','4','1');    % sub sample images by hours
     subSample = str2double(buttonSub);
     dateRange = 5:subSample:21;    % date range between 5am and 11pm
     for ii=1:length(photos)         % for each photo
         for kk = 1:length(dateRange)   % for each interval 
            if strcmp(photos(ii,1).date(13:14),num2str(dateRange(kk),'%02.f'))
                temp = [temp; photos(ii,1)];
            end
         end
     end 
     photos = temp;
end
%% Measure eroding feature in each image
%%% User clicks position of feature in image. Code calculates distance
%%% and records associated time the image was taken.
cd(folder_name)
nfiles = length(photos);
ef   = nan(1,nfiles);
time = nan(1,nfiles);

f = figure;
for ii= 1:length(photos)
    cla;
    I = imread(photos(ii,1).name);
    imshow(I)
    hFig = gcf;
    hAx  = gca;
    %set(hFig,'units','normalized','outerposition',[0 0 1 1]);
    set(hAx,'Unit','normalized','Position',[0 0 1 1]);
    hold on
    %plot([xline(1) xline(2)],[vertAdjust vertAdjust],'-.r');
    plot([base(1) edge(1)],[base(2) edge(2)],'-r')
    plot([xline(1) xline(2)],[(vertAdjust+300) (vertAdjust+300)],'--k');
    % plot([xline(1) xline(2)],[base(2) base(2)],'-r'); % Plots the transect line
    plot(base(1),base(2),'ob');
    title([sprintf('Click where the line crosses the shoreline or eroding feature. (%i of %d) ',ii,nfiles),photos(ii,1).date])
    set(f,'name',([sprintf('Click where the line crosses the shoreline or eroding feature. (%i of %d) ',ii,nfiles),photos(ii,1).date]))
    [efpixvals,YTest] = ginput(1);% Eroding Feature x value
    if YTest < (vertAdjust+300)         % if user clicks below black line, omits that photo
        ef(ii)       = (base(1) - efpixvals(1))/pix2m ;
        time(ii)     = datenum(photos(ii,1).date); 
    elseif ii > 1
        ef(ii)       = ef(ii-1);
        time(ii)     = time(ii-1);
    end
end

ef = -(ef)*seaward; % make all values seaward of seaward stake positive

close(f)
%% plot it

figure;
plot (datenum(time),ef,':ob');
ylabel('Distance from Seaward Stake (meters)');
grid on
datetick('x',1,'keepticks')
% xmin = time(1);
% xmax = time(end);
% set(gca, 'Xlim', [xmin xmax])
% 
% choice=questdlg('Which date format?',...
%     'Date Format','dd-mmm-yyyy HH:MM:SS','dd-mmm-yyyy','mmmyy','dd-mmm-yyyy HH:MM:SS');
% 
% if strcmp(choice,'dd-mmm-yyyy HH:MM:SS')==1;
%       dateFormat=0; %dd-mmm-yyyy HH:MM:SS
% elseif  strcmp(choice,'dd-mmm-yyyy');
%     dateFormat=1; %dd-mmm-yyyy
% elseif strcmp(choice,'mmmyy');
%     dateFormat=12; %mmmyy
% else dateFormat=0;
% end
% set(gca, 'Xlim', [xmin xmax]);
% datetick('x',dateFormat);


[file_name2, file_path2]=uiputfile(strcat(calib(1:7),photo_dir(1).folder(end-16:end),'.txt'),'Save Results');
cd(file_path2)

dateTime = datestr(time);
A = [ef', time'];

save(file_name2,'A','-ascii');


% fid = fopen(file_name2,'w');
% fprintf(fid, '%s%2.2f,%f,%s\n', 'Distance', A(:,1));
% fprintf(fid, '%s%2.2f,%f,%s\n', 'Datenum', A(:,2));
% fprintf(fid, '%s%2.2f,%f,%s\n', 'Datetime', A(:,3))
 fclose 'all';

end
