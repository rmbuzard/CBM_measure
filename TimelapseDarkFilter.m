function [list_files] = TimelapseBlurFilter(photos)
%%% TimelapseBlurFilter

%%% Removes images that are blurry or dark
%%% Really, it looks at the histogram of the image and removes those with a
%%% distribution characteristic of a blurry or dark image

% requires System Identification Toolbox
if license('test','identification_toolbox') == 0
    error('Error: missing System Identification Toolbox. Go to Apps > Get More Apps to download and install toolbox.')
end

numPhotos = length(photos);

% Get statistics on the images
s = clock;
h = [];
for ii = 1:numPhotos

    img = imread(photos(ii).name);  % read in photo
    imgGray = rgb2gray(img);        % turn to grayscale
    H = mean(imgGray(:));           % read histogram mean
    if H > 70                       % if average value is greater than 70
        h = [h;ii];                 % store index number
    end
    % every 10 photos estimate time left
    if sum(ii == 1:10:numPhotos)
        is = etime(clock,s);
        esttime = is * numPhotos / ii - is;
    end
    clc
    disp(strcat(['Detecting dark photos: ',...
        num2str(ii), ' of ', num2str(numPhotos),...
        ' remaining time = ',...
        num2str(floor(esttime/60)),':',...
        num2str(round(rem(esttime,60)),'%02u')]))
end
clc

% store in output variable
list_files = photos(h);
end
