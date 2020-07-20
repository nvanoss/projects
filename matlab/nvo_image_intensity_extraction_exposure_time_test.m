clear all
disp 'hi jude!'
close all

%% 
processing_folder = cd;

% Set and save data folder
cd '../../Reefat' %/Data Files'
data_folder = cd;

% Vector containing file numbering
%times = [5, 10, 20, 50, 100, 200, 300, 500];
times = [100:20:240];

% Enter number of image files
[~,N] = size(times);

%%

% Create data vectors
notGreen = 2:1:N+1; % no 1 or 0 values to avoid mistaken booleans
greenData = zeros(1,N);
greyData = zeros(1,N);
green_shift = zeros(1,N);
grey_shift = zeros(1,N);

for k = 0:1:N-1
    i = k + 1;
    % Get data from each file
	filename = sprintf('./Esposure_Time_Correction/Exposure Time Correction (1e4-1e8)_03-03-20/1e8/%i ms_Binning 3 by 3_1e8.png', times(i));
    %'./EOT_Fluorescence_10um_Filter/Fluorescence/1e8/5 ms_Binning 3 by 3_%02i.png', k);   
    A = importdata(filename);
    %B = imread(filename); % same output as line above
    grey = rgb2gray(A);

    % Separate color values from 3D RGB array
    red = A(:,:,1); %red
    green = A(:,:,2); %green
    blue = A(:,:,3); %blue

    % Returns index of non-zero elments stepping rows then columns
    % [1 3] array index values
    % [2 4]
    yellow = find(red);
    cyan = find(blue);

    % Check if image is pure green
    if ~(isempty(yellow) && isempty(cyan))
        notGreen(i) = true;
    else
        notGreen(i) = false;
    end

    % Sum green and grey data
    intensity_by_column = sum(green);
    greenData(i) = sum(intensity_by_column);
    intensity_by_column = sum(grey);
    greyData(i) = sum(intensity_by_column);
    
    % Calculate shift per reading
    if i > 1
        green_shift(i) = greenData(i) - greenData(k);
        grey_shift(i) = greyData(i) - greyData(k);
    else
        green_shift(i) = 0;
        grey_shift(i) = 0;
    end
    
end

% Construct data array for output
headers = ["Exposure Time, ms"; "0 for all green image"; ...
    "Green Intensity"; "Grey Intensity"; ...
    "Green Intensity per ms Exposure Time"; ...
    "Grey Intensity per ms Exposure Time"];

% Calculate intensity change
% Plot intensity shift relative to initial intensity

green_intensity_per_ms_exposure = greenData./times;
grey_intensity_per_ms_exposure = greyData./times;

figure(1)
subplot(1,2,1)
plot(times, greenData, 'g*-')
hold on
plot(times, greyData, 'k*-')
legend('Green', 'Grey', 'location', 'northwest')
title('Intensity per Reading')
xlabel('Exposure Time, ms')
ylabel('Intensity, Saturation at 255 for Green & 150 for Grey')

subplot(1,2,2)
hold on
plot(times, green_intensity_per_ms_exposure, 'g*-')
plot(times, grey_intensity_per_ms_exposure, 'k*-')
legend('Green', 'Grey', 'location', 'northwest');
title('Intensity per ms Exposure Time')
xlabel('Exposure Time, ms')
ylabel('Intensity per ms Exposure')

data = [times; notGreen; greenData; greyData; ...
    green_intensity_per_ms_exposure; grey_intensity_per_ms_exposure];
data = [headers, data];
columnData = transpose(data);

% Return to processing folder
cd (processing_folder);

%D = im2double(green); % converts to double and normalizes

% convert the intensity to double to make sure the entire table 
% is a double, which allows x and y to be greater than 255
%data(i) = double(A(:,:,2));
