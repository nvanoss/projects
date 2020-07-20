clear all
disp 'hi jude!'
close all

%% 
processing_folder = cd;

% Set and save data folder
cd '../../Reefat/Data Files'
data_folder = cd;

% Enter number of image files
N = 11;

%%

% Create data vectors
notGreen = 2:1:N+1; % no 1 or 0 values to avoid mistaken booleans
greenData = zeros(1,N);
greyData = zeros(1,N);
green_shift = zeros(1,N);
grey_shift = zeros(1,N);

i = 1;
for k = 0:1:N-1
    
    % Get data from each file
	filename = sprintf('./HIV Paper/Flow_Through_03-08-20(from beads converted on 11-25-19)/Fluorescence/5 ms_Binning 3 by 3_%02i.png', k);   
    A = importdata(filename);
    %B = imread(filename); % same output as line above
    grey = rgb2gray(A);

    % Separate color values from 3D RGB array
    red = A(:,:,1); %red
    green = A(:,:,2); %green
    blue = A(:,:,3); %blue

    % Returns index of non-zero elments stepping rows then columns
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
    
    i = i + 1;
end

% Construct data array for output
headers = ["Index"; "0 for all green image"; "Green Intensity"; ...
    "Grey Intensity"; "Cumulative Green Intensity Change"; ...
    "Green Shift per Reading"; "Cumulative Grey Intensity Change"; ...
    "Grey Shift per Reading"];

% Calculate intensity change
% Plot intensity shift relative to initial intensity

green_intensity_change = greenData - greenData(1);
grey_intensity_change = greyData - greyData(1);

figure(1)
subplot(1,3,1)
plot(greenData, 'g*-')
hold on
plot(greyData, 'k*-')
legend('Green', 'Grey', 'location', 'southeast')
title('Intensity per Reading')
xlabel('Reading #')
ylabel('Intensity, Saturation per Pixel at 255 for Green & 150 for Grey')
subplot(1,3,2)
hold on
plot(green_intensity_change, 'g*-')
plot(grey_intensity_change, 'k*-')
legend('Green', 'Grey', 'location', 'southeast');
title('Cumulative Intensity Change')
xlabel('Reading #')
ylabel('Intensity Shift from Initial')
subplot(1,3,3)
hold on
plot(green_shift, 'g*-')
plot(grey_shift, 'k*-')
legend('Green', 'Grey', 'location', 'southeast')
title('Intensity Increase per Reading')
xlabel('Reading #')
ylabel('Intensity Shift from Preceeding Reading')

data = [0:1:N-1; notGreen; greenData; greyData; green_intensity_change; ...
    grey_intensity_change; green_shift; grey_shift];
data = [headers, data];
columnData = transpose(data);

% Return to processing folder
cd (processing_folder);

%D = im2double(green); % converts to double and normalizes

% convert the intensity to double to make sure the entire table 
% is a double, which allows x and y to be greater than 255
%data(i) = double(A(:,:,2));
