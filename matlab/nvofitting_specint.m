clear all
disp 'hi jude!'
close all

%% Setup Varaibles

% Save processing folder
processing_folder = cd;

% Set and save data folder
cd '../../Reefat/Data Files/HIV Paper/Flow_Through_03-08-20(from beads converted on 11-25-19)'
data_folder = cd;

% Enter number of spectra files
N = 11;

% Enter number of concentrations
n = 8;

% Enter number of data sets
S = 3;
    
% Line ABOVE where numerical data starts
firstDataLine = 17;

% Enter wavelength range
fixedOffset = 70;
first = 1535; % 600nm index
last = 2320; % 800nm index

%%

% Generate gaussian fittings and determine which display 
% a forward trend of peak shift

% Get an array of colors to plot
colors = lines(N);

% Initialize spectral integration array
specInt = zeros(S,n,N);
normSpecInt = zeros(S,n,N);
fixedSpecInt = zeros(S,n,N);
fixedNormSpecInt = zeros(S,n,N);

% Applying 1, 2, 3, 4, and 5 degree gauss fits
%figure(1)
for g = 1:1:S
    
    switch g
        case 1
            cd (data_folder)
            cd './Beads_from_11-25-19/EOT'
        case 2
            cd (data_folder)
            cd './Beads_from_03-11-20/EOT'
        otherwise
            cd (data_folder)
            cd './Beads_from_03-13-20/EOT'
    end
    
    for h = 1:1:n
        
            figure()
            title(sprintf('Experiment %i, Concentration %iF', g, h))
        
        for k = 0:1:N-1

            i = k + 1;

            % Get the data from each data file
            filename = sprintf('./%iF/%iF_%02i.txt', h, h, i);

            A = importdata(filename,'\t',firstDataLine);
            X = A.data(:,1);
            Y = A.data(:,2);
            Yn = Y/max(Y);

            % Consider only data in wavelength range
            peakX = X(first:last);
            peakY = Y(first:last);

            [~,v] = min(peakY);
            [~,u] = max(peakY);
            [s,~] = size(X);

            % Initialize Fixed Spectral Integration Range and Tracking Arrays
            if i == 1
                offset = first + u;
                if h == 1 && g == 1
                    range = u - v - fixedOffset;
                    intShift = zeros(S,n,range+1);
                    intShiftNorm = zeros(S,n,range+1);
                    int = zeros(range+1,1);
                    int0 = zeros(range+1,1);
                    int0norm = zeros(range+1,1);
                    two93 = zeros(range+1,1);
                end
            end

            % Set fixed range for spectral integration
            fixedX = X(offset-range:offset);
            fixedY = Y(offset-range:offset);
            fixedYn = Yn(offset-range:offset);

                hold on
                plot(fixedX,fixedY)

            % Spectral Integration Calculation on raw data
            for j = 1:1:range
                int(j) = fixedX(j+1) - fixedX(j);
                % Calculate Small Spectral Shift
                if i == N
                    intShift(g,h,:) = int(j);
                    intShiftNorm(g,h,:) = int(j);
                    two93(:,1) = intShift(g,h,:);
                    intShift(g,h,:) = two93.*((fixedY - int0)./int0).^2;
                    two93(:,1) = intShiftNorm(g,h,:);
                    intShiftNorm(g,h,:) = two93.*((fixedYn - int0norm)./int0norm).^2;
                end
                int(j) = int(j) * fixedY(j);
                if i == 1
                    int0(j) = fixedY(j);
                    int0norm(j) = fixedYn(j);
                end
            end
            fixedSpecInt(g,h,i) = sum(int);
                     
            % Spectral Integration Calculation on normalized data
            for j = 1:1:range
                int(j) = X(offset-range+j+1) - X(offset-range+j);
                int(j) = int(j) * Yn(offset-range+j);
            end
            fixedNormSpecInt(g,h,i) = sum(int);
            
        end
    end
end

% Return to processing folder
cd (processing_folder);

%% Create data array and plots

data = [fixedSpecInt; fixedNormSpecInt; specInt; normSpecInt];

intShift(:,:,range+1) = [];
intShiftNorm(:,:,range+1) = [];
intShiftSum = sqrt(sum(intShift,3));
intShiftNormSum = sqrt(sum(intShiftNorm,3));
intShiftAvg = sum(intShiftSum);
intShiftNormAvg = sum(intShiftNormSum);

% Average Spectral Integration per Concentration across Data Sets
fixedSpecIntAvg = sum(fixedSpecInt)./S;
fixedSpecIntAvgShift = fixedSpecIntAvg(1,:,1) - fixedSpecIntAvg(1,:,N);
fixedNormSpecIntAvg = sum(fixedNormSpecInt)./S;
fixedNormSpecIntAvgShift = fixedNormSpecIntAvg(1,:,1) - fixedNormSpecIntAvg(1,:,N);

% Plot Spectral Integration shift relative to initial reading
figure()
subplot(2,2,1)
plot(fixedSpecIntAvgShift, 'r*-')
title('Fixed Window Averaged Spectral Integration Shift vs Concentration Level')
xlabel('Concentration, F')
ylabel('Spectral Integration')

subplot(2,2,2)
plot(fixedNormSpecIntAvgShift, 'r*-')
title('Fixed Window Average Noramilzed Spectral Integration Shift vs Concentration Level')
xlabel('Concentration, F')
ylabel('Spectral Integration')

subplot(2,2,3)
plot(intShiftAvg, 'r*-')
title('Small Average Spectral Integration Shift vs Concentration Level')
xlabel('Concentration, F')
ylabel('Spectral Integration')

subplot(2,2,4)
plot(intShiftNormAvg, 'r*-')
title('Small Average Normalized Spectral Integration Shift vs Concentration Level')
xlabel('Concentration, F')
ylabel('Spectral Integration')

rangeStart = X(offset - range);
rangeEnd = X(offset);
chartTitle = sprintf('Fixed Integration Window from %3.0f nm to %3.0f nm', rangeStart, rangeEnd);
sgtitle(chartTitle)

% Make data folder active folder for convenient saving of data
cd (data_folder);
% Edit location in line below to where data is saved
cd ('./HIV Paper/Flow_Through_03-08-20(from beads converted on 11-25-19)/Beads_from_03-11-20/EOT');
% Apply breakpoint on line below to store data in the location above easily
cd (processing_folder);