clear all
disp 'hi jude!'
close all

%% Setup Varaibles
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

% Enter number of fittings + 1 for raw data
f = 6;

% Line ABOVE where numerical data starts
firstDataLine = 17;

% Enter (arbitrary) wavelength range
first = 1875;%1925; % 700nm index
last = 2270;%2320; % 800nm index

% get an array of colors to plot
colors = lines(N);

% create array to store centroid data
cent = [2,S,n,f,N];

%% Centroid data extraction and storage

% top variable for experiment date variation
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

    % variable for concentration variation
    for h = 1:1:n

        figure()
        sgtitle(sprintf('Experiment %i Concentration %iF',g,h))

        %applying 1, 2, 3, 4, and 5 degree gauss fits
        for k = 0:N-1:N-1

            i = k + 1;

            % Get the data from each data file
            filename = sprintf('./%iF/%iF_%02i.txt', h, h, i);

            A = importdata(filename,'\t',firstDataLine);
            X = A.data(:,1);
            Y = A.data(:,2);

            % Consider only data in wavelength range
            [~,u] = max(Y);

            % Set Fixed Spectral Integration Range
            if i == 1
                offset = u + 100;
                range = 200;
            end

            peakX = X(offset-range:offset+range);
            peakY = Y(offset-range:offset+range);
            Y_fitted = zeros(2*range+1,5);

            %access processing functions
            cd (processing_folder);

            %generate the model function
            g1 = gauss1(peakX,peakY); 
            g2 = gauss2(peakX,peakY);
            g3 = gauss3(peakX,peakY);
            g4 = gauss4(peakX,peakY);
            g5 = gauss5(peakX,peakY);

            %evaluate the model function
            Y_fitted(:,1) = g1(peakX);
            Y_fitted(:,2) = g2(peakX);
            Y_fitted(:,3) = g3(peakX);
            Y_fitted(:,4) = g4(peakX);
            Y_fitted(:,5) = g5(peakX);

            %normalize
        %     Yn = Y/max(Y);
        %     Y_fitted_3 = Y_fitted_3/max(Y_fitted_3);
        %     Y_fitted_4 = Y_fitted_4/max(Y_fitted_4);
        %     Y_fitted_5 = Y_fitted_5/max(Y_fitted_5);

            %make polygons
            poly0 = polyshape(peakX,peakY);
            %warning('off','last')
            poly1 = polyshape(peakX,Y_fitted(:,1));
            poly2 = polyshape(peakX,Y_fitted(:,2));
            poly3 = polyshape(peakX,Y_fitted(:,3));
            poly4 = polyshape(peakX,Y_fitted(:,4));
            poly5 = polyshape(peakX,Y_fitted(:,5));

            hold on
            subplot(2,3,1)
            hold on
            plot(poly0)
            title('No Fit')
            subplot(2,3,2)
            hold on
            plot(poly1)
            title('1st Degree')
            subplot(2,3,3)
            hold on
            plot(poly2)
            title('2nd Degree')
            subplot(2,3,4)
            hold on
            plot(poly3)
            title('3rd Degree')
            subplot(2,3,5)
            hold on
            plot(poly4)
            title('4th Degree')
            subplot(2,3,6)
            hold on
            plot(poly5)
            title('5th Degree')

            %save the centroid of each array
            [cent(1,g,h,1,i),cent(2,g,h,1,i)] = centroid(poly0);
            [cent(1,g,h,2,i),cent(2,g,h,2,i)] = centroid(poly1);
            [cent(1,g,h,3,i),cent(2,g,h,3,i)] = centroid(poly2);
            [cent(1,g,h,4,i),cent(2,g,h,4,i)] = centroid(poly3);
            [cent(1,g,h,5,i),cent(2,g,h,5,i)] = centroid(poly4);
            [cent(1,g,h,6,i),cent(2,g,h,6,i)] = centroid(poly5);

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
            
        end
    end
end

%% Shift Calculations

centShift = zeros(S,n,f);
centShift(:,:,:) = sqrt((cent(1,:,:,:,N)-cent(1,:,:,:,1)).^2 + (cent(2,:,:,:,N)-cent(2,:,:,:,1)).^2);
e1 = zeros(n,f);
e1(:,:) = centShift(1,:,:);
e2 = zeros(n,f);
e2(:,:) = centShift(2,:,:);
e3 = zeros(n,f);
e3(:,:) = centShift(3,:,:);

%% Plots

figure()
hold on
plot(e1)
title('11/25 Beads Centroid Shift vs Concentration')
legend('No Fit','1st Order','2nd Order','3rd Order','4th Order','5th Order')
figure()
hold on
plot(e2)
title('3/11 Beads Centroid Shift vs Concentration')
legend('No Fit','1st Order','2nd Order','3rd Order','4th Order','5th Order')
figure()
hold on
plot(e3)
title('3/13 Beads Centroid Shift vs Concentration')
legend('No Fit','1st Order','2nd Order','3rd Order','4th Order','5th Order')

cd (processing_folder)
