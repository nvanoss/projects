clear all
close all

% Enter Specifications:  EDIT THESE 3 LINES
N = 12;                                     % Number of devices
wafer_id = '2';                             % Wafer number
id_offset = 11;                             % ID number of first device

% Store processing folder
processing_folder = cd;

% Set and save data folder
cd '../VNA S11 Data'
data_folder = cd;

% Create matrix to store relevant values
B = zeros(N,8);
header = ["S11 @ 9.2 MHz", "S11 @ 11.8 MHz", "S11 @ 13.2 MHz", ...
    "Z @ 11.8 MHz", "max transfer freq, MHz", ...
    "max S11 value", "max index", "S21 @ 9.2 MHz"];

% Import and plot data from each file

for k = 1:1:N
    
    % Generate file name and path
    device_id = sprintf('%02i', k + id_offset - 1);
	filename = strcat('./', wafer_id, '-', device_id, '.csv');
    plot_title = strcat({'Wafer #'}, wafer_id, {', Device #'}, device_id);
    
    % Import data
    A = importdata(filename);
    
    % Create data vectors
    freq = A(:,1)./1000000;
    s11 = A(:,2);
    s21 = A(:,4);
    
    % Plot
    figure
    plot(freq,s11,freq,s21)
    title(plot_title)
    xlabel('Frequency, MHz')
    ylabel('dB Relative to Input Signal')
    legend('S11: Reflection at Input', ...
        'S21: Transmission at Output','Location','south')
    
    % Harvest Data
    u = find(A(:,1)>=9100000,1);
    v = find(A(:,1)>=11700000,1);
    w = find(A(:,1)>=13100000,1);
    
    % Store 9.2, 11.8, 13.2, min freq, min, and min index in array
    B(k,1) = A(u,2);
    B(k,2) = A(v,2);
    s11 = (10^(B(k,2)/10));
    B(k,3) = A(w,2);
    B(k,4) = 50*(1+s11)/(1-s11);
    B(k,6) = min(A(u:w,2));
    [mins,dummy] = find(A(:,2)==B(k,6),N);
    
    B(k,8) = A(u,4);
    
    if size(mins) == 1
        B(k,7) = mins;
        B(k,5) = A(mins,1)/1000000;
    else
        B(k,7) = 0;
    end
    
end

str = string(B);
data = [header; str];

% Return to data processing folder
cd (processing_folder);
