%% Close and Clear all
close all           % close figures
clear all           % clear variables

%% IDT parameters and input frequency
fw = 0.000083;      % finger width
lambdaf = 4*fw;     % Rayleigh SAW wavelength calculated from finger width
f = 11830000;       % frequency
w = 2*pi*f;         % angular frequency
aper = 15e-3;       % aperture

%% LiNbO3 Properties 
% http://www.roditi.com/SingleCrystal/LiNbO3/liNBO3-Properties.html
csaw = 3931;        % Leaky SAW Phase Velocity Real Component 3980;
vr = 3994;          % Rayleigh SAW Phase Velocity 3882;
density = 4647;     % LiNbO3 density
ksq = 0.055;        % LiNbO3 electromechanical coupling

%% Substrate Properties 
pf = 998;           % fluid density
vw = 1500;          % fluid longitudinal wave velocity - Shiokawa %1498;
mu = 0.001;         % shear viscosity
mup = 0.00247;      % bulk viscosity

%% SAW velocity, wavelength, amplitude, and Rayleigh Angle Equation
vl = csaw + 1i*67.7;% Leaky SAW Complex Velocity 3884 + 1i*70.7; % Shiokawa
lambda = vr/f;          % wavelength calculated from velocity
kr = w/real(vl);        % real component of Leaky SAW propagation constant
A = 0.000000001;        % SAW amplitude
theta = asin(vw/vr)*360/(2*pi);    % Rayleigh angle

%% Attenuation and wave propagation constants
lsaw = density*vr*lambda/(pf*vw);       % attenuation length
kr2 = 2*pi/lambda;                      % wavelength based calculation of Leaky SAW propagation constant
kr3 = w/vl;                             % complex Leaky SAW propagation constant
a = sqrt(1-(vr/vw)^2); % **** correct!  % attenuation constant alpha
al = 1i*a;                              % multiply by i to get real number

%% Area of the cross section of the wave perpendicular to wave motion
y = -0.01:0.0001:0.01;              % consider 5mm to either side
yamp = A * exp(-(lsaw^-1)*abs(y));  % peak SAW amplitude attenuated along y-axis
yarea = sum(yamp*(y(2)-y(1)));      % integrated cross sectional area

%% IDT impedance data
Zc = 182;               % resistive impedance (equal to resistance)
Zr = 20; %25;           % capacitive impedance
C = 1/(Zc*w);           % effective capacitance

%% Power calculations
% Simulation voltage variable
V = 25;
% estimated acoustic power for shear wave
P1 = 0.5 * density * yarea * w^2 * A^2 * vr;
% Hashimoto acoustic power equation
P2 = (w * C/aper * ksq^2 * V^2);
% Electrical consumption of device 1 modeled after Bourquin and Cooper IDT
P3 = (V^2) / sqrt(Zc^2 + Zr^2);
% Q-factor
Q = 1 / (w*Zr*C);
% Bourquin and Cooper claimed 0.72 dB attenuation in electrical losses
point72 = 10*log10(P2/(P3)); % -7dB // *(Q/(Q+1))

%% Imaginary component of propagation constant
% ki = -1370; % Shiokawa
ki = -lsaw^-1;      % imaginary component of Leaky SAW propagation constant
kl = kr + 1i * ki;  % complete complex Leaky SAW propagation constant
wl = kl * vl;       % angular frequency of Leaky SAW

%% Generate spacial and time ranges
P = 100;                    % number of points per wavelength
N = lsaw/lambda;            % number of wavelengths in attenuation length
points = N*P;               % total number of points
points = round(points);     % make integer value
lattice = 1:1:points;       % generate 1-D lattice
x = linspace(0,lsaw,points);% generate x variable
t = x./real(vl);            % generate time variable
z = zeros(points);          % initialize loop arrays
Fx = zeros(points);
Fz = zeros(points);
F1s = zeros(points);
Fsaw = zeros(points);
Fintt = 1:1:points;         % initialize loop vectors for integrated force
Fints = 1:1:points;

%% Generate Leaky SAWs with different initial phase conditions
% each phase angle will fill a row, columns represent x distance travelled
for i = 1:1:points % initial phase angles
    % Create different attenuation profiles for different phase angles
    z(i,:) = A*sin(i*2*pi/points+2*pi.*lattice/P).*exp(2*ki*x);
    % X-component of streaming force
    Fx(i,:) = -pf*(1+al^2)*(A*w)^2*ki*exp(2*(ki*x+al*ki*z(i,:)));
    % Z-component of streaming force
    Fz(i,:) = -pf*(1+al^2)*(A*w)^2*al*ki*exp(2*(ki*x+al*ki*z(i,:)));
    % Sum component forces
    F1s(i,:) = sqrt(Fx(i,:).^2+Fz(i,:).^2);
    % Unified force equation
    Fsaw(i,:) = -pf*(1+al^2)^(3/2)*(A*w)^2*ki*exp(2*(ki*x+al*ki*z(i,:)));
    % Time integrated streaming force
    Fintt(i) = sum(Fsaw(i,:)*(t(2)-t(1)));  % unified force per unit time
    % Time integrated streaming force
    Fints(i) = sum(F1s(i,:)*(t(2)-t(1)));   % component force per unit time
end

%% Unified and Component force equivalence check
check = Fsaw - F1s;         % equal to 10e-12

%% Check that total force is directed at Rayleigh angle
rayleigh = atan(Fx./Fz)*360/(2*pi); % ****** correct!

%% Amplitude envelope and volume calculations
% volume of amplitude envelope
volume = sum(yarea*exp(2*ki*x)*(x(2)-x(1)));
Z = zeros(length(y),points);        % 3D x-y grid to generate z values
vol = zeros(points,points);         % 3D x-y grid to generate volume values
for i = 1:1:points
    % apply y-axis and x-axis attenuation to initial amplitude for envelope
    Z(:,i) = A*exp(-(lsaw^-1)*abs(y))*exp(2*ki*x(i));
    % generate volume values with a column for each phase angle
    vol(:,i) =  yarea * z(i,:) / A; % rows contain equivalent x locations
end

%% Average volume over time
volav = mean(sum(abs(vol)));    % // order of mean and sum operations null

%% Redefine new Amplitude variable with values varying by 1pm from 1pm to 1nm
Amp = 0.0000000001:0.0000000001:0.000000001;% vary amplitude for force plot

%% Force extraction
Finttav = mean(Fintt);              % 11.5 mN @ 1nm A
Fintsav = mean(Fints);              % 11.5 mN @ 1nm A
Fintsav = Fintsav * Amp.^2 / (A^2); % Amplitude replacement
Finttav = Finttav * Amp.^2 / (A^2);
Fav = sum(mean(Fsaw));
Fvol = Fav * volav * Amp.^3 / (A^3); % 17.4 mN @ 1nm A

% check that sqrt force varies with A => force varies linearly with power
pwr1 = sqrt(Finttav)./Amp;                  % constant
pwr2 = sqrt(Fvol)./Amp;                     % constant

% estimated acoustic power for shear wave
pwr3 = 0.5 * density * yarea * sum(abs(sin(2*pi*lattice(1:P)/P))) * w^2 * A^2 * vr;

%% Chassis terminal velocity calculation
Cd = 1.45;%:0.10:1.15;                          % drag coefficient
area = 0.002368666;                             % footprint of chassis rev3
width = 0.03;                                   % ideal width of chassis
chavol = area * lsaw;                           % ideal submersion volume
mass = pf * chavol;                             % mass based on fluid disp
acceleration = Finttav(7:8) / mass;             % acceleration from force
time = 0:0.01:1;                                % time lattice
u = [time*acceleration(1);time*acceleration(2)];% velocity lattice
Fd = 0.5 * pf * u.^2 * Cd * width * lsaw;       % drag force
t7 = Fd(1,:) - Finttav(7);                      
t8 = Fd(2,:) - Finttav(8);                      % subtract Fd from Fsaw
[~,t7t] = min(abs(t7));                         
[~,t8t] = min(abs(t8));                         % find index of Fd=Fsaw
u7 = u(1,t7t);
u8 = u(2,t8t);                                  % find velocity of index

%% Power calculation sketchbook
distance = round(lambda/2,6);
depth = 0.000001:0.000001:distance;
attenuation = (depth(2)-depth(1))*exp(depth);
volcheck = 1.7 / (density * w^2 * A^2 * f);
%depth = sum(attenuation);
%p0 = density * yarea * depth * w^2 * A^2 * f;
% P1 = 0.5 * density * yarea * w^2 * A^2 * vr;
% W = kg*m^2/s^3
% px = kg/m^3 * m^3 * 

%% Plots
% yamp plot
figure()
plot(y,yamp)
xlabel('Distance along y-axis, m')
ylabel('Leaky SAW Amplitude, m')
title('SAW Attenuation across y-axis')

% Time Integrated Force vs Peak SAW Amplitude plot
figure()
plot(Amp*1e9,Finttav*1e3)
xlabel('Peak SAW Amplitude, nm')
ylabel('SAW Streaming Force, mN')
title('SAW Streaming Force vs SAW Amplitude')

% Time Averaged Force vs Peak SAW Amplitude plot
figure()
plot(Amp*1e9,Fvol*1e3);
xlabel('Peak SAW Amplitude, nm')
ylabel('Time Averaged Streaming Force, mN')
title('SAW Streaming Force vs SAW Amplitude')

% Leaky SAW Amplitude Envelope Mesh
figure()
mesh(x,y,Z)
xlabel('x, distance into fluid along axis of propagation')
ylabel('y, horizontal distance from center of SAW beam')
zlabel('z, amplitude of SAW')
title('Leaky SAW Ampitude Envelope')

% Leaky SAW Propagation Animated Line Plot
figure()
X = zeros(length(y),points);    % generate x-y grid
for k = 1:1:length(y)
    X(k,:) = x;
end

% Plot Stationary Wave
g = animatedline('Color','b','LineWidth',0.1,'LineStyle',':','MaximumNumPoints',points*length(y));
axis([0,lsaw*1e3,y(1)*1e3,y(length(y))*1e3,-A*1e9,A*1e9])
xlabel('x, distance into fluid along axis of propagation, mm','FontSize',16)
ylabel('y, horizontal distance from center of SAW beam, mm','FontSize',16)
zlabel('z, amplitude of SAW, nm','FontSize',16)
%title('Leaky SAW Propagation Simulation')
for i = 1:1:points
    addpoints(g,X(:,i)*1e3,y*1e3,Z(:,i)*cos(i*2*pi/P)*1e9);
    drawnow
end

% Plot Propagating Wave
clearpoints(g)
h = animatedline('Color','b','LineWidth',0.1,'LineStyle',':','MaximumNumPoints',points*length(y));
axis([0,lsaw*1e3,y(1)*1e3,y(length(y))*1e3,-A*1e9,A*1e9])
for k = 1:1:P
    clearpoints(h)
    for i = 1:1:points
        addpoints(h,X(:,i)*1e3,y*1e3,Z(:,i)*cos(k*8*pi/P-i*2*pi/P)*1e9);
    end
    drawnow
end
