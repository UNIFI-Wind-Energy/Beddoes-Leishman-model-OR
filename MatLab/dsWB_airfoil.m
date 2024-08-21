clc
clear
close all

addpath BL

% ------------------------------------------------------------------------ INPUT

% pitch law

A0 = 14;                                                                    % harmonic pitch law average value [deg]
A1 = 10;                                                                    % harmonic pitch law amplitude [deg]

k = 0.077;                                                                  % reduced frequency [-]

% inflow conditions

M = 0.1;                                                                    % freestream Mach number
T0 = 298.15;                                                                % freestream temperature [K]

% airfoil properties

chord = 0.457;                                                              % airfoil chord [m]
x_AC = 0.25;                                                                % airfoil aerodynamic center [x/c] - depends on Mach

% numerical set-up

N_rev = 10;                                                                 % number of pitching cycles
N = 180;                                                                    % number of timesteps per cycle

formulation = 'incompressible';                                             % formulation of the attached flow module: incompressible | compressible
fMode = 'raw';                                                              % handling of the f function: 'fit' use Leishman exponential fitting | 'raw' use data from static polars
vortexModule = 'on';                                                        % activate LEV module: 'on' | 'off'
timeConstantsMod = 'on';                                                    % activate modification of time constants due to LEV: 'on' | 'off'
secondaryVortex ='on';                                                      % activate secondary vortex shedding: 'on' | 'off'

% input data

file_validation = 'reference data/14+10_k0077_M01.txt';                     % reference data
file_constants = 'polarData/S809_constants.txt';                            % BL model constants
file_polar = 'polarData/S809_Re1000k.txt';                           % airfoil polar data

% ------------------------------------------------------------------------






%% initialization

% input validation data

validationData = importdata(file_validation);

AOA_val = validationData(:,1);
CL_val = validationData(:,2);
CD_val = validationData(:,3);
CM_val = validationData(:,4);

CN_val =  CL_val .* cosd(AOA_val) + CD_val .* sind(AOA_val);
CC_val =  CL_val .* sind(AOA_val) - CD_val .* cosd(AOA_val);

% input calibration data

calibrationData = importdata(file_constants).data;

% input polar data

polarData = load(file_polar);

AOA = polarData(:,1);
CL = polarData(:,2);
CD = polarData(:,3);
CM = polarData(:,4);

CN =  CL .* cosd(AOA) + CD .* sind(AOA);
CC =  CL .* sind(AOA) - CD .* cosd(AOA);

% derived parameters

a = sqrt(1.4*287*T0);                                                       % freestream speed of sound [m/s] 
V = M*a;                                                                    % freestream velocity [m/s]

omega = (2*k*V)/chord;                                                      % pitching pulsation [1/rad]
T = 2*pi/omega;                                                             % pitching period [s]


%% Pitch simulation

% generate pitch law

t = linspace(0, N_rev*T, N_rev*N);
aoa_f= deg2rad(A0) + deg2rad(A1) * sin(omega*t);
aoa_rate = deg2rad(A1) * omega * cos(omega*t);
theta_rate = aoa_rate;
h_rate = aoa_rate * 0.0;

% initialize vectors

cl=zeros(1,length(t),1);
cd=zeros(1,length(t),1);
cn=zeros(1,length(t),1);
ct=zeros(1,length(t),1);
cm=zeros(1,length(t),1);
f_lag=zeros(1,length(t),1);
tv=zeros(1,length(t),1);
comp=zeros(length(t),9);
bl=zeros(length(t),6);
state = zeros(1, 28);

dt = t(2)-t(1);

% run pitching simulation

for i=1:length(t)

    [cn(i), ct(i), cl(i), cd(i), cm(i), f_lag(i), tv(i), comp(i,:), bl(i,:), state] = BL(aoa_f(i), aoa_rate(i), theta_rate(i), h_rate(i), V, M, dt, chord, x_AC, calibrationData, polarData, formulation, fMode, timeConstantsMod, vortexModule, secondaryVortex, state);

end




%% plot data

range_plot = (N_rev - 1) * T <= t & t <= N_rev*T;

% figure 1 - simulation history

fig=figure(1);

fig.Units='normalized';
fig.Position=[-0.0005    0.0380    0.3495    0.8833];

xlimits = [(N_rev-1) N_rev];

subplot(8,1,1)
plot(t/T, rad2deg(aoa_f),'-k')
hold on
title('AoA')
xlim(xlimits)
ylim([A0-A1-5 A0+A1+5])
grid on

subplot(8,1,2)
plot(t/T, aoa_rate,'-k')
title('aoarate')
hold on
xlim(xlimits)
ylim([min(aoa_rate) max(aoa_rate)])
grid on

subplot(8,1,3)
plot(t/T, f_lag,'-k')
hold on
title('f')
xlim(xlimits)
ylim([0 1])
grid on

subplot(8,1,4)
plot(t/T, tv/calibrationData(33),'-k')
hold on
title('tv/Tvl')
xlim(xlimits)
ylim([0 3])
grid on

subplot(8,1,5)
plot(t/T, bl(:,1),'-k')
hold on
title('Tf')
xlim(xlimits)
ylim([0 20])
grid on

subplot(8,1,6)
plot(t/T, bl(:,2),'-k')
hold on
title('Tv')
xlim(xlimits)
ylim([0 10])
grid on

subplot(8,1,7)
plot(t/T, cn,'-k')
title('CN')
xlim(xlimits)
ylim([0 2]);
grid on

subplot(8,1,8)
plot(t/T, ct,'-k')
hold on
title('CC')
xlabel('t/T [-]')
xlim(xlimits)
ylim([-0.1 0.5])
grid on

% figure 2 - CN, CC

aoa_plot = rad2deg(aoa_f(range_plot));

cl = cl(range_plot);
cd = cd(range_plot);
cm = cm(range_plot);
cn = cn(range_plot);
ct = ct(range_plot);

fig2 = figure(2);

fig2.Units='normalized';
fig2.Position=[0.3578    0.6269    0.4036    0.2926];

subplot(121)
plot(AOA, CN, '-o', aoa_plot, cn, '-k', AOA_val, CN_val, '-r')
hold on
title('CN')
xlabel('AoA [deg]')

subplot(122)
plot(AOA, CC,'-o', aoa_plot, ct, '-k', AOA_val, CC_val, '-r')
title('CC')
xlabel('AoA [deg]')

% figure 3 - CL, CD, CM

fig3 = figure(3);

fig3.Units = 'normalized';
fig3.Position = [0.7729    0.0380    0.2000    0.8833];

xlimits = [-5 30];

subplot(311)
plot(AOA, CL, '-o', aoa_plot, cl, '-k', AOA_val, CL_val, '-r')
xlim(xlimits)
hold on
title('CL')

subplot(312)
plot(AOA, CD, '-o', aoa_plot, cd, '-k', AOA_val, CD_val, '-r')
title('CD')
xlim(xlimits)
hold on

subplot(313)
plot(AOA, CM, '-o', aoa_plot, cm, '-k', AOA_val, CM_val, '-r')
hold on
title('CM')
xlabel('AoA [deg]')
xlim(xlimits)

legend('static data', 'BL model','EXP','Location','SouthWest','fontSize',14)


%% output to file

P = [aoa_plot' cn' ct' cd' cm'];

dlmwrite(sprintf("Output/results_%d+%d_k%g_M%g.txt", A0,A1,k,M), P);

exportgraphics(fig3, sprintf("Output/results_%d+%d_k%g_M%g.pdf",A0,A1,k,M),'Resolution',300)
