clc
clear
close all

addpath BL

% ------------------------------------------------------------------------ INPUT

% inflow conditions

M = 0.1;                                                                    % freestream Mach number
rho = 998.2;
V0 = 2;

% turbine properties

R = 1;
chord = 0.25;                                                               % airfoil chord [m]
x_AC = 0.25;                                                                % airfoil aerodynamic center [x/c] - depends on Mach

TSR = 4;

% numerical set-up

N_rev = 10;                                                                 % number of pitching cycles
N = 180;                                                                    % number of timesteps per cycle

formulation = 'incompressible';                                             % formulation of the attached flow module: incompressible | compressible
fMode = 'fit';                                                              % handling of the f function: 'fit' use Leishman exponential fitting | 'raw' use data from static polars
vortexModule = 'on';                                                        % activate LEV module
timeConstantsMod = 'on';                                                    % activate modification of time constants due to LEV
secondaryVortex ='on';                                                      % activate secondary vortex shedding

% input data

file_validation = 'reference data/14+10_k0077_M01.txt';                     % reference data
file_constants = 'polarData/S809_constants.txt';                            % BL model constants
file_polar = 'polarData/S809_Re1000k_smooth.txt';                           % airfoil polar data

% ------------------------------------------------------------------------

%% initialization

% input polar and calibration data

validationData = importdata(file_validation);                         

calibrationData = importdata(file_constants).data;        

polarData = load(file_polar);                                              
 
AOA = polarData(:,1);
CL = polarData(:,2);
CD = polarData(:,3);
CM = polarData(:,4);

% Derived parameters

a = sqrt(1.4*287*298.15);

omega = V0/R*TSR;

%% Derived quantities 

theta = linspace(0, 2*pi*N_rev, N*N_rev);

aoa_f = -atan( sin(theta)./(TSR+cos(theta)) );

Vrel = sqrt( (V0*cos(theta)+omega*R).^2 + (V0*sin(theta)).^2);

aoa_rate = -omega*(TSR*cos(theta) + 1)./(TSR^2 + 2*TSR*cos(theta) + 1);

theta_rate = -omega;

h_rate = aoa_rate*0.0;

t = theta./omega;


%% Simulation set-up

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

dt=t(2)-t(1);


%% Dynamic stall simulation

for i=1:length(t)
            
    [cn(i), ct(i), cl(i), cd(i), cm(i), f_lag(i), tv(i), comp(i,:), bl(i,:), state] = BL(aoa_f(i), aoa_rate(i), theta_rate, h_rate(i), Vrel(i), M, dt,chord, x_AC, calibrationData, polarData, formulation, fMode, timeConstantsMod, vortexModule, secondaryVortex, state);

end


% Compute turbine forces

CN = cl .* cos(aoa_f) + cd .* sin(aoa_f);
CT = cl .* sin(aoa_f) - cd .* cos(aoa_f);

FN = 0.5*rho*CN.*Vrel.^2;
FT = 0.5*rho*CT.*Vrel.^2;


%% plot data

range_plot = (N_rev - 1) * 2*pi <= theta & theta <= N_rev*2*pi;

% figure 1

fig=figure(1);

fig.Units='normalized';
fig.Position=[-0.0005    0.0380    0.3495    0.8833];

x_plot = theta;
xlimits = [(N_rev-1) N_rev] * 2*pi;

subplot(8,1,1)
plot(x_plot, rad2deg(aoa_f),'-k')
hold on
title('AoA')
xlim(xlimits)
ylim([min(rad2deg(aoa_f))-5 max(rad2deg(aoa_f))+5])
grid on

subplot(8,1,2)
plot(x_plot, aoa_rate,'-k')
title('aoarate')
hold on
xlim(xlimits)
ylim([min(aoa_rate) max(aoa_rate)])
grid on

subplot(8,1,3)
plot(x_plot, f_lag,'-k')
hold on 
title('f')
xlim(xlimits)
ylim([0 1])
grid on

subplot(8,1,4)
plot(x_plot, tv/calibrationData(23),'-k')
hold on 
title('tv/Tvl')
xlim(xlimits)
ylim([0 3])
grid on

subplot(8,1,5)
plot(x_plot, bl(:,1),'-k')
hold on 
title('Tf')
xlim(xlimits)
ylim([0 20])
grid on

subplot(8,1,6)
plot(x_plot, bl(:,2),'-k')
hold on 
title('Tv')
xlim(xlimits)
ylim([0 10])
grid on

subplot(8,1,7)
plot(x_plot, cn,'-k')
title('CN')
xlim(xlimits)
ylim([0 2]);
grid on

subplot(8,1,8)
plot(x_plot, ct,'-k')
hold on 
title('CC')
xlabel('theta [°]')
xlim(xlimits)
ylim([-0.1 0.5])
grid on

% figure 2

aoa_plot = rad2deg(aoa_f(range_plot));

cl = cl(range_plot);
cd = cd(range_plot);
cm = cm(range_plot);
cn = cn(range_plot);
ct = ct(range_plot);

fig2 = figure(2);

fig2.Units = 'normalized';
fig2.Position = [0.3578    0.6269    0.4036    0.2926];

subplot(121)
plot(aoa_plot, cn,'-k', validationData(:,1), validationData(:,2), '-r')
hold on 
title('CN')
xlabel('AoA [deg]')

subplot(122)
plot(aoa_plot, ct,'-k')
title('CC')
xlabel('AoA [deg]')

% figure 3

fig3 = figure(3);

fig3.Units = 'normalized';
fig3.Position = [0.7729    0.0380    0.2000    0.8833];

xlimits = [-10 30];

subplot(311)
plot(aoa_plot, cl, '-k')
xlim(xlimits)
hold on 
title('CL')

subplot(312)
plot(aoa_plot, cd, '-k', validationData(:,3), validationData(:,4), '-r')
title('CD')
xlim(xlimits)
hold on

subplot(313)
plot(aoa_plot, cm, '-k', validationData(:,5), validationData(:,6), '-r')
hold on
title('CM')
xlabel('AoA [deg]')
xlim(xlimits)

legend('BL model','EXP','Location','SouthWest','fontSize',14)

% figure 4 - turbine loads

figure(4)

subplot(2,1,1)

plot(rad2deg(theta), FN,'-k')
hold on
plot([0 360],[0 0],'-.k')
xlim([0 360])
title('FN')

subplot(2,1,2)

plot(rad2deg(theta),FT,'-k')
hold on
plot([0 360],[0 0],'-.k')
xlim([0 360])
title('FT')

%% output to file

% P = [aoa_plot' cn' ct' cd' cm'];
% 
% dlmwrite(sprintf("Output/results_%d+%d_k%g_M%g.txt", A0,A1,k,M),P);
% 
% exportgraphics(fig3, sprintf("Output/results_%d+%d_k%g_M%g.pdf",A0,A1,k,M),'Resolution',300)



