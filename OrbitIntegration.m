clc
clear all
close all

%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

% Load kernel
cspice_furnsh ( METAKR );

planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN','301','5','VENUS','4','6'};

observer = 'EARTH';% or 339

% global G;
% G = 6.67e-20; % km % or -17


%% Ephemeris from SPICE

% Define initial epoch for a satellite
initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2030 NOV 21 11:22:23.659';% NOV! %'2030 DEC 28 00:03:25.693'; %'2030 DEC 28 00:03:25.693';%'2030 NOV 21 11:22:23.659';
%'2030 DEC 28 00:03:25.693'; % 7 months

initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );

step = 3600; %86400; %86400 3600 - every hour

% Create et time vector
et_vector = initial_et:step:end_et;

energy = zeros(3, length(et_vector));  % 1 row Kinetic, 2 row Potential, 3 row - Total Mechanical

energy_ab4 = zeros(3, length(et_vector));

% Satellite initial position w.r.t the Earth center
initial_state = [-561844.307770134;-1023781.19884100;-152232.354717768;0.545714129191316;-0.288204299060291;-0.102116477725135]; 

% Create a structure for a satellite
sat = create_sat_structure(initial_state);

% Get initial states for calculating initial energy
[earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init] = create_structure( planets_name_for_struct, initial_et, observer);


%% ODE Integration
% Case without influence from other planets
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
options87 = odeset('RelTol',1e-3,'AbsTol',1e-6);
global influence;
influence = zeros(3,2);
pressure = 1; %0 if no solar pressure needed

tic
orbit = ode45(@(t,y) force_model(t,y),et_vector,initial_state,options);    
toc

tic 
[orbit_ab8, tour] = adambashforth8(@force_model,et_vector,initial_state, length(et_vector), step);
toc

tic 
[orbit_rkv89, tourrkv] = RKV89(@force_model,et_vector,initial_state, length(et_vector), step);
toc

tic 

[tour1, orbit_ode87] = ode87(@(t,y) force_model(t,y),et_vector,initial_state, options87);   
toc
orbit_ode87 = orbit_ode87';


%% Mechanical Energy

% First calculate the initial energies
% b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
% [init_total, init_kinetic, init_potential] = calculate_energy(b);
% Initial_energy = init_total;
% Initial_kinetic = init_kinetic;
% Initial_potential = init_potential;



% 
load('irassihalotime.mat', 'Date')
load('irassihalogmat.mat', 'Gmat')

%% Plotting

figure(1)
subplot(1,2,1)
view(3)
grid on
hold on
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'r')% 
%plot3(orbit_ab8(1,:),orbit_ab8(2,:),orbit_ab8(3,:),'g')
subplot(1,2,2)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'b')% 

figure(3)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'b');% Reference
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'r');% RK45
plot3(orbit_ab8(1,:),orbit_ab8(2,:),orbit_ab8(3,:),'g'); % ABM8
plot3(orbit_rkv89(1,:),orbit_rkv89(2,:),orbit_rkv89(3,:),'m'); % RKV89
plot3(orbit_ode87(1,:),orbit_ode87(2,:),orbit_ode87(3,:),'y'); % RK87

%% Plots info
figure(1)
title('Integrated ephemeris of a satellite w.r.t the Earth, 3D');
subplot(1,2,1)
legend('Integrated Orbit RK', 'Integrated Orbit AB4');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
subplot(1,2,2)
legend('GMAT orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(3)
title('Reference vs Integration');
legend('Reference','RK45','ABM8', 'RKV89', 'RK87');
xlabel('x');
ylabel('y');
grid on



%cspice_kclear;