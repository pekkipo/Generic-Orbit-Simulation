clc
clear all
close all

%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

% Load kernel
cspice_furnsh ( METAKR );

planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN','301','5','VENUS','4','6'};

observer = 'EARTH';% or 339

global G;
G = 6.67e-20; % km % or -17


%% Ephemeris from SPICE

% Define initial epoch for a satellite
initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2030 DEC 28 00:03:25.693'; % 7 months

initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );

step = 86400;

% Create et time vector
et_vector = initial_et:step:end_et;

energy = zeros(3, length(et_vector));  % 1 row Kinetic, 2 row Potential, 3 row - Total Mechanical

% Satellite initial position w.r.t the Earth center
initial_state = [-561844.307770134;-1023781.19884100;-152232.354717768;0.545714129191316;-0.288204299060291;-0.102116477725135]; 

% Create a structure for a satellite
sat = create_sat_structure(initial_state);
% field1 = 'name'; value1 = 'Satellite';
% field2 = 'x'; value2 = initial_state(1);
% field3 = 'y'; value3 = initial_state(2);
% field4 = 'z'; value4 = initial_state(3);
% field5 = 'vx'; value5 = initial_state(4);
% field6 = 'vy'; value6 = initial_state(5);
% field7 = 'vz'; value7 = initial_state(6);
% field8 = 'mass'; value8 = 6000;
% field9 = 'GM'; value9 = 0;
% field10 = 'coords'; value10 = [initial_state(1);initial_state(2);initial_state(3)];
% sat = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6, field7,value7, field8,value8, field9,value9, field10,value10);

% Get initial states for calculating initial energy
[earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init] = create_structure( planets_name_for_struct, initial_et, observer);


%% ODE Integration
% Case without influence from other planets
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
global influences;
influences = zeros(6,7);
pressure = 1; %0 if no solar pressure needed

tic
orbit = ode45(@(t,y) force_model(t,y,observer),et_vector,initial_state,options);    
toc

%% ODE integration looped

looped_orbit = zeros(6, length(et_vector));

tic
for p = 1:length(et_vector)
    
    % Create structure for one point
    [Point_Earth, Point_Sun, Point_Moon, Point_Jupiter, Point_Venus, Point_Mars, Point_Saturn] = create_structure( planets_name_for_struct, et_vector(p), observer);
    str_planets = [Point_Earth; Point_Sun; Point_Moon; Point_Jupiter; Point_Venus; Point_Mars; Point_Saturn];
    
    if p == 1   
    looped_orbit(:,p) = initial_state;
    %point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(1) et_vector(2)],initial_state,options);
    energy(1,1) = 0;%init_kinetic; % 0;
    energy(2,1) = 0;%init_potential;%0;
    energy(3,1) = 0; %init_total;%0;
    elseif and(p > 1, p < length(et_vector))
    new_initial_state = looped_orbit(:,p-1); 
    point = ode45(@(t,y) looped_force_model(t,y,str_planets, pressure),[et_vector(p) et_vector(p+1)],new_initial_state,options);     % usual - n and n+1
    looped_orbit(:,p) = point.y(:,length(point.x));
    elseif p == length(et_vector)
    new_initial_state = looped_orbit(:,p-1);
    point = ode45(@(t,y) looped_force_model(t,y,str_planets, pressure),[et_vector(p-1) et_vector(p)],new_initial_state,options);     
    looped_orbit(:,p) = point.y(:,length(point.x));
    end
end
toc
% Transpose for convenience
looped_orbit = looped_orbit';

%% Mechanical Energy

% First calcilate the initial energies
b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
[init_total, init_kinetic, init_potential] = calculate_energy(b);
Initial_energy = init_total;
Initial_kinetic = init_kinetic;
Initial_potential = init_potential;

% Calculate for each step
for epoch = 1:length(et_vector)
    
    % Create a structure for the satellite
    sat_at_this_time = create_sat_structure(orbit.y);
    % Information about planets at a given epoch
    [earth, sun, moon, jupiter, venus, mars, saturn] = create_structure( planets_name_for_struct, et_vector(epoch), observer);
    bodies = [sat_at_this_time, earth, sun, moon, jupiter, venus, mars, saturn];
    [total, kinetic, potential] = calculate_energy(bodies);
    energy(1,epoch) = kinetic;
    energy(2,epoch) = potential;
    energy(3,epoch) = total;
end

%% Plotting

figure(1)
view(3)
grid on
hold on
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'r')% 
plot3(looped_orbit(:,1),looped_orbit(:,2),looped_orbit(:,3),'b') % loop

figure(2)
view(2)
grid on
hold on
plot(et_vector(1,:), energy(1,:), 'r');
plot(et_vector(1,:), energy(2,:), 'g');
plot(et_vector(1,:), energy(3,:), 'b');

%% Plots info
figure(1)
title('Integrated ephemeris of a satellite w.r.t the Earth, 3D');
legend('Integrated Orbit', 'Integrated with a loop');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(2)
title('Simplified energy');
%legend('Kietic', 'Potential', 'Total mechanical energy');
legend('Kinetic energy','Potential energy','Total energy');
xlabel('x');
ylabel('y');
grid on

%cspice_kclear;