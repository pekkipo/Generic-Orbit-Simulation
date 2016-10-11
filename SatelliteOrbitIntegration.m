clc
clear all
close all


%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

% Load kernel
cspice_furnsh ( METAKR );


% Now including the SUN. EARTH is zero but I ll leave for possible future
% needs
planets = {'SUN','MERCURY','VENUS','EARTH','4','5','6','7','8','301'}; % FOR SPEZER function
planets_names = {'SUN','MERCURY','VENUS','EARTH','MARS','JUPITER','SATURN','URANUS','NEPTUNE','MOON'}; % for getting GM's

observer = 'EARTH';% or 339

global G;
G = 6.67e-11;

%% Ephemeris from SPICE


planets_initial_states = {
    planets_names{1},planets_names{2},planets_names{3},planets_names{4},planets_names{5},planets_names{6},planets_names{7},planets_names{8},planets_names{9},planets_names{10};
    zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1);
    };

% Define initial epoch for a satellite
initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2031 DEC 28 00:03:25.693'; % 7 months

initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );

% Create et time vector
et_vector = initial_et:86400:end_et;

% Satellite initial position w.r.t the Earth center
initial_state = [-561844.307770134;-1023781.19884100;-152232.354717768;0.545714129191316;-0.288204299060291;-0.102116477725135]; 

% Get initial positions of planets at the initial epoch
for i=1:length(planets)
    planets_initial_states{2,i} = cspice_spkezr ( planets{i}, initial_et, 'J2000', 'NONE', observer );
end

% CELL ARRAY
ephemerises = {
     planets{1},planets{2},planets{3},planets{4},planets{5},planets{6},planets{7},planets{8},planets{9},planets{10}; % 1 row
     zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6),zeros(length(et_vector),6);
     };
 

  % FILL IN THE EPHEMERISES
 for n = 1:length(ephemerises)
     ephemerises{2,n} = cspice_spkezr ( planets{n}, et_vector, 'J2000', 'NONE', observer );
      % transpose for convenience
     % ephemerises{2,n} = ephemerises{2,n}';
 end
 
 % Structurs stuff
    planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN','301','5','VENUS','4','6'};
    % First create planets structs
    
    for pl=1:length(planets_name_for_struct) 
    
    field1 = 'name';  value1 = planets_name_for_struct{1,pl};
    field2 = 'x';  value2 = zeros(1,length(et_vector));
    field3 = 'y';  value3 = zeros(1,length(et_vector));
    field4 = 'z';  value4 = zeros(1,length(et_vector));
    field5 = 'vx';  value5 = zeros(1,length(et_vector));
    field6 = 'vy';  value6 = zeros(1,length(et_vector));
    field7 = 'vz';  value7 = zeros(1,length(et_vector));
    field8 = 'mass'; value8 = [];
    field9 = 'GM'; value9 = [];
    field10 = 'coords'; value10 = [];
        if pl == 1
        Earth = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Earth.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Earth.mass = (Earth.GM * 10^9)/G;
        Earth_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Earth.x = Earth_coords(1,:);
        Earth.y = Earth_coords(2,:);
        Earth.z = Earth_coords(3,:);
        Earth.vx = Earth_coords(4,:);
        Earth.vy = Earth_coords(5,:);
        Earth.vz = Earth_coords(6,:);
        Earth.coords = Earth_coords;
        elseif pl == 2
        Sun = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Sun.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Sun.mass = (Sun.GM * 10^9)/G;  
        Sun_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Sun.x = Sun_coords(1,:);
        Sun.y = Sun_coords(2,:);
        Sun.z = Sun_coords(3,:);
        Sun.vx = Sun_coords(4,:);
        Sun.vy = Sun_coords(5,:);
        Sun.vz = Sun_coords(6,:);
        Sun.coords = Sun_coords;
        elseif pl == 3
        Moon = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Moon.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Moon.mass = (Moon.GM * 10^9)/G;
        Moon_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Moon.x = Moon_coords(1,:);
        Moon.y = Moon_coords(2,:);
        Moon.z = Moon_coords(3,:);
        Moon.vx = Moon_coords(4,:);
        Moon.vy = Moon_coords(5,:);
        Moon.vz = Moon_coords(6,:);
        Moon.coords = Moon_coords;
        elseif pl == 4
        Jupiter = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Jupiter.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Jupiter.mass = (Jupiter.GM * 10^9)/G;
        Jupiter_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Jupiter.x = Jupiter_coords(1,:);
        Jupiter.y = Jupiter_coords(2,:);
        Jupiter.z = Jupiter_coords(3,:);
        Jupiter.vx = Jupiter_coords(4,:);
        Jupiter.vy = Jupiter_coords(5,:);
        Jupiter.vz = Jupiter_coords(6,:);
        Jupiter.coords = Jupiter_coords;
        elseif pl == 5
        Venus = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Venus.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Venus.mass = (Venus.GM * 10^9)/G;
        Venus_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Venus.x = Venus_coords(1,:);
        Venus.y = Venus_coords(2,:);
        Venus.z = Venus_coords(3,:);
        Venus.vx = Venus_coords(4,:);
        Venus.vy = Venus_coords(5,:);
        Venus.vz = Venus_coords(6,:);
        Venus.coords = Venus_coords;
        elseif pl == 6
        Mars = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Mars.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Mars.mass = (Mars.GM * 10^9)/G;
        Mars_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Mars.x = Mars_coords(1,:);
        Mars.y = Mars_coords(2,:);
        Mars.z = Mars_coords(3,:);
        Mars.vx = Mars_coords(4,:);
        Mars.vy = Mars_coords(5,:);
        Mars.vz = Mars_coords(6,:);
        Mars.coords = Mars_coords;
        elseif pl == 7
        Saturn = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Saturn.GM = cspice_bodvrd( planets_name_for_struct{2,pl}, 'GM', 1 );
        Saturn.mass = (Saturn.GM * 10^9)/G;
        Saturn_coords = cspice_spkezr ( planets_name_for_struct{2,pl}, et_vector, 'J2000', 'NONE', observer );
        Saturn.x = Saturn_coords(1,:);
        Saturn.y = Saturn_coords(2,:);
        Saturn.z = Saturn_coords(3,:);
        Saturn.vx = Saturn_coords(4,:);
        Saturn.vy = Saturn_coords(5,:);
        Saturn.vz = Saturn_coords(6,:);
        Saturn.coords = Saturn_coords;
        end
    end


    
%% Integation Part                  

% Get GM's of planets
planets_gms = zeros(1,10);
for i = 1:length(planets_gms)
    planets_gms(1,i) = cspice_bodvrd( planets_names{i}, 'GM', 1 );
end


%% Mechanical Energy
global energy;
energy = zeros(3, length(et_vector)); % 1 row Kinetic, 2 row Potential, 3 row - Total Mechanical


% Initial energy (now only considering Earth, Sun and the satellite)
% satellite mass is negligible, so = 0 thus we get simplified equation for
% mechanical energy

% Getting mass of planets in kg requires multiplying GM from spice by 10^9
% and dividing by G = 6.67e-11 m3kg-1s-2

Sun_mass = (planets_gms(1,1) * 10^9)/G;
Earth_mass = (planets_gms(1,4) * 10^9)/G;

r_sun_earth = sqrt(ephemerises{2,1}(1,1)^2 + ephemerises{2,1}(2,1)^2 + ephemerises{2,1}(3,1)^2); % earth coords = 0
potential_energy = (G/2)*( (Earth_mass*(Sun_mass/r_sun_earth)) + (Sun_mass*(Earth_mass/r_sun_earth)) );

global Initial_energy;
kinetic_energy = Sun_mass*(sqrt(ephemerises{2,1}(4,1)^2 + ephemerises{2,1}(5,1)^2 + ephemerises{2,1}(6,1)^2))^2;
Initial_energy = kinetic_energy - potential_energy;




%% ODE Integration

% Case without influence from other planets
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% USUAL ODE
% [tour, orbit_usual] = ode45(@(t,y) new_sat_force_model(t,y,planets_gms, 0, observer),et_vector,initial_state,options);

% ODE WITH ITERATIONS AND DIFFERENT FUNCTION

global influences;
influences = zeros(6,7);


pressure = 1; %0 if no solar pressure needed

orbit = zeros(6, length(et_vector));
str_orbit = zeros(6, length(et_vector));

tic
for n = 1:length(et_vector)
    
    % Data for all planets
    planetspoints = zeros(6, length(planets));
   
    
    for m = 1:length(ephemerises)
        onepoint = ephemerises{2,m}(:,n); % the whole column of n-th column, ie data for one point
        planetspoints(:,m) = onepoint;
    end
    
    % could have avoided the stuff above just changing planetspoints to
    % ephemeris and function from iter to usual sat force model
    % but I'll leave it like it is now..'cause why not huh

    if n == 1   
    orbit(:,n) = initial_state;
    %point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(1) et_vector(2)],initial_state,options);
    energy(1,1) = 0;
    energy(2,1) = 0;
    energy(3,1) = 0;
    elseif and(n > 1, n < length(et_vector))
    new_initial_state = orbit(:,n-1); 
    point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(n) et_vector(n+1)],new_initial_state,options);     % usual - n and n+1
    orbit(:,n) = point.y(:,length(point.x));
    elseif n == length(et_vector)
    new_initial_state = orbit(:,n-1);
    point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(n-1) et_vector(n)],new_initial_state,options);     
    orbit(:,n) = point.y(:,length(point.x));
    end


% Take the last step (most precise) evaluation and set it as a point of an
% orbit

%orbit(:,n) = point.y(:,length(point.x));

end
toc

tic
for p = 1:length(et_vector)
    
    %str_planetspoints = zeros(1, length(planets_name_for_struct));
    
    
    % A bit lame but - Re-create structure but this time each value is only
    % for one point
    [Point_Earth, Point_Sun, Point_Moon, Point_Jupiter, Point_Venus, Point_Mars, Point_Saturn] = create_structure( planets_name_for_struct, et_vector(p), observer);
    str_planets = [Point_Earth; Point_Sun; Point_Moon; Point_Jupiter; Point_Venus; Point_Mars; Point_Saturn];
    % could have avoided the stuff above just changing planetspoints to
    % ephemeris and function from iter to usual sat force model
    % but I'll leave it like it is now..'cause why not huh
    
    if p == 1   
    str_orbit(:,p) = initial_state;
    %point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(1) et_vector(2)],initial_state,options);
    energy(1,1) = 0;
    energy(2,1) = 0;
    energy(3,1) = 0;
    elseif and(p > 1, p < length(et_vector))
    new_initial_state = str_orbit(:,p-1); 
    str_point = ode45(@(t,y) struct_iter_new_sat_force_model(t,y,str_planets, p, pressure),[et_vector(p) et_vector(p+1)],new_initial_state,options);     % usual - n and n+1
    str_orbit(:,p) = str_point.y(:,length(str_point.x));
    elseif p == length(et_vector)
    new_initial_state = str_orbit(:,p-1);
    str_point = ode45(@(t,y) struct_iter_new_sat_force_model(t,y,str_planets, p, pressure),[et_vector(p-1) et_vector(p)],new_initial_state,options);     
    str_orbit(:,p) = str_point.y(:,length(str_point.x));
    end
end
toc
% Transpose for convenience
orbit = orbit';
str_orbit = str_orbit';

% Plotting

figure(1)
view(3)
grid on
hold on
plot3(orbit(:,1),orbit(:,2),orbit(:,3),'r')% loop
% plot3(orbit_usual(:,1),orbit_usual(:,2),orbit_usual(:,3),'b') % Usual ODE

figure(2)
view(2)
grid on
hold on
%plot(et_vector(1,:), energy(1,:), 'r');
%plot(et_vector(1,:), energy(2,:), 'g');
plot(et_vector(1,:), energy(3,:), 'b');

%% Plots info
figure(1)
title('Integrated ephemeris of a satellite w.r.t the Earth, 3D');
legend('Integrated Orbit with a loop');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(2)
title('Simplified energy');
%legend('Kietic', 'Potential', 'Total mechanical energy');
legend('Total mechanical energy');
xlabel('x');
ylabel('y');
grid on

%cspice_kclear;