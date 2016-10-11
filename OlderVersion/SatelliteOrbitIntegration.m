clc
clear all
close all


%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';
% Now including the SUN. EARTH is zero but I ll leave for possible future
% needs
planets = {'SUN','MERCURY','VENUS','EARTH','4','5','6','7','8','301'}; % FOR SPEZER function
planets_names = {'SUN','MERCURY','VENUS','EARTH','MARS','JUPITER','SATURN','URANUS','NEPTUNE','MOON'}; % for getting GM's

observer = 'EARTH';% or 339

%% Ephemeris from SPICE
% Load kernel
cspice_furnsh ( METAKR );

planets_initial_states = {
    planets_names{1},planets_names{2},planets_names{3},planets_names{4},planets_names{5},planets_names{6},planets_names{7},planets_names{8},planets_names{9},planets_names{10};
    zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1);
    };

% Define initial epoch for a satellite
initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2032 DEC 28 00:03:25.693'; % 7 months

initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );

% Create et time vector
et_vector = initial_et:86400/4:end_et;

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

%% Integation Part                  

% Get GM's of planets
planets_gms = zeros(1,10);
for i = 1:length(planets_gms)
    planets_gms(1,i) = cspice_bodvrd( planets_names{i}, 'GM', 1 );
end


%% Mechanical Energy
global energy;
energy = zeros(3, length(et_vector)); % 1 row Kinetic, 2 row Potential, 3 row - Total Mechanical
global G;
G = 6.67e-17;


% Initial energy (now only considering Earth, Sun and the satellite)
% satellite mass is negligible, so = 0 thus we get simplified equation for
% mechanical energy

% Getting mass of planets in kg requires multiplying GM from spice by 10^9
% and dividing by G = 6.67e-11 m3kg-1s-2

Sun_mass = (planets_gms(1,1) * 10^3)/G;
Earth_mass = (planets_gms(1,4) * 10^3)/G;

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
    %orbit(:,n) = initial_state;
    point = ode45(@(t,y) iter_new_sat_force_model(t,y,planets_gms,planetspoints, n, pressure),[et_vector(1) et_vector(2)],initial_state,options);
    energy(1,1) = 0;
    energy(2,1) = 0;
    energy(3,1) = 0;
    orbit(:,n) = point.y(:,length(point.x));
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

% Transpose for convenience
orbit = orbit';

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
plot(et_vector(1,:), energy(1,:), 'r');
plot(et_vector(1,:), energy(2,:), 'g');
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
