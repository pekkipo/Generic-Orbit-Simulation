
clc
clear all
close all

%% Define local variables
METAKR = 'planetsorbitskernels.txt';
moon = '301';
earth = 'EARTH';
sun = 'SUN';
mars = '4';
rev_amount = 1;
cspice_furnsh ( METAKR );
observer = earth;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
G=6.67e-20;

initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2032 JUN 19 00:03:25.693';
% Convert to UTC
initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );
step = 86400;
%rev_period = 86400*27; % 27 - Moon around Earth, 365 - Earth around Sun, 687 - Mars around Sun
%end_et_days = initial_et + rev_periods(chosen_body_index);%225*86400; % for one Venus revolution

% Create et time vector
et_vector = initial_et:step:end_et;%end_et;
% Instead of my for loop
initial_state = cspice_spkezr ( moon, initial_et, 'J2000', 'NONE', observer ); % need it for integration

% Ephemeris
ephemeris = cspice_spkezr ( moon, et_vector, 'J2000', 'NONE', observer );

% Integrated orbit
    % RK
    orbit = ode45(@(t,y) verification_force_model(t,y),et_vector,initial_state,options);
    % AB 
    [orbit_ab4, tour] = adambashforth8(@verification_force_model,et_vector,initial_state, length(et_vector), step);

% Difference
%difference1 = orbit.y - ephemeris;
difference2 = orbit_ab4(:,1:length(orbit_ab4)-1) - ephemeris;


% Energy
planets_simplified = {'EARTH', 'SUN', 'MOON';'EARTH', 'SUN', '301'};
energy = zeros(3, length(et_vector)); % RK
energy1 = zeros(3, length(et_vector)); % AB
% Get initial states for calculating initial energy
[earth_init, sun_init, moon_init] = simplified_create_structure( planets_simplified, initial_et, observer);

dkin_dpot_rk = zeros(1, length(et_vector));
dkin_dpot_ab = zeros(1, length(et_vector));

% First calculate the initial energies
b = [earth_init, sun_init, moon_init];
%[init_total, init_kinetic, init_potential] = calculate_energy(b);
[init_total, init_kinetic, init_potential] = simple_energy(moon_init, earth_init);
Initial_energy = init_total;
Initial_kinetic = init_kinetic;
Initial_potential = init_potential;



% Calculate for each step
for epoch = 1:length(et_vector)
    
    % Moon structure
     field1 = 'name';  value1 = 'Moon';
    field2 = 'x';  value2 = [];
    field3 = 'y';  value3 = [];
    field4 = 'z';  value4 = [];
    field5 = 'vx';  value5 = [];
    field6 = 'vy';  value6 = [];
    field7 = 'vz';  value7 = [];
    field8 = 'mass'; value8 = [];
    field9 = 'GM'; value9 = [];
    field10 = 'coords'; value10 = [];
    
    Moon = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
    Moon.GM = cspice_bodvrd( '301', 'GM', 1 );
    Moon.mass = (Moon.GM)/G;
    Moon_coords = cspice_spkezr ( '301', et_vector(epoch), 'J2000', 'NONE', observer );
    Moon.x = Moon_coords(1);
    Moon.y = Moon_coords(2);
    Moon.z = Moon_coords(3);
    Moon.vx = Moon_coords(4);
    Moon.vy = Moon_coords(5);
    Moon.vz = Moon_coords(6);
    Moon.coords = Moon_coords(1:3);
    %
        
    [Earth, Sun] = simplified_create_structure( planets_simplified, et_vector(epoch), observer);
    bodies = [Earth, Sun, Moon];
   % [total, kinetic, potential] = calculate_energy(bodies);
   [total, kinetic, potential] = simple_energy(Moon, Earth);
    kin1 = kinetic - Initial_kinetic;
    pot1 = potential - Initial_potential;
    tot1 = total - Initial_energy;
    energy(1,epoch) = kin1;
    energy(2,epoch) = pot1;
    energy(3,epoch) = tot1;
    
    dkin_dpot_rk(1,epoch) = kin1 - pot1;
end

for epoch1 = 1:length(et_vector)
    
    % Moon structure
     field1 = 'name';  value1 = 'Moon';
    field2 = 'x';  value2 = [];
    field3 = 'y';  value3 = [];
    field4 = 'z';  value4 = [];
    field5 = 'vx';  value5 = [];
    field6 = 'vy';  value6 = [];
    field7 = 'vz';  value7 = [];
    field8 = 'mass'; value8 = [];
    field9 = 'GM'; value9 = [];
    field10 = 'coords'; value10 = [];
    
    Moon1 = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
    Moon1.GM = cspice_bodvrd( '301', 'GM', 1 );
    Moon1.mass = (Moon1.GM)/G;
    Moon_coords1 = cspice_spkezr ( '301', et_vector(epoch1), 'J2000', 'NONE', observer );
    Moon1.x = Moon_coords1(1);
    Moon1.y = Moon_coords1(2);
    Moon1.z = Moon_coords1(3);
    Moon1.vx = Moon_coords1(4);
    Moon1.vy = Moon_coords1(5);
    Moon1.vz = Moon_coords1(6);
    Moon1.coords = Moon_coords1(1:3);
    %
    
    [Earth1, Sun1] = simplified_create_structure( planets_simplified, et_vector(epoch1), observer);
    bodies1 = [Earth1, Sun1, Moon1];
    %[total1, kinetic1, potential1] = calculate_energy(bodies1);
    [total1, kinetic1, potential1] = simple_energy(Moon1, Earth1);
    kin2 = kinetic1 - Initial_kinetic;
    pot2 = potential1 - Initial_potential;
    tot2 = total1 - Initial_energy;
    energy1(1,epoch1) = kin2;
    energy1(2,epoch1) = pot2;
    energy1(3,epoch1) = tot2;
    
    dkin_dpot_ab(1,epoch1) = kin2 - pot2;
end


figure(1)
subplot(1,2,1)
view(3)
grid on
hold on
plot3(ephemeris(1,:),ephemeris(2,:),ephemeris(3,:),'r') % Ephemeris from SPICE
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'b')% Ephemeris from integration
subplot(1,2,2)
view(3)
grid on
hold on
plot3(ephemeris(1,:),ephemeris(2,:),ephemeris(3,:),'r') % Ephemeris from SPICE
plot3(orbit_ab4(1,:),orbit_ab4(2,:),orbit_ab4(3,:),'b')

figure(2)
view(3)
grid on
hold on
plot3(ephemeris(1,:),ephemeris(2,:),ephemeris(3,:),'r') % Ephemeris from SPICE
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'g')%
plot3(orbit_ab4(1,:),orbit_ab4(2,:),orbit_ab4(3,:),'b')

figure(3)
view(2)
hold on
plot(et_vector(1,:),difference2(1,:), 'r');
plot(et_vector(1,:),difference2(2,:), 'g');
plot(et_vector(1,:),difference2(3,:), 'b');

figure(4)
subplot(1,2,1)
view(2)
hold on
plot(et_vector(1,:), energy(1,:), 'r');
plot(et_vector(1,:), energy(2,:), 'g');
plot(et_vector(1,:), energy(3,:), 'b');
subplot(1,2,2)
view(2)
hold on
plot(et_vector(1,:), energy1(1,:), 'r');
plot(et_vector(1,:), energy1(2,:), 'g');
plot(et_vector(1,:), energy1(3,:), 'b');

figure(5)
view(2)
hold on
plot(et_vector(1,:),dkin_dpot_rk(1,:), 'r');
plot(et_vector(1,:),dkin_dpot_ab(1,:), 'b');





%% Plots info
figure(1)
subplot(1,2,1)
title('SPICE vs RK');
legend('SPICE','Integrated');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
subplot(1,2,2)
title('SPICE vs AB');
legend('SPICE','Integrated');
xlabel('x');
ylabel('y');
zlabel('z');

figure(2)
title('Merged Graph');
legend('SPICE','RK','AB');
xlabel('x');
ylabel('y');
zlabel('z');
grid on


figure(4)
subplot(1,2,1)
title('Energy RK');
legend('kinetic','potential','total');
xlabel('E');
ylabel('t');
grid on
subplot(1,2,2)
title('Energy AB');
legend('kinetic','potential','total');
xlabel('E');
ylabel('t');
grid on



