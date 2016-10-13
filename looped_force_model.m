function yp = looped_force_model( t,y0, bodies, pressure_on )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% pressure_on: 1 - including solar pressure, 0 - without solar pressure
 % bodies - vector of structures
 % Create a structure for the satellite
% field1 = 'name'; value1 = 'Satellite';
% field2 = 'x'; value2 = y0(1);
% field3 = 'y'; value3 = y0(2);
% field4 = 'z'; value4 = y0(3);
% field5 = 'vx'; value5 = y0(4);
% field6 = 'vy'; value6 = y0(5);
% field7 = 'vz'; value7 = y0(6);
% field8 = 'mass'; value8 = 6000;
% field9 = 'GM'; value9 = 0;
% field10 = 'coords'; value10 = [y0(1);y0(2);y0(3)];
% sat = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6, field7,value7, field8,value8, field9,value9, field10,value10);

sat = create_sat_structure(y0);
earth = bodies(1);
sun = bodies(2);
moon = bodies(3);
jupiter = bodies(4);
venus = bodies(5);
mars = bodies(6);
saturn = bodies(7);

%% Accelerations due to:

% GRAVITY

% y0 - satellite, rows: x y z vx vy vz
% Radiuses between the body and the satellite
R_earth = sqrt((sat.x - earth.x)^2 + (sat.y - earth.y)^2 +  (sat.z - earth.z)^2);
R_sun = sqrt((sun.x - sat.x)^2 + (sun.y - sat.y)^2 +  (sun.z - sat.z)^2);
R_moon = sqrt((moon.x - sat.x)^2 + (moon.y - sat.y)^2 +  (moon.z - sat.z)^2);
R_jupiter = sqrt((jupiter.x - sat.x)^2 + (jupiter.y - sat.y)^2 +  (jupiter.z - sat.z)^2);
R_venus = sqrt((venus.x - venus.x)^2 + (venus.y - sat.y)^2 +  (venus.z - sat.z)^2);
R_mars = sqrt((mars.x - sat.x)^2 + (mars.y - sat.y)^2 +  (mars.z - sat.z)^2);
R_saturn = sqrt((saturn.x - sat.x)^2 + (saturn.y - sat.y)^2 +  (saturn.z - sat.z)^2);
% 
% Radiuses between celestial bodies
R_earth_sun = sqrt((sun.x - earth.x)^2 + (sun.y - earth.y)^2 +  (sun.z - earth.z)^2);
R_earth_moon = sqrt((moon.x - earth.x)^2 + (moon.y - earth.y)^2 +  (moon.z - earth.z)^2);
R_earth_jupiter = sqrt((jupiter.x - earth.x)^2 + (jupiter.y - earth.y)^2 +  (jupiter.z - earth.z)^2);
R_earth_venus = sqrt((venus.x - earth.x)^2 + (venus.y - earth.y)^2 +  (venus.z - earth.z)^2);
R_earth_mars = sqrt((mars.x - earth.x)^2 + (mars.y - earth.y)^2 +  (mars.z - earth.z)^2);
R_earth_saturn = sqrt((saturn.x - earth.x)^2 + (saturn.y - earth.y)^2 +  (saturn.z - earth.z)^2);

% Earth is a primary body here

earth_influence = -((earth.GM*(sat.coords - earth.coords))/(R_earth)^3);
sun_influence = (sun.GM*(((sun.coords - sat.coords)/R_sun^3) -  ((sun.coords - earth.coords)/R_earth_sun^3)));
moon_influence = (moon.GM*(((moon.coords - sat.coords)/R_moon^3) -  ((moon.coords - earth.coords)/R_earth_moon^3)));
jupiter_influence = (jupiter.GM*(((jupiter.coords - sat.coords)/R_jupiter^3) -  ((jupiter.coords - earth.coords)/R_earth_jupiter^3)));
venus_influence = (venus.GM*(((venus.coords - sat.coords)/R_venus^3) -  ((venus.coords - earth.coords)/R_earth_venus^3)));
mars_influence = (mars.GM*(((mars.coords - sat.coords)/R_mars^3) -  ((mars.coords - earth.coords)/R_earth_mars^3)));
saturn_influence = (saturn.GM*(((saturn.coords - sat.coords)/R_saturn^3) -  ((saturn.coords - earth.coords)/R_earth_saturn^3)));


a_earth_sat =  earth_influence + sun_influence + moon_influence + jupiter_influence + venus_influence + mars_influence + saturn_influence;


% influences(:,1) = earth_influence;
% influences(:,2) = sun_influence;
% influences(:,3) = moon_influence;
% influences(:,4) = jupiter_influence;
% influences(:,5) = venus_influence;
% influences(:,6) = mars_influence;
% influences(:,7) = saturn_influence;

%Energy
% r_sun_earth = sqrt(sun.x^2 + sun.y^2 + sun.z^2); % earth coords = 0
% potential_energy = (G/2)*( (earth.mass*(sun.mass/r_sun_earth)) + (sun.mass*(earth.mass/r_sun_earth)) );
% 
% kinetic_energy = sun.mass*(sqrt(sun.vx^2 + sun.vy^2 + sun.vz^2))^2;
% 
% epoch_energy = kinetic_energy - potential_energy;
% 
% epoch_energy = epoch_energy - Initial_energy;
% 
% energy(1,epoch) = kinetic_energy;
% energy(2,epoch) = potential_energy;
% energy(3,epoch) = epoch_energy;

%% Deal with summations
% b - bodies
b = [sat, earth, sun, moon, jupiter, venus, mars, saturn]; % N = 8
N = length(b);

% Total kinetic energy
% kinetic = zeros(1,N);
% for k = 1:length(N)
%     kinetic_for_one_body = (b(k).mass*(sqrt(b(k).vx^2 + b(k).vy^2 + b(k).vz^2))^2)/2;
%     kinetic(1,k) = kinetic_for_one_body;
% end
% total_kinetic = sum(kinetic);
% 
% % Total potential energy
% potential = zeros(1,N);
% for i = 1:length(N)
%     second_terms = zeros(1,N);
%     for j = 1:length(N) %-1
%         if j ~= i
%             r_ij = sqrt((b(j).x - b(i).x)^2 + (b(j).y - b(i).y)^2 + (b(j).z - b(i).z)^2);
%             second_terms(1,j) = b(j).mass / r_ij;
%         else
%             second_terms(1,j) = 0;
%         end
%     end
%     total_second_terms = sum(second_terms);
%     
%     potential_for_one_body = b(i).mass*total_second_terms;
%     potential(1,i) = potential_for_one_body;
% end
% 
% total_potential = (G/2)*sum(potential);
% 
% 
% epoch_energy = total_kinetic - total_potential;

% [epoch_energy, total_kinetic, total_potential] = calculate_energy(b);
% 
% epoch_kinetic = total_kinetic - Initial_kinetic;
% epoch_potential = total_potential - Initial_potential;
% 
% epoch_energy = epoch_energy - Initial_energy;
% 
% 
% energy(1,epoch) = epoch_kinetic;
% energy(2,epoch) = epoch_potential;
% energy(3,epoch) = epoch_energy;

%% SOLAR PRESSURE
A = 264; % m2
refl = 0.5; % -
Crefl = 1+refl; % -
m = 6500; %kg
AU = 149*10^6; %km
radius = 151.5*10^6; %km

if pressure_on == 1
solar_a = -(4.56*10^-6)*Crefl*A/m*(radius/radius^3)*AU^2;
elseif pressure_on == 0
solar_a = 0;
end

%% Total Acceleration for a given planet
yp=zeros(6,1);
yp(1)=y0(4);
yp(2)=y0(5);
yp(3)=y0(6);

yp(4)= a_earth_sat(1);
yp(5)= a_earth_sat(2);
yp(6)= a_earth_sat(3);

% yp(4)=sun_influence(1) + total(1) + solar_a;
% yp(5)=sun_influence(2) + total(2) + solar_a;
% yp(6)=sun_influence(3) + total(3) + solar_a;



end
