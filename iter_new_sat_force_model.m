function yp = iter_new_sat_force_model( t,y0, gm_planets, planets_ephemeris, epoch, pressure_on )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global influences;
global energy;
global Initial_energy;
global G;

% pressure_on: 1 - including solar pressure, 0 - without solar pressure
planets = planets_ephemeris;

GM_SUN = gm_planets(1,1);
GM_MERCURY = gm_planets(1,2);
GM_VENUS = gm_planets(1,3);
GM_EARTH = gm_planets(1,4);
GM_MARS = gm_planets(1,5);
GM_JUPITER = gm_planets(1,6);
GM_SATURN = gm_planets(1,7);
GM_URANUS = gm_planets(1,8);
GM_NEPTUNE = gm_planets(1,9);
GM_MOON = gm_planets(1,10);


%% Accelerations due to:
% This time I will consider 3 body system - Sun, Earth, Satellite
% Later I will add other sources
% GRAVITY

% y0 - satellite, rows: x y z vx vy vz

% SUN Nr 1
R_sun = sqrt((planets(1,1) - y0(1))^2 + (planets(2,1) - y0(2))^2 +  (planets(3,1) - y0(3))^2);
% EARTH Nr 4
R_earth = sqrt((y0(1) - planets(1,4))^2 + (y0(2) - planets(2,4))^2 +  (y0(3) - planets(3,4))^2);
% MOON Nr 10
R_moon = sqrt((planets(1,10) - y0(1))^2 + (planets(2,10) - y0(2))^2 +  (planets(3,10) - y0(3))^2);
% JUPITER Nr 6
R_jupiter = sqrt((planets(1,6) - y0(1))^2 + (planets(2,6) - y0(2))^2 +  (planets(3,6) - y0(3))^2);
% VENUS Nr 3
R_venus = sqrt((planets(1,3) - y0(1))^2 + (planets(2,3) - y0(2))^2 +  (planets(3,3) - y0(3))^2);
% MARS Nr 5
R_mars = sqrt((planets(1,5) - y0(1))^2 + (planets(2,5) - y0(2))^2 +  (planets(3,5) - y0(3))^2);
% SATURN Nr 7
R_saturn = sqrt((planets(1,7) - y0(1))^2 + (planets(2,7) - y0(2))^2 +  (planets(3,7) - y0(3))^2);
% 
R_earth_sun = sqrt((planets(1,1) - planets(1,4))^2 + (planets(2,1) - planets(2,4))^2 +  (planets(3,1) - planets(3,4))^2);
R_earth_moon = sqrt((planets(1,10) - planets(1,4))^2 + (planets(2,10) - planets(2,4))^2 +  (planets(3,10) - planets(3,4))^2);
R_earth_jupiter = sqrt((planets(1,6) - planets(1,4))^2 + (planets(2,6) - planets(2,4))^2 +  (planets(3,6) - planets(3,4))^2);
R_earth_venus = sqrt((planets(1,3) - planets(1,4))^2 + (planets(2,3) - planets(2,4))^2 +  (planets(3,3) - planets(3,4))^2);
R_earth_mars = sqrt((planets(1,5) - planets(1,4))^2 + (planets(2,5) - planets(2,4))^2 +  (planets(3,5) - planets(3,4))^2);
R_earth_saturn = sqrt((planets(1,7) - planets(1,4))^2 + (planets(2,7) - planets(2,4))^2 +  (planets(3,7) - planets(3,4))^2);

% Earth is a primary body here

earth_influence = -(GM_EARTH/(R_earth)^3)*(y0 - planets(:,4));
sun_influence = (GM_SUN*(((planets(:,1) - y0)/R_sun^3) -  ((planets(:,1) - planets(:,4))/R_earth_sun^3)));
moon_influence = (GM_MOON*(((planets(:,10) - y0)/R_moon^3) -  ((planets(:,10) - planets(:,4))/R_earth_moon^3)));
jupiter_influence = (GM_JUPITER*(((planets(:,6) - y0)/R_jupiter^3) -  ((planets(:,6) - planets(:,4))/R_earth_jupiter^3)));
venus_influence = (GM_VENUS*(((planets(:,3) - y0)/R_venus^3) -  ((planets(:,3) - planets(:,4))/R_earth_venus^3)));
mars_influence = (GM_MARS*(((planets(:,5) - y0)/R_mars^3) -  ((planets(:,5) - planets(:,4))/R_earth_mars^3)));
saturn_influence = (GM_SATURN*(((planets(:,7) - y0)/R_saturn^3) -  ((planets(:,7) - planets(:,4))/R_earth_saturn^3)));


a_earth_sat =  earth_influence + sun_influence + moon_influence + jupiter_influence + venus_influence + mars_influence + saturn_influence;


influences(:,1) = earth_influence;
influences(:,2) = sun_influence;
influences(:,3) = moon_influence;
influences(:,4) = jupiter_influence;
influences(:,5) = venus_influence;
influences(:,6) = mars_influence;
influences(:,7) = saturn_influence;

% Energy
Earth_mass = (GM_EARTH * 10^9)/G;
Sun_mass = (GM_SUN * 10^9)/G;
% Moon_mass  = (GM_EARTH * 10^9)/G;
% Jupiter_mass = (GM_EARTH * 10^9)/G;
% Venus_mass = (GM_EARTH * 10^9)/G;
% Mars_mass = (GM_EARTH * 10^9)/G;

r_sun_earth = sqrt(planets(1,1)^2 + planets(2,1)^2 + planets(3,1)^2); % earth coords = 0
potential_energy = (G/2)*( (Earth_mass*(Sun_mass/r_sun_earth)) + (Sun_mass*(Earth_mass/r_sun_earth)) );
kinetic_energy = Sun_mass*(sqrt(planets(4,1)^2 + planets(5,1)^2 + planets(6,1)^2))^2;

epoch_energy = kinetic_energy - potential_energy;

epoch_energy = epoch_energy - Initial_energy;

energy(1,epoch) = kinetic_energy;
energy(2,epoch) = potential_energy;
energy(3,epoch) = epoch_energy;


% SOLAR PRESSURE
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