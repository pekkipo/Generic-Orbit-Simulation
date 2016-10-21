function yp = verification_force_model( t,y0  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

observer = 'SUN';
planets_simplified = {'EARTH', 'SUN', 'MOON';'EARTH', 'SUN', '301'};
[earth, sun, moon] = simplified_create_structure( planets_simplified, t, observer);

G=6.67e-20;

% GRAVITY
% Earth is a primary body here!

% Radiuses between the body and the satellite
R_e_moon = sqrt((moon.x - earth.x)^2 + (moon.y - earth.y)^2 +  (moon.z - earth.z)^2);
R_sun_moon = sqrt((sun.x - moon.x)^2 + (sun.y - moon.y)^2 +  (sun.z - moon.z)^2);


% Radiuses between celestial bodies
R_earth_sun = sqrt((sun.x - earth.x)^2 + (sun.y - earth.y)^2 +  (sun.z - earth.z)^2);

% Earth is a primary body here

earth_influence = -(G*(earth.mass + moon.mass)/(R_e_moon)^3)*(moon.coords - earth.coords);
sun_influence = (sun.GM*(((sun.coords - moon.coords)/R_sun_moon^3) -  ((sun.coords - earth.coords)/R_earth_sun^3)));

a_moon =  earth_influence + sun_influence;

%% Total Acceleration for a given planet
yp=zeros(6,1);
yp(1)=y0(4);
yp(2)=y0(5);
yp(3)=y0(6);

yp(4)= a_moon(1);
yp(5)= a_moon(2);
yp(6)= a_moon(3);

end
