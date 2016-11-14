function solar_a = srp(type, earth, sun)
% SRP solar radiation pressure
    % Simplified, considers perependicular light

% type = 0 - L2env paper formula
% type = 1 - Montebruck formula
% USE - 1, Montenbruck p 75, equation 3.75 . r(.) - geocentrical position of the sun

% Constants
A = 264; % m2
refl = 0.5; % -
Crefl = 1+refl; % -
m = 6500; %kg
% AU would change with time, so better do this:
AU = sqrt((earth.x - sun.x)^2 + (earth.y - sun.y)^2 + (earth.z - sun.z)^2);
%c = 299792458; % m/s
c = 299792; % km/s
P0 = -(4.56*10^(-6)); % N/m2 = kg/m*s2


% GET INFO ABOUT THE sun
 r_vector = sun.coords;
 r3 = (sqrt((sun.x)^2 + (sun.y)^2 +  (sun.z)^2))^3; 
 
% unit_vector = [sat.vx; sat.vy; sat.vz]/(sqrt(sat.vx^2 + sat.vy^2 + sat.vz^2));

if type == 0
    solar_a = (1340*A*Crefl)/(m*c); % 1340 flux at L2 in w/m2 = kg/s3
   % solar_a = solar_a*unit_vector;
elseif type == 1
    % divide by 10^3 as P0 value in k/m*s2, while I need kilometers
    solar_a = (P0/(10^3))*Crefl*(A/m)*(r_vector/r3)*(AU^2);
    % Works! Kinda..now the difference in orbits is quite big
end
end

