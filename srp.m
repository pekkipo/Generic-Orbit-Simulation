function solar_a = srp(type)
% SRP solar radiation pressure
    % Simplified, considers perependicular light

% type = 0 - L2env paper formula
% type = 1 - Montebruck formula

% Constants
A = 264; % m2
refl = 0.9; % -
Crefl = 1+refl; % -
m = 6000; %kg
AU = 149*10^6; %km
radius = 151.5*10^6; %km
c = 299792458; % m/s
P0 = (4.56*10^-6); % N/m2 = kg/m*s2
fluxatl2 = 1340; % w/m2 = kg/s3

if type == 0
solar_a = (fluxatl2*A*Crefl)/(m*c);
else
solar_a = -P0*Crefl*A/m*(radius/radius^3)*(AU^2);
end
end

