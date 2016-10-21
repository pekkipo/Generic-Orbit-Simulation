function [total, kinetic, potential] = simple_energy(b, earth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

kinetic = (b.mass/2*(sqrt(b.vx^2 + b.vy^2 + b.vz^2))^2);
potential = earth.GM*b.mass/sqrt((earth.x - b.x)^2 + (earth.y - b.y)^2 + (earth.z - b.z)^2);

total = kinetic - potential;

end

