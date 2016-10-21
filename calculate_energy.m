function [epoch_energy, total_kinetic, total_potential] = calculate_energy(b)
%Calculate energy for one epoch
%   Takes in a vector of structures describing bodies and calculate total
%   mechanical energy for one point

%global G;
G = 6.67e-20; % Changing from -17 to -20 makes kinetic and potenital of the order
N = length(b);

% Total kinetic energy
kinetic = zeros(1,N);
for k = 1:N
    kinetic_for_one_body = (b(k).mass*(sqrt(b(k).vx^2 + b(k).vy^2 + b(k).vz^2))^2)/2;
    kinetic(1,k) = kinetic_for_one_body;
end
%disp(kinetic)
total_kinetic = sum(kinetic);

% Total potential energy
potential = zeros(1,N);
second_terms = zeros(1,N);
    for i = 1:N
    %second_terms = zeros(1,N);
     for j = 1:N %-1
         if j ~= i
            r_ij = sqrt((b(j).x - b(i).x)^2 + (b(j).y - b(i).y)^2 + (b(j).z - b(i).z)^2);
            second_terms(1,j) = b(j).mass / r_ij;
            %disp(second_terms(1,j));
         else
            second_terms(1,j) = 0;
            %disp(second_terms(1,j));
         end
     end
    total_second_terms = sum(second_terms);
    
    potential_for_one_body = b(i).mass*total_second_terms;
    potential(1,i) = potential_for_one_body;
    end
disp(potential)
total_potential = (G/2)*sum(potential);


epoch_energy = total_kinetic - total_potential;

end