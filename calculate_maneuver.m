function [ deltaV ] = calculate_maneuver(init_guess)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Set initial state

    % Initial guess
    dV = init_guess;

    deltaV = fsolve(@evaluation, dV);
    disp(deltaV);

end

