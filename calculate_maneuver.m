function [ deltaV ] = calculate_maneuver()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Set initial state

    % Initial guess
    dV = [0.001; -0.003; -0.005];

    deltaV = fsolve(@evaluation, dV);
    disp(deltaV);

end

