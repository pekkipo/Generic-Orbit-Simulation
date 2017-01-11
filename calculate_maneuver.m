function [ deltaV ] = calculate_maneuver(init_guess)
    % Initial guess
    dV = init_guess;
    deltaV = fsolve(@evaluation, dV);
    %deltaV = newtonraphson(@evaluation, dV);
    disp(deltaV);
end

