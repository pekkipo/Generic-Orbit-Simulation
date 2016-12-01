function [ deltaV, isFound ] = MyNewtonRaphson(dV, tstar, ystar )
%[deltaV, isFound] = NewtonRaphson(dV, epochs(length(epochs)),orbit_rkv89_emb(length(orbit_rkv89_emb)));
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        % 3 independent variables deltaVX, deltaVY, deltaVZ
        % 2 independent variables Vx and Vz which are supposed to be zero
        % +- tolerance
        
        tolerance = e-8;

        dV = [dV(4);dV(5);dV(6)];
        Vx = ystar(4);
        Vz = ystar(6);
        derivatives = force_model(tstar, ystar);
        % ystar(4) - Vx of the state at which y=0, ystar(6) - Vz at that
        % epoch
        J = [ystar(4)/derivatives(4), ystar(4)/derivatives(5), ystar(4)/derivatives(6); ystar(6)/derivatives(4), ystar(6)/derivatives(5), ystar(6)/derivatives(6)];
        J = pinv(J);
        
        deltaV = dV - J*([Vx Vz] - [0 0]);
        
        l_bound = 0 - tolerance;
        r_bound = 0 + tolerance;
        
        if (Vx > l_bound && Vx < r_bound) && (Vz > l_bound && Vz < r_bound)
            isFound = true;
        end
end

