function [ deltaV, isFound ] = MyNewtonRaphson(dV, tstar, ystar )
%[deltaV, isFound] = NewtonRaphson(dV, epochs(length(epochs)),orbit_rkv89_emb(length(orbit_rkv89_emb)));
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        % 3 independent variables deltaVX, deltaVY, deltaVZ
        % 2 dependent variables Vx and Vz which are supposed to be zero
        % +- tolerance
        
        tolerance = 1e-8;

        dV = [dV(4);dV(5);dV(6)];
        Vx = ystar(4);
        Vz = ystar(6);
        d_Vx = 0.00001;
        d_Vy = 0.00001;
        d_Vz = 0.00001;
        second_y = ystar + [0;0;0;d_Vx;d_Vy;d_Vz];
        derivatives = force_model(tstar, ystar);
        first_state = ystar;
        first_dV = [first_state(4);first_state(5);first_state(6)];
        derivatives2 = force_model(tstar, second_y);
        second_state = ystar + [0;0;0;derivatives2(4);derivatives2(5);derivatives2(6)];
        second_dV = [second_state(4);second_state(5);second_state(6)];
        % ystar(4) - Vx of the state at which y=0, ystar(6) - Vz at that

            J11 = (second_dV(1)  - first_dV(1))/d_Vx;
            J12 = (second_dV(1)  - first_dV(1))/d_Vy;
            J13 = (second_dV(1)  - first_dV(1))/d_Vz;
            J21 = (second_dV(3)  - first_dV(3))/d_Vx;
            J22 = (second_dV(3)  - first_dV(3))/d_Vy;
            J23 = (second_dV(3)  - first_dV(3))/d_Vz;


        J = [J11, J12, J13; J21, J22, J23];
        %J = pinv(J);
        J=J';
        
        state_with_old_dv = [Vx;Vz];
        desired_state = [0;0];
        deltaV = dV - J*(state_with_old_dv - desired_state);
        %deltaV = dV - J;
        
        l_bound = 0 - tolerance;
        r_bound = 0 + tolerance;
        
        if (Vx > l_bound && Vx < r_bound) && (Vz > l_bound && Vz < r_bound)
            isFound = true;
            
        else
            isFound = false;
        end
        
        disp(deltaV);
        disp(Vx);
        disp(Vz);
        
        deltaV = [0;0;0;deltaV(1);deltaV(2);deltaV(3)];
end

%         J11 = ((Vx+derivatives(4)) - Vx)/derivatives(4);
%         J12 = ((Vx+derivatives(4)) - Vx)/derivatives(5);
%         J13 = ((Vx+derivatives(4)) - Vx)/derivatives(6);
%         J21 = ((Vz+derivatives(6)) - Vz)/derivatives(4);
%         J22 = ((Vz+derivatives(6)) - Vz)/derivatives(5);
%         J23 = ((Vz+derivatives(6)) - Vz)/derivatives(6);
