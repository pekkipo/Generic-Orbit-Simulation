function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


    % Final Time
    final_epoch = 10000.140558185330e+006;

  % Set initial state

 

R0 = [-1.445778814591190e+06;9.392434662856446e+05;2.375365526234356e+04];
 V0 = [0.006242577552817;-0.001461621210794;3.374369574140406e-04];
% Initial Time
init_epoch = 9.824219618479762e+08;

    phi0 = reshape(eye(6), 36, 1);
    
    init_state = [R0; V0+dV; phi0];
    
    % choose integrator - ONE
rkv_89emb = true;
rkv_89 = false;
ode_45 = false;
ode_113 = false;
ode_87 = false;
abm_8 = false;

    %% RKV89
    if rkv_89emb 
    [t, y0state] = full_rkv89emb_maneuvers(@full_force_model, init_epoch, init_state);
    
    %desiredX = abs(y0state(1)) - 3.+05; % vary this value depending on the maneuver 
    %ystar = [desiredX;y0state(4);y0state(6)];

    ystar = [y0state(4);y0state(6)];
    
    % simplified version. When only X coordinate need to be adjusted
   % desiredX = abs(y0state(1)) - 1.9e+05;
    %ystar = desiredX;
%     
    disp(ystar);
    end
    %%
    
    %% ODE45
    if ode_45
    global ode45_l2_state;
    % Have to make it global as event_handler doesn't provide an
    % opportunity to output other values
    
     options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
     solution = ode45(@simplified_force_model_srp,[init_t final_epoch],init_state,options);
       
     disp(ode45_l2_state);
     ystar = [ode45_l2_state(4);ode45_l2_state(5)];
    end
end

