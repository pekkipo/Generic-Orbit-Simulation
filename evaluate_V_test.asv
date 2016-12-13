function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


    % Final Time
    final_epoch = 10000.140558185330e+006;

  % Set initial state

    R0 = [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005];  
 V0 = [5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];
% Initial Time
init_epoch = 958.910668311133e+006;


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
    ystar = [y0state(4);y0state(6);];
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

