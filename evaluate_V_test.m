function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


    % Final Time
    final_epoch = 10000.140558185330e+006;

  % Set initial state
R0 = [1441893.99780414;-714789.877498610;-698421.856243787];
V0 = [-0.00181347353676022;-0.00775519839178187;-0.00267104564885495];
init_epoch = 9.982928485909765e+08;


init_state = [R0; V0+dV];
    
    % choose integrator - ONE
rkv_89emb = true;
rkv_89 = false;
ode_45 = false;
ode_113 = false;
ode_87 = false;
abm_8 = false;

% Last arugment for integrator: 
% true - if intermediate maneuver after 1.5
% false - if usual stopping condition

    %% RKV89
    if rkv_89emb 
    [t, y0state] = simple_rkv89emb_maneuvers(@simplified_force_model, init_epoch, init_state);
    
    ystar = [y0state(4);y0state(6)];
    
    disp(ystar);
    disp('maneuver');
    disp(dV);
    end
    %% ODE87
    if ode_87
       [t, y0state, output_state, y0state_E] = full_ode87(@full_force_model, [init_epoch final_epoch], init_state); 
        ystar = [y0state(4);y0state(6)];
         disp(ystar);
    disp('maneuver');
    disp(dV);
    end
    
    
    
    %% ODE45
    if ode_45

     options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
     solution = ode45(@full_force_model,[init_epoch final_epoch],init_state,options);
     
      L2point = cspice_spkezr('392', solution.xe, 'J2000', 'NONE', '399');
      conv_state = solution.ye;
      conv_state(1:6) = solution.ye(1:6) - L2point;
      xform = cspice_sxform('J2000','L2CENTERED', solution.xe);
      L2state = xform*conv_state(1:6);
            

     ystar = [L2state(4);L2state(6)];
     
    disp(ystar);
    disp('maneuver');
    disp(dV);
    end
    
    if ode_113

   
    
     options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
     solution = ode113(@full_force_model,[init_epoch final_epoch],init_state,options);
       
    L2point = cspice_spkezr('392', solution.xe, 'J2000', 'NONE', '399');
      conv_state = solution.ye;
      conv_state(1:6) = solution.ye(1:6) - L2point;
      xform = cspice_sxform('J2000','L2CENTERED', solution.xe);
      L2state = xform*conv_state(1:6);

    
                 
     %disp(ode45_l2_state);
     ystar = [L2state(4);L2state(6)];
     
    disp(ystar);
    disp('maneuver');
    disp(dV);
    end
    
end

