function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


    % Final Time
    final_epoch = 10000.140558185330e+006;

  % Set initial state
R0 = [1485367.12031441;-657428.270129797;-680699.461857011];
V0 = [-0.0177617504518854;-0.000775359754553838;-0.00145931403670204];
init_epoch = 9.669336658152890e+08;

    phi0 = reshape(eye(6), 36, 1);
    
    init_state = [R0; V0+dV; phi0];
    
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
    [t, y0state,output_state, last_point_in_E] = simple_rkv89emb_maneuvers(@simplified_force_model, init_epoch, init_state);
    
    %desiredX = abs(y0state(1)) - 3.6361e+05; % vary this value depending on the maneuver 
    
    %scalar_pos = norm([-1.446510599792428e+06;9.382527780557508e+05;2.085037199136087e+04]); % for dV3
    %desiredX = abs(norm([last_point_in_E(1);last_point_in_E(2);last_point_in_E(3)])) - scalar_pos;% dV1: 5.62e+05; % EME
    % Now dV3
    
%     desX = abs(last_point_in_E(1)) - 1.44651e+06;%dv3 1.44e+06; % Reference coordinates of simplified + SRP
%     desY = abs(last_point_in_E(2)) - 9.382527e+05;%dv3 9.38e+05;
%     desZ = abs(last_point_in_E(3)) - 2.08503e+04;%dv3 2.08e+04;
    
%     desX = last_point_in_E(1);%dv3 1.44e+06; % Reference coordinates of simplified + SRP
%     desY = last_point_in_E(2);%dv3 9.38e+05;
%     desZ = last_point_in_E(3);%dv3 2.08e+04;
% 
% 
%     tolerance = 45000;
%     
%     % dV2 norm
%     dvnorm = norm([5.573172181394189e+05;7.817329418631776e+05; 6.197496009955091e+05]);
% 
%     desiredX = dvnorm - norm([desX;desY;desZ]);
%     
%       if (desiredX < tolerance && desiredX > -tolerance)
%        desiredX = 0;
%      end
%     
    %desiredX = [desX;desY;desZ];

    
    %ystar = [y0state(4);y0state(6);desiredX];
    

    ystar = [y0state(4);y0state(6)];
    
    % simplified version. When only X coordinate need to be adjusted
   % desiredX = abs(y0state(1)) - 1.9e+05;
    %ystar = desiredX;
%     
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
      phi = reshape(solution.ye(7:end), 6, 6);
      phi = xform*phi*xform^(-1);
      phi = reshape(phi, 36,1);
      L2state = [L2state; phi];
                 
     %disp(ode45_l2_state);
     ystar = [L2state(4);L2state(6)];
    end
    
    if ode_113

   
    
     options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
     solution = ode113(@full_force_model,[init_epoch final_epoch],init_state,options);
       
    L2point = cspice_spkezr('392', solution.xe, 'J2000', 'NONE', '399');
      conv_state = solution.ye;
      conv_state(1:6) = solution.ye(1:6) - L2point;
      xform = cspice_sxform('J2000','L2CENTERED', solution.xe);
      L2state = xform*conv_state(1:6);
      phi = reshape(solution.ye(7:end), 6, 6);
      phi = xform*phi*xform^(-1);
      phi = reshape(phi, 36,1);
      L2state = [L2state; phi];
                 
     %disp(ode45_l2_state);
     ystar = [L2state(4);L2state(6)];
     
     disp(ystar);
    end
    
end

