function [ystar] = new_evaluate_V( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


  % Set initial state
R0 = [-1447522.71235131;938357.823129112;22691.6506246949];
V0 = [0.00640971096544824;-0.000955550545351669;0.000543973359516699];
init_epoch = 9.824283500767865e+08;
final_epoch = 9.902741246238446e+08;

%global output_state;

    phi0 = reshape(eye(6), 36, 1);
    init_state = [R0; V0; phi0];
    
    
    [t, y0state] = ode87_test_y(@full_force_model, [init_epoch final_epoch], init_state);

    ystar = [y0state(2);y0state(4);y0state(6)];
    
    
      
    disp(ystar);
    disp('maneuver');
    disp(dV);
 
end


% 
%     [t1, y0state1, output1, last_point_in_E1] = ode87_test_y(@full_force_model, [init_epoch1 final_epoch1], init_state1);
%     output_state = [output_state, output1];
%     init_state2 = [last_point_in_E1(1:3); last_point_in_E1(4:6)+dV2; phi0];
%     
%     [t2, y0state2, output2, last_point_in_E2] = ode87_test_y(@full_force_model, [init_epoch2 final_epoch2], init_state2);
%     output_state = [output_state, output2];
%     init_state3 = [last_point_in_E2(1:3); last_point_in_E2(4:6)+dV3; phi0];
%     
%     [t3, y0state, output3, last_point_in_E2] = ode87_test_y(@full_force_model, [init_epoch3 final_epoch3], init_state3);
%     output_state = [output_state, output3];
%     
%     init_epoch1 = 9.824283500767865e+08;
% final_epoch1 = init_epoch1 + 2.61525818235270e+006;
% 
% init_epoch2 = final_epoch1;
% final_epoch2 = init_epoch2 + 2.61525818235270e+006; 
% 
% init_epoch3 = final_epoch2;
% final_epoch3 = init_epoch3 + 2.61525818235270e+006; 

