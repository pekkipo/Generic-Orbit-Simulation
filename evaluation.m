function [ystar] = evaluation(dV)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global R0;
global V0;
global start_time;
final_epoch = 1000000.140558185330e+006;


phi0 = reshape(eye(6), 36, 1);
init_state = [R0; V0+dV; phi0];

% BE SURE THAT THE FORCE MODEL USED HERE IS THE SAME AS THE ONE IN THE
% INTEGRATION

[t, y0state] = rkv89emb_maneuvers(@full_force_model, start_time, init_state);
    
    ystar = [y0state(4);y0state(6)];

    disp(ystar);

end

