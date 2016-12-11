function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    

    % Initial Time
    init_epoch = 9.982927336197942e+08;
    % Final Time
    final_epoch = 10000.140558185330e+006;

    R0 = [1.441868581740170e+06; -7.148156992129891e+05; -6.984326839322323e+05];
    V0 = [-0.001817314659897;-0.007755542122498;-0.002671479865356];
    phi0 = reshape(eye(6), 36, 1);
    
    init_state = [R0; V0+dV; phi0];

    [t, y0state] = rkv89emb_maneuvers(@simplified_force_model_srp, init_epoch, init_state);
    
    ystar = [y0state(4);y0state(6)];

    disp(ystar);

end

