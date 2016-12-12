function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


    % Final Time
    final_epoch = 10000.140558185330e+006;

  % Set initial state

     R0 = [5.470091051396191e+05;7.746809140774258e+05;6.136761030644372e+05];
     V0 =  [-0.561328688405743;0.285354232749995;0.124723965335411];
% Initial Time
    init_epoch = 9.748741253169044e+08;


    phi0 = reshape(eye(6), 36, 1);
    
    init_state = [R0; V0+dV; phi0];

    [t, y0state] = full_rkv89emb_maneuvers(@full_force_model, init_epoch, init_state);
    
    
    ystar = [y0state(4);y0state(6);];
    
    disp(ystar);

end

