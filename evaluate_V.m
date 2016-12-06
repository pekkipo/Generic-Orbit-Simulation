function [gv, d_gv] = evaluate_V( V0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    % The very beginning
    init_epoch = 9.589106683111327e+08; 
    %9.747668814452398e+08;
    final_epoch = 9.912611163179582e+08;%9.903366994711639e+08;

    % MY NEW VALUES
    %positions = [5.772997863348321e+05; 7.789010371518550e+05; 6.176535354084261e+05];
    % Very beginning
    positions = [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005];
    phi0 = [1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1];
    
    init_state = [positions; V0; phi0];

    y0state = rkv89emb_maneuvers(@force_model_maneuvers, [init_epoch final_epoch], init_state, 2, true);
    disp(y0state([4,6]));
    
    %last_ind = length(orbit);
    
    gv = [y0state(4);y0state(6)];
    reshaped_y = reshape(y0state(7:end),6,6);
    %d_gv = [y0state(28)'; y0state(42)'];
    %d_gv = [y0state(25:30); y0state(37:42)];
    d_gv = [reshaped_y(4,4:6);reshaped_y(6,4:6)];
    
    % different approach
%     M = y0state(7:42);
%     M = reshape(M, 6,6);
%     gv = M*init_state(1:6) - [y0state(1);y0state(2);y0state(3);0;y0state(5);0] + y0state(1:6); 
    
    
end

