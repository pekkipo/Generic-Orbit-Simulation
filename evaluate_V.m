function [ gv, d_gv ] = evaluate_V( V0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    gv = zeros(2,1);
    d_gv = zeros(2,1);
    
    % for now define manually
    %init_state = [1.120509091323728e+06; 1.885186194584079e+03; 2.569197970885182e+05; ...
    %    -0.011126300564490; 0.394098240322462; 0.001333729434840];
    init_state = V0;
    init_epoch = 9.747626128418571e+08;
    final_epoch = 9.903366994711639e+08;
    

    phi0 = [1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1];
    %eye(6) would do the same :)
    init_state = [init_state; phi0];
    
    [epochs, orbit, lastState_E] = rkv89emb(@force_model_maneuvers, [init_epoch final_epoch], init_state, 2, true);
    
    last_ind = length(orbit);
    
    gv = [orbit(4, last_ind); orbit(6, last_ind)];
    d_gv = [orbit(28, last_ind); orbit(42, last_ind)];
    
end

