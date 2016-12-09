function [ystar] = evaluate_V_test( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    % The very beginning
    % Initial Time
    init_epoch = 9.669336468097093e+08;
    % Final Time
    final_epoch = 999.140558185330e+006;


%  1   R0 = [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005];  
%     V0 = [5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];
 R0 = [1.485359168304580e+06;  -6.574308768747580e+05; -6.807000270951350e+05];  
 V0 = [-0.017761982577070;  -7.753583934581932e-04;  -0.001459319442820];
    phi0 = reshape(eye(6), 36, 1);
    
    init_state = [R0; V0+dV; phi0];

    [t, y0state] = rkv89emb_maneuvers(@force_model_maneuvers, [init_epoch final_epoch], init_state);
    
    ystar = [y0state(4);y0state(6)];
    %reshaped_y = reshape(y0state(7:end),6,6);
    %Phi = [reshaped_y(4,4:6);reshaped_y(6,4:6)];

    disp(ystar);

end

