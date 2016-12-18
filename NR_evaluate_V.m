function [ystar] = NR_evaluate_V( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


  % Set initial state
R0 = [1.48529760318229e+006;-6.57082469514976e+005;-6.80554647691743e+005];  
V0 = [-17.6089143212389e-003;-790.596318730594e-006;-1.41794858234358e-003];
init_epoch = 9.669342245618324e+08;
final_epoch = 9.747786247088768e+08;

init_state = [R0; V0+dV];
    

rkv_89emb = false;
ode_87 = true;



    %% RKV89
    if rkv_89emb 
    [t, y0state, output, stateE] = full_rkv89emb_maneuvers(@full_force_model, init_epoch, init_state);
    
    ystar = [y0state(4);y0state(6)];
    disp(ystar);
    disp('maneuver');
    disp(dV);
    end
    
    %% ODE87
    if ode_87
    [t, y0state, output, stateE] = ode87_test_y(@full_force_model, [init_epoch final_epoch], init_state); 
       
    desX = stateE(1) - 5.7180234e+05;
    desY = stateE(2) - 7.745195e+05;
    desZ = stateE(3) - 6.1742255e+05;
    desVX = stateE(4) + 0.542851805921651;
    desVY = stateE(5) - 0.286090687941076;
    desVZ = stateE(6) - 0.124778406579352;
    
    ystar = [desX;desY;desZ;desVX;desVY;desVZ];
    disp(ystar);
    disp('maneuver');
    disp(dV);
    disp('t');
    disp(t(end));
    end

end

