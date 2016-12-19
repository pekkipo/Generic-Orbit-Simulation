function [ystar] = NR_evaluate_V( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


  % Set initial state
R0 = [-6.05487505154821e+005;-1.01012602798113e+006;-1.48976351540980e+005];
V0 = [511.424488127811e-003;-288.436177303014e-003;-126.196496781444e-003];
init_epoch = 9.902741246238446e+08;
final_epoch = 9.982927336197942e+08;

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
       
    desX = stateE(1) - 1.441867e+06;
    desY = stateE(2) + 7.14815699e+05;
    desZ = stateE(3) + 6.9843268e+05;
    desVX = stateE(4) + 0.001817314659897;
    desVY = stateE(5) + 0.007755542122498;
    desVZ = stateE(6) + 0.002671479865356;
    
    ystar = [desX;desY;desZ;desVX;desVY;desVZ];
    disp(ystar);
    disp('maneuver');
    disp(dV);
    disp('t');
    disp(t(end));
    end

end

