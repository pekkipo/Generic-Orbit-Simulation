function [ystar] = NR_evaluate_V( dV )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    


  % Set initial state
R0 = [1441867.00000083;-714815.699000797;-698432.680000360];
V0 = [-0.00180416368562985;-0.00778595175744212;-0.00264024278703223];
init_epoch = 9.982927336197942e+08;
final_epoch = 1.006143146045806e+09;

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
       
    desX = stateE(1) - 6.198677e+05;
    desY = stateE(2) - 7.5667763e+05;
    desZ = stateE(3) - 6.068734737121682e+05;
    desVX = stateE(4) + 0.517107765090152;
    desVY = stateE(5) - 0.320457712270583;
    desVZ = stateE(6) - 0.139325366614006;
    
     -51.9914783815620e+000
   -3.53528859175276e+000
   -3.50182675651740e+000
   -19.2677041412903e-006
   -44.7378293937928e-006
   -812.276175737869e-009
    
    
    ystar = [desX;desY;desZ;desVX;desVY;desVZ];
    disp(ystar);
    disp('maneuver');
    disp(dV);
    disp('t');
    disp(t(end));
    end

end

