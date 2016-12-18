
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

global G;
G = 6.673e-20;
global L2frame;
L2frame = true;
global rkv89emb_lastpiece;
rkv89emb_lastpiece = false;
global RKV_89_emb_check
RKV_89_emb_check = false;

% Set initial state

R0 = [1.48529760318229e+006;-6.57082469514976e+005;-6.80554647691743e+005];  
V0 = [-17.6089143212389e-003;-790.596318730594e-006;-1.41794858234358e-003];
init_epoch = 9.669342245618324e+08;
final_epoch = 9.747786247088768e+08;
init_state = [R0; V0];

% Initial guess
dV = [-7.803777280688135e-04; 0.001854569833090; -0.007247538179753];

  %deltaV = fsolve(@evaluate_V_test, dV);
  
  [deltaV, resnorm, F, exitflag, output, jacob] = newtonraphson(@NR_evaluate_V, dV);
  
    disp(deltaV);
    
    Init_state = init_state;
    Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;

% choose integrator - ONE
rkv_89emb = false;
ode_87 = true;


% CHECK THE FORCE MODEL USED

%% RKV89
if rkv_89emb 
[t, y0state, output_state, y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_epoch , Init_state);
end
%% ODE87
if ode_87
    [t, y0state, output_state, y0state_E] = ode87_test_y(@full_force_model, [init_epoch final_epoch] , Init_state);
end

% Graphical check of the orbit part
figure
hold on
plot3(output_state(1,:),output_state(2,:),output_state(3,:),'r','LineWidth',2)

