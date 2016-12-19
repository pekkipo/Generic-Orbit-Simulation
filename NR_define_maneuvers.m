
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
R0 = [-6.05487505154821e+005;-1.01012602798113e+006;-1.48976351540980e+005];
V0 = [511.424488127811e-003;-288.436177303014e-003;-126.196496781444e-003];
init_epoch = 9.902741246238446e+08;
final_epoch = 9.982927336197942e+08;

init_state = [R0; V0];

% Initial guess
dV = [-0.001625936670348; -0.003125208256016; -0.008088501084076];

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

