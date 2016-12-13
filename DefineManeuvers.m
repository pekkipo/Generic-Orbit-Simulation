
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

% Set initial state

R0 = [5.573172181394189e+05;7.817329418631776e+05; 6.197496009955091e+05];
V0 = [-0.551355109432290; 0.282329029512973;0.123359752969536];
% Initial Time
initial_time = 9.748534057646916e+08;
% Final Time
final_time = 10000.140558185330e+006;


phi0 = reshape(eye(6), 36, 1);
init_state = [R0; V0; phi0];


global G;
G = 6.673e-20;
global L2frame;
L2frame = true;
global checkrkv89_emb
checkrkv89_emb = false;
global rkv89emb_lastpiece;
rkv89emb_lastpiece = false;
% Initial guess
dV = [ 26.5303022604513e-003;1.23533361008047e-003;4.00576011107553e-003];

%options = optimoptions('fsolve','TolFun', 1e-4, 'TolX', 1e-4);
deltaV = fsolve(@evaluate_V_test, dV);
disp(deltaV);

Init_state = init_state;
Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;


% choose integrator - ONE
rkv_89emb = true;
rkv_89 = false;
ode_45 = false;
ode_113 = false;
ode_87 = false;
abm_8 = false;

% CHECK THE FORCE MODEL USED

%% RKV89
if rkv_89emb 
[t, y0state, output_state, y0state_E] = full_rkv89emb_maneuvers(@full_force_model, initial_time , Init_state);
end
%% ODE87
if ode_87
    [t, y0state, output_state, y0state_E] = ode87(@full_force_model, initial_time , Init_state);
end

%% ODE45
if ode_45
 options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
 solution = ode45(@simplified_force_model_srp,[initial_time final_time],Init_state,options);
 epochs = solution.x;
 orbit = solution.y;
 output_state = EcenToL2frame( orbit, epochs );
end
% Init_state is the new value 
% y0state is the state, from which the next integration should start


% Graphical check of the orbit part
figure
hold on
plot3(output_state(1,:),output_state(2,:),output_state(3,:),'r','LineWidth',2)

