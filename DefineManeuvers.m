
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

% Set initial state

R0 = [1485367.12031441;-657428.270129797;-680699.461857011];
V0 = [-0.0177617504518854;-0.000775359754553838;-0.00145931403670204];
init_epoch = 9.669336658152890e+08;
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
global RKV_89_emb_check
RKV_89_emb_check = false;
% Initial guess
dV = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 

%options = optimoptions('fsolve','TolFun', 1e-4, 'TolX', 1e-4);

% Use corrector
% usecor = true;
% if usecor
%     deltaV = fsolve(@evaluate_V_test, dV);
%     disp(deltaV);
%     Init_state = init_state;
%     Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;
% else
%     Init_state = init_state;
%     deltaV = [13.2530134946924e-003; -16.2338801932272e-003; 4.06482679813973e-003];
%     Init_state(4:6,:) = Init_state(4:6,:) + deltaV;   
% end
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
[t, y0state, output_state, y0state_E] = simple_rkv89emb_maneuvers(@simplified_force_model, init_epoch , Init_state);
end
%% ODE87
if ode_87
    [t, y0state, output_state, y0state_E] = full_ode87(@full_force_model, [init_epoch final_time] , Init_state);
end

%% ODE45
if ode_45
 options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
 solution = ode45(@full_force_model,[init_epoch final_time],Init_state,options);
 epochs = solution.x;
 orbit = solution.y;
 output_state = EcenToL2frame( orbit, epochs );
end

if ode_113
 options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
 solution = ode113(@full_force_model,[init_epoch final_time],Init_state,options);
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

