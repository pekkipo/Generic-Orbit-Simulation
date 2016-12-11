
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

% Set initial state

 R0 = [1.441868581740170e+06; -7.148156992129891e+05; -6.984326839322323e+05];
 V0 = [-0.001817314659897;-0.007755542122498;-0.002671479865356];
% Initial Time
initial_time = 9.982927336197942e+08;
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

% Initial guess
dV = [-0.001625936670348; -0.003125208256016; -0.008088501084076];


deltaV = fsolve(@evaluate_V_test, dV);
disp(deltaV);

Init_state = init_state;
Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;
[t, y0state, output_state, y0state_E] = rkv89emb_maneuvers(@simplified_force_model_srp, initial_time , Init_state);

% Init_state is the new value 
% y0state is the state, from which the next integration should start


% Graphical check of the orbit part
figure
hold on
plot3(output_state(1,:),output_state(2,:),output_state(3,:),'r','LineWidth',2)

