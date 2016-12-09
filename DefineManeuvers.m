
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

% Set initial state
% Insertion point
%1 R0 = [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005];  
% V0 = [5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];
 R0 = [1.485359168304580e+06;  -6.574308768747580e+05; -6.807000270951350e+05];  
 V0 = [-0.017761982577070;  -7.753583934581932e-04;  -0.001459319442820];

% Initial Time
initial_time = 9.669336468097093e+08;
% Final Time
final_time = 999.140558185330e+006;

phi0 = reshape(eye(6), 36, 1);
init_state = [R0; V0; phi0];


global G;
G = 6.673e-20;
global L2frame;
L2frame = true;
global checkrkv89_emb
checkrkv89_emb = false;

% Initial guess
dV = [0; 0; -0.005];


deltaV = fsolve(@evaluate_V_test, dV);
disp(deltaV);

Init_state = init_state;
Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;
[t, y0state, output_state, y0state_E] = rkv89emb_maneuvers(@force_model_maneuvers, [initial_time final_time] , Init_state);

% Init_state is the new value 
% y0state is the state, from which the next integration should start


% Graphical check of the orbit part
figure
hold on
plot3(output_state(1,:),output_state(2,:),output_state(3,:),'r','LineWidth',2)

