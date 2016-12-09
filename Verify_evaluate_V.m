%% Verify evaluate_V

% Initial state
R0 = [6.245032181291690e+05;   1.570803138086074e+06;  2.223197192596576e+05];  
V0 = [-2.387498156748843e-02;   1.696536087288626e-02;   7.399166018905870e-03];

% Initial Time
initial_time = 975.5469091830604e+06;

% Final Time
final_time = 993.140558185330e+006;

phi0 = reshape(eye(6), 36, 1);
init_state = [R0; V0; phi0];


% Ephemeris time vectors
METAKR = which('planetsorbitskernels.txt');
cspice_furnsh ( METAKR );

% EME_et = cspice_str2et(datestr(Edate2));

global G;
G = 6.673e-20;
global L2frame;
L2frame = true;
global checkrkv89_emb
checkrkv89_emb = false;

dV = [0; 0; -0.005];

deltaV = fsolve(@evaluate_V_test, dV);

Init_state = init_state;
Init_state(4:6,:) = Init_state(4:6,:)+ deltaV;
[t, y0state, output_state] = rkv89emb_maneuvers(@force_model_maneuvers, [initial_time final_time] , Init_state);

figure
hold on
plot3(output_state(1,:),output_state(2,:),output_state(3,:),'r','LineWidth',2)

