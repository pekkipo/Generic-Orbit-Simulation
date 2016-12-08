
METAKR = 'planetsorbitskernels.txt';

cspice_furnsh ( METAKR );
planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN','301','5','VENUS','4','6'};
observer = 'EARTH';% or 339

global G;
G = 6.673e-20;
global L2frame;
L2frame = true;

initial_state = [5.795985038263178e+05; 7.776779586882917e+05;...
6.171179196351578e+05; -0.538364883921726; 0.286800406339146;0.125771126285189];
initial_et = 9.747626128418571e+08;

sat = create_sat_structure(initial_state);
% Get initial states for calculating initial energy
[earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init] = create_structure( planets_name_for_struct, initial_et, observer);

% Pos + velocites
%X0 = [5.795985038263178e+05; 7.776779586882917e+05;6.171179196351578e+05; -0.538364883921726; 0.286800406339146;0.125771126285189];

% Only velocities
 %V0 = [-0.58364883921726; 0.286800406339146;0.125771126285189];
 %V0 = [-5.399272545222726e-001;   2.861191946127703e-001;   1.254733378780861e-001];
% MY NEW VALUES
%V0 = [-0.538669578263083; 0.286257511925448; 0.125184841442128];
% v after three months
%V0 = [-0.017291522800783; 0.003697404618536; -0.007326275507767];
% Initial values before maneuver
 %V0 = [5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];
 %V0 = [5.457150405565953e-001;  -2.882041405454675e-001;  -1.021163220453542e-001];
 % after maneuver 5.457150405565953e-001  -2.882041405454675e-001  -1.021163220453542e-001
 
 % V0 after 3 months
 %V0 = [-0.017291522800783; 0.003697404618536; -0.007326275507767];

 % After 6 months with initial maneuver
 %V0 = [-0.537777527563825; 0.287834967935762; 0.126887706287377];
 % my vels after 6 months:
 %V0 = [-0.538669578247981; 0.286257511952463; 0.125184841471309];
 % with their maneuver
 %V0 = [-0.534973222258811; 0.281547765267124; 0.139797011011215];
 % my values [-0.58364883921726; 0.286800406339146;0.125771126285189];
% luisa [-5.399272545222726e-001;   2.861191946127703e-001;   1.254733378780861e-001];
% Add STM
%phi0 = [1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1];
%phi0 = eye(6);
%phi0 = reshape(phi, 36,1);
%X0 = [X0; phi0];

% Trying the orbit form the beginning
%V0 = [5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];

V0 =   [-26.1894854280253e-003;17.3026218026468e-003;9.61777015464231e-003];

%V0 = [0;0;0];
    
%options=optimoptions(@fsolve, 'Display', 'iter-detailed', 'Jacobian', 'on', 'TolFun', 1e-9);

%options=optimoptions(@fsolve, 'Algorithm', 'Levenberg-Marquardt','Display', 'iter-detailed', 'TolFun', 1e-3,'Jacobian', 'on','TolX', 1e-3);
%options=optimoptions(@fsolve, 'Display', 'iter-detailed', 'Jacobian', 'on');
options=optimoptions(@fsolve, 'Algorithm', 'Levenberg-Marquardt','Display', 'iter-detailed', 'TolFun', 1e-7,'Jacobian', 'on','TolX', 1e-6);

V = fsolve(@evaluate_V, V0, options);
%[gv, d_gv] = evaluate_V( V0 );

deltaV = V - V0;





