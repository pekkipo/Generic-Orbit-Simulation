clc
clear all

METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

%% Settings
RKV_89 = false;
RKV_89_emb = true;
ABM8 = false;
ODE113 = false;
ODE45 = false;
ODE87 = false;
COWELL = false;

reverse_check = false;

global L2frame;
L2frame = true;


global RKV_89_check;
global RKV_89_emb_check;
global ABM8_check;
global ODE113_check;
global ODE45_check;
global ODE87_check;
global COWELL_check;


% initially set to false. Changes automatically to true below in the next if statement

RKV_89_check = false;
RKV_89_emb_check = false;
ABM8_check = false;
ODE113_check = false;
ODE45_check = false;
ODE87_check = false;
COWELL_check = false;

% Force model type
   % this is changed within the force_model function
   % change it here manually for displaying on the graphs
model = 'Simplified+SRP';

%% Load kernel
cspice_furnsh ( METAKR );
planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN',...
    '301','5','VENUS','4','6'};
observer = 'EARTH';

global G;
G = 6.673e-20;

%% Setting up some values and structures
% Satellite initial position w.r.t the Earth center
initial_state =  [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005;...
    5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];

% initial_epoch = cspice_str2et ( 2030 MAY 22 00:03:25.693 ); USE IF ET EPOCH NOT KNOWN

initial_epoch = 958.910668311133e+006; % 22 May 2030
sat = create_sat_structure(initial_state);


%% Integration
if RKV_89_emb
    
        orbit_RKV_89 = [];
        totalepochs_rkv89 = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        % Maneuvers applied. Maneuvers for different integrators and models
        % are kept in the file MANEUVERS.TXT or can be calculated in
        % Calculate_Maneuvers.m script
        dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
        dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [];
        dV6 = [];
        deltaVs = {dV1;dV2;dV3;dV4;dV5;dV6};
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 4;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            maneuver = deltaVs{maneuver_number};
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            phi0 = reshape(eye(6), 36, 1);
            init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit_rkv89_emb, y0state_E] = rkv89emb_maneuvers(@force_model, [init_t final_point] , init_state);

            orbit_rkv_89 = [totalorbit_rkv89, orbit_rkv89_emb];
            totalepochs_rkv89 = [totalepochs_rkv89, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
end


