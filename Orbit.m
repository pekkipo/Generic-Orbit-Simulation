clc
clear all

%% Load kernels
METAKR = 'planetsorbitskernels.txt';
cspice_furnsh ( METAKR );
planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN',...
                           '301','5','VENUS','4','6'};
observer = 'EARTH';

%% Initial values
% Satellite initial position and velocity w.r.t the Earth center
initial_state =  [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005;...
                  5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];

initial_epoch = 958.910668311133e+006;
sat = create_sat_structure(initial_state);

%% Integration
orbit = [];
epochs = [];
final_point = 10000.747744831378550e+08; % The integrator should not reach it. 
complete = false;
init_t = initial_epoch;
init_state = initial_state;


% Shows the consecutive number of the maneuver applied
maneuver_number = 1;

% Number of required integrations. One integration - approximately 3 months
n_integrations = 6;
n = 1;

deltaVs = zeros(3,n_integrations);
        

global R0;
global V0;
global start_time;
     
 while ~complete
            
     % Automatic maneuver calculation
     % These global values are passed into the evalutaion.m function
     
     R0 = init_state(1:3);
     V0 = init_state(4:6);
     start_time = init_t;
     
     initial_guess = [0;0;0]; 
     
     % Calculate the maneuver
     deltaV = calculate_maneuver(initial_guess);
     % Add the maneuver to the array
     deltaVs = [deltaVs; deltaV];
     
     init_state(1:6) = init_state(1:6) + [0;0;0;deltaV(1);deltaV(2);deltaV(3)];

     [t, y0state, orbit_part, y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t , init_state);
     % t - epochs of current orbit part
     % y0state - state where y = 0 in L2 frame
     % orbit_part - positions and velocities of the sat
     % y0state_E - state where y=0 in L2 frame expressed in Earth-centered
     % coordinates
     
     orbit = [orbit, orbit_part];
     epochs = [epochs, t];
     
     % Change the values for the next orbit part integration
     init_t = epochs(end);
     init_state = y0state_E;
            
     maneuver_number = maneuver_number + 1;
            
     n = n+1;
     if n > n_integrations 
        complete = true; 
     end
 
 end
        
 
 % Write array with maneuvers into a data file
 save('automatic_maneuvers.mat','deltaVs');
 

%% PLOTTING
figure(1)
view(3)
grid on
hold on
plot3(0,0,0,'*r'); % nominal L2 point
plot3(orbit(1,:),orbit(2,:),orbit(3,:),'b'); % orbit
set(gca,'fontsize',16)

%% PLOTS INFO
figure(1)
title('HALO orbit around L2 SEM.');
legend('Nominal L2 point', 'Full model','Location','northeast');
set(gca,'fontsize',16)
xlabel('X, km');
ylabel('Y, km');
zlabel('Z, km');        
        
        
        
        