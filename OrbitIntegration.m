
clc
clear all
close all

%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

% Settings
full_mission = false; 
RKV_89 = true;
ABM = false;
RK45 = false;
PD78 = false;

if not(full_mission)
    load('irassihalotime.mat', 'Date');
    load('irassihalogmat.mat', 'Gmat');
    
    
else
    load('IRASSIFullMissionDate.mat', 'Date');
    load('IRASSIFullMission.mat', 'Gmat');
end

% Load kernel
cspice_furnsh ( METAKR );
planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN','301','5','VENUS','4','6'};
observer = 'EARTH';% or 339

% global G;
% G = 6.67e-20; % km % or -17


%% Ephemeris from SPICE
% Define initial epoch for a satellite
initial_utctime = '2030 MAY 22 00:03:25.693'; 
end_utctime = '2030 NOV 21 11:22:23.659';% NOV! %'2030 DEC 28 00:03:25.693'; %'2030 DEC 28 00:03:25.693';%'2030 NOV 21 11:22:23.659';
%'2030 DEC 28 00:03:25.693'; % 7 months
initial_et = cspice_str2et ( initial_utctime );
end_et = cspice_str2et ( end_utctime );
%step = 86400/10; %86400; %86400 3600 - every hour

if not(full_mission)
   et_vector = zeros(1,length(Date));
   for d=1:length(Date)
        utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
        et_vector(d) = cspice_str2et (utcdate);
   end
else
    et_vector = zeros(1,11621);
     for d=3245:1:14866-1
    utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
    et_vector(d-3244) = cspice_str2et (utcdate);
    end
end

disp(length(et_vector));


%% Setting up some values and structures
% Satellite initial position w.r.t the Earth center
initial_state = [-561844.307770134;-1023781.19884100;-152232.354717768;0.545714129191316;-0.288204299060291;-0.102116477725135]; 
% Create a structure for a satellite
sat = create_sat_structure(initial_state);
% Get initial states for calculating initial energy
[earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init] = create_structure( planets_name_for_struct, initial_et, observer);



%% Check influences
global influence;
influence = zeros(3,2);


%% INTEGRATION PART
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% ODE45
if RK45 == true
tic
orbit = ode45(@(t,y) force_model(t,y),et_vector,initial_state,options);    
toc
end

% Adams-Bashforth-Moulton Predictor-Corrector
if ABM == true
tic 
[orbit_ab8, tour] = adambashforth8(@force_model,et_vector,initial_state, length(et_vector));
toc
end

% Runge-Kutta-Verner 8(9)
simpleRKV89 = false;
embedded_estimation = true;
if RKV_89 == true
    
    if simpleRKV89 == true
       [orbit_rkv89, tourrkv] = RKV89(@force_model,et_vector,initial_state, length(et_vector));
    end
    if embedded_estimation == true
    tic
        orbit_rkv89_emb(:,1) = initial_state;
        next_step = 60; % initial value for next_step.
        final = false;
        n = 1;
        epochs(1) = et_vector(1);
        while not(final)
                [state, newstep, last] = rkv(@force_model,epochs(n),orbit_rkv89_emb(:,n), next_step, et_vector(length(et_vector)));
                next_step = newstep;
                final = last;

        n=n+1;
        epochs(n) = epochs(n-1) + next_step;
        orbit_rkv89_emb(:,n) = state;
        
        if n == 2 || n == 3
            disp(next_step);
        end
        
        end
        toc
    end
    % tic 
    % [orbit_rkv89, tourrkv] = rkv(@force_model,et_vector(1), et_vector(length(et_vector)),initial_state, 60);
    % toc
end

% Prince Dormand 7(8)
if PD78 == true
options87 = odeset('RelTol',1e-13,'AbsTol',1e-13, 'MaxStep',2700,'InitialStep',60);
[tour1, orbit_ode87] = ode45(@(t,y) force_model(t,y),et_vector,initial_state, options87);
orbit_ode87 = orbit_ode87';
end
toc


%% Maneuvers Insertion

% Integrate with maneuver inserted
%n_et = [9120, 14866]; % 2 for now
% n_et = 9120-3244; % 3244 - number of epoch before HALO orbit starts
% maneuver1 = [-0.02263165253058913;0.02267983525317713;-0.001364259283054504]; 
% %global t_at_etvector;
% %t_at_etvector = et_vector(n_et);
% for i=1:length(et_vector)
%     [orbit_rkv89, tourrkv] = RKV89(@force_model,et_vector,initial_state, length(et_vector));
%     if i == n_et
%         orbit_rkv89(4,i) = orbit_rkv89(4) + maneuver(1);
%         orbit_rkv89(5,i) = orbit_rkv89(5) + maneuver(2);
%         orbit_rkv89(6,i) = orbit_rkv89(6) + maneuver(3);
%     end
% end 

% for i=n_et+1:length(et_vector)
%     [orbit_rkv89, tourrkv] = RKV89(@force_model,et_vector,orbit_rkv89(6,n_et), length(et_vector));
%     if i == n_et
%         orbit_rkv89(4,i) = orbit_rkv89(4) + maneuver(1);
%         orbit_rkv89(5,i) = orbit_rkv89(5) + maneuver(2);
%         orbit_rkv89(6,i) = orbit_rkv89(6) + maneuver(3);
%     end
% end 




%% The differences
difference_rkv89emb = abs(Gmat(:,1:5875) - orbit_rkv89_emb(:,1:5875));
%difference_ab8 = abs(Gmat - orbit_ab8);
%difference_rkv89 = abs(Gmat - orbit_rkv89);


%% Total Energy checks

% energy = zeros(3, length(et_vector));  % 1 row Kinetic, 2 row Potential, 3 row - Total Mechanical
% 
% energy_ab4 = zeros(3, length(et_vector));
% First calculate the initial energies
% b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
% [init_total, init_kinetic, init_potential] = calculate_energy(b);
% Initial_energy = init_total;
% Initial_kinetic = init_kinetic;
% Initial_potential = init_potential;

%% Plotting

figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'b');% Reference
if RK45 == true
plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'r');% RK45
end
if ABM == true
plot3(orbit_ab8(1,:),orbit_ab8(2,:),orbit_ab8(3,:),'g'); % ABM8
end
if RKV_89 == true
    if simpleRKV89 == true
    plot3(orbit_rkv89(1,:),orbit_rkv89(2,:),orbit_rkv89(3,:),'c'); % RKV89
    end
plot3(orbit_rkv89_emb(1,:),orbit_rkv89_emb(2,:),orbit_rkv89_emb(3,:),'m'); % RKV89 with real error estimate
%plot3(orbit_rkv89(1,:),orbit_rkv89(2,:),orbit_rkv89(3,:),'c');
end
if PD78 == true
plot3(orbit_ode87(1,:),orbit_ode87(2,:),orbit_ode87(3,:),'y'); % RK87
end

figure(2)
grid on
hold on
plot(et_vector(1,1:5875),difference_rkv89emb(1,1:5875),et_vector(1,1:5875),difference_rkv89emb(2,1:5875),et_vector(1,1:5875),difference_rkv89emb(3,1:5875) );% Reference

figure(3)
grid on
hold on
plot(et_vector,difference_rkv89,et_vector,difference_rkv89,et_vector,difference_rkv89);% Reference

%% Plots info
figure(3)
title('Reference vs Integration');
legend('Reference','RK45','ABM8', 'RKV89', 'RKV89 embedded');
xlabel('x');
ylabel('y');
grid on


%cspice_kclear;