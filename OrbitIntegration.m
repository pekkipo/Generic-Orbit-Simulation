
clc
clear all
close all

%% Define local variables
METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

%% Settings
full_mission = false; % full mission or just a test part before the first maneuver
one_revolution = true; % only one maneuver applied % if false then all mission till the end
starting_from_earth = false; % mission with leop phase. Leave it false always!
RKV_89 = false;
ABM = true;
RK45 = true;
PD78 = true;
apply_maneuvers = false;
check_energy = false;

if not(full_mission)
    load('irassihalotime.mat', 'Date');
    load('irassihalogmat.mat', 'Gmat');
       
else
    load('IRASSIFullMissionDate.mat', 'Date');
    load('IRASSIFullMission.mat', 'Gmat');
end

%% Load kernel
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
ABMstep = 2700; %86400; %86400 3600 - every hour

abm_et_vector = initial_et:ABMstep:end_et;

if not(full_mission)
   et_vector = zeros(1,length(Date));
   for d=1:length(Date)
        utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
        et_vector(d) = cspice_str2et (utcdate);
   end
else
    if one_revolution == true
        et_vector = zeros(1,11621);
        for d=3245:1:14866-1
        utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
        et_vector(d-3244) = cspice_str2et (utcdate);
        end
    else
        if ~starting_from_earth
            et_vector = zeros(1,length(Date));
            for d=3245:length(Date)
            utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
            et_vector(d-3244) = cspice_str2et (utcdate);
            end
        else
            et_vector = zeros(1,length(Date));
            for d=1:length(Date)
            utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
            et_vector(d) = cspice_str2et (utcdate);
            end
        end
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

%% Maneuvers Insertion
% In HALO orbit phase
%              [ 21 Nov 2030 09:06:53.955 ; 19 May 2031 17:30:30.286; 19
%              Nov 2031 07:20:00.739; 16 May 2032 09:02:25.529; 
%              15 Nov 2032 20:12:11.617; 13 May 2033 12:06:40.408;
%              12 Nov 2033 20:30:57.484;10 May 2034 11:17:41.967
global epochs_numbers;
global maneuvers;

if starting_from_earth == true
    epochs_numbers = [9120; 14866; 20741; 26472; 32344; 38067; 43935; 49652]; % 8 maneuvers
else
    epochs_numbers = [9120; 14866; 20741; 26472; 32344; 38067; 43935; 49652];
    dif = 3244; % difference in epochs numbers. Halo phase starts at 3245. Which is 1 in my case
    for h = 1:length(epochs_numbers)
        epochs_numbers(h) = epochs_numbers(h) - dif;
    end
end   
if apply_maneuvers == true
    maneuver1 = [0.003696355989169846;-0.004709746685339394;0.01461216953990576]; 
    % my maneuver with corrections
    %maneuver1 = [0.003696355989169846-0.0016;-0.004709746685339394-0.0013;0.01461216953990576-0.0011];
    maneuver2 = [-0.004873280356119337;-0.007500302117829953;0.01748835216221812];          
    maneuver3 = [0.004395508963826083;-0.006574683170090312;0.01163890236306115];          
    maneuver4 = [-0.003729004790675886;-0.002912961186885862;0.01277290066887374];        
    maneuver5 = [0.003390651811125526;-0.005644399141577779;0.01143305771064858];           
    maneuver6 = [-0.003048218009434191;-0.003950944322589558;0.01099888048023886];        
    maneuver7 = [0.001693792364818979;-0.002281349967081013;0.008762453816533339];        
    maneuver8 = [-0.002822612198624273;-0.005352768508478975;0.01193703854590986];
    
    maneuvers = {maneuver1,maneuver2,maneuver3,maneuver4,maneuver5,maneuver6,maneuver7,maneuver8};
else
    zeromaneuvers = zeros(3,1);
    
    maneuvers = {zeromaneuvers,zeromaneuvers,zeromaneuvers,zeromaneuvers,zeromaneuvers,zeromaneuvers,zeromaneuvers,zeromaneuvers};
end




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
[orbit_ab8, tour] = adambashforth8(@force_model,abm_et_vector,initial_state, length(abm_et_vector));
toc
end

% Runge-Kutta-Verner 8(9)
simpleRKV89 = true;
embedded_estimation = false;
tic
if RKV_89 == true
    
    if simpleRKV89 == true
       %[orbit_rkv89, tourrkv] = RKV89(@force_model,et_vector,initial_state, length(et_vector));
       [orbit_rkv89, tourrkv] = RKV89_2(@force_model,et_vector,initial_state, length(et_vector));
    end
    if embedded_estimation == true
    
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
        orbit_rkv89_emb(:,n) = state;
        
        % Add maneuver if this is required epoch
        for k = 1:length(epochs_numbers)
            if n == epochs_numbers(k)
                % If this epoch is one of the epoch presented in maneuvers
                % array - add dV to its components
                applied_maneuver = maneuvers{k};
                orbit_rkv89_emb(4,n) = orbit_rkv89_emb(4,n) + applied_maneuver(1);
                orbit_rkv89_emb(5,n) = orbit_rkv89_emb(5,n) + applied_maneuver(2);
                orbit_rkv89_emb(6,n) = orbit_rkv89_emb(6,n) + applied_maneuver(3);
                
                next_step = 60; % Change next_step to 60 as if I started the integration from the beginning
                
            end
        end
        
        
        epochs(n) = epochs(n-1) + next_step;
        
        
            if n == 2 || n == 3
                disp(next_step);
            end
        end
        
    end
end
toc

% Prince Dormand 7(8)
if PD78 == true
    options87 = odeset('RelTol',1e-13,'AbsTol',1e-13, 'MaxStep',2700,'InitialStep',60);
    orbit_ode87 = zeros(6, length(et_vector));
    orbit_ode87(:,1) = initial_state;
    [pd78_et_vector, orbit_ode87] = ode87(@(t,y) force_model(t,y),[et_vector(1) et_vector(5879)], orbit_ode87(:,1), options87);
    orbit_ode87 = orbit_ode87';
    pd78_et_vector = pd78_et_vector';
    %     options87 = odeset('RelTol',1e-13,'AbsTol',1e-13, 'MaxStep',2700,'InitialStep',60);
%     orbit_ode87 = zeros(6, length(et_vector));
%     orbit_ode87(:,1) = initial_state;
%     for n = 1:length(et_vector)-1
%         [tour1, state] = ode87(@(t,y) force_model(t,y),[et_vector(n) et_vector(n+1)], orbit_ode87(:,n), options87);
%         state = state';
%         state = state(:,size(state,2));
%         orbit_ode87(:,n+1) = state;
%     end

end
%% Checking accuracy of integrators

%Reverse method
if RKV_89 == true
    et_vector_reversed = fliplr(et_vector);
    [orbit_rkv89_reversed, tourrkv] = RKV89_2(@force_model,et_vector_reversed,orbit_rkv89(:,length(orbit_rkv89)), length(et_vector_reversed));
    rkv89_conditions_difference = abs(fliplr(orbit_rkv89_reversed) - orbit_rkv89);
    rkv89_flp = fliplr(orbit_rkv89_reversed);
    rkv89_initial_value_difference = abs(rkv89_flp(:,1) - orbit_rkv89(:,1));
    disp('difference RKV_89');
    disp(rkv89_initial_value_difference);
end

if ABM == true
    abm_et_vector_reversed = fliplr(abm_et_vector);
    [orbit_ab8_reversed, tour] = adambashforth8(@force_model,abm_et_vector_reversed,orbit_ab8(:,length(orbit_ab8)), length(abm_et_vector_reversed));
    abm_conditions_difference = abs(fliplr(orbit_ab8_reversed) - orbit_ab8);
    abm_flp = fliplr(orbit_ab8_reversed);
    abm_initial_value_difference = abs(abm_flp(:,1) - orbit_ab8(:,1));
    disp('difference ABM');
    disp(abm_initial_value_difference);  
end

if RK45 == true
    et_vector_reversed = fliplr(et_vector);
    orbit_reversed = ode45(@(t,y) force_model(t,y),et_vector_reversed,orbit.y(:,length(orbit.y)),options);  
    rk_conditions_difference = abs(fliplr(orbit_reversed.y) - orbit.y(:,1:length(orbit_reversed.y)));
    rk_flp = fliplr(orbit_reversed.y);
    rk_initial_value_difference = abs(rk_flp(:,1) - orbit.y(:,1));
    disp('difference RK45');
    disp(rk_initial_value_difference);  
end

if PD78 == true
    pd78_vector_reversed = fliplr(pd78_et_vector);
    [tour1, orbit_ode87_reversed] = ode87(@(t,y) force_model(t,y),[pd78_vector_reversed(1) pd78_vector_reversed(length(pd78_vector_reversed))],orbit_ode87(:,length(orbit_ode87)), options87);
    orbit_ode87_reversed = orbit_ode87_reversed';
    pd78_conditions_difference = abs(fliplr(orbit_ode87_reversed) - orbit_ode87(:,1:length(orbit_ode87_reversed)));
    pd78_flp = fliplr(orbit_ode87_reversed);
    pd78_initial_value_difference = abs(pd78_flp(:,1) - orbit_ode87(:,1));
    disp('difference PD78');
    disp(pd78_initial_value_difference);  
end

%% Create a table with results

Integrators = {'RKV89';'ABM';'RK45';'PD78'};
if ~RKV_89 
    rkv89_initial_value_difference = zeros(6,1);
end
if ~ABM
    abm_initial_value_difference = zeros(6,1);
end
if ~RK45 
    rk_initial_value_difference = zeros(6,1);
end
if ~PD78 
    pd78_initial_value_difference = zeros(6,1);
end
Init_diffs = [rkv89_initial_value_difference;abm_initial_value_difference;rk_initial_value_difference;pd78_initial_value_difference];
%x_diffs = [rkv89_initial_value_difference;abm_initial_value_difference;rk_initial_value_difference;pd78_initial_value_difference ];
dX = [rkv89_initial_value_difference(1);abm_initial_value_difference(1);rk_initial_value_difference(1);pd78_initial_value_difference(1)];
dY = [rkv89_initial_value_difference(2);abm_initial_value_difference(2);rk_initial_value_difference(2);pd78_initial_value_difference(2)];
dZ = [rkv89_initial_value_difference(3);abm_initial_value_difference(3);rk_initial_value_difference(3);pd78_initial_value_difference(3)];
dVX = [rkv89_initial_value_difference(4);abm_initial_value_difference(4);rk_initial_value_difference(4);pd78_initial_value_difference(4)];
dVY = [rkv89_initial_value_difference(5);abm_initial_value_difference(5);rk_initial_value_difference(5);pd78_initial_value_difference(5)];
dVZ = [rkv89_initial_value_difference(6);abm_initial_value_difference(6);rk_initial_value_difference(6);pd78_initial_value_difference(6)];
dX_scalar = [sqrt(dX(1)^2+dY(1)^2+dZ(1)^2);sqrt(dX(2)^2+dY(2)^2+dZ(2)^2);sqrt(dX(3)^2+dY(3)^2+dZ(3)^2);sqrt(dX(4)^2+dY(4)^2+dZ(4)^2)];
dVX_scalar = [sqrt(dVX(1)^2+dVY(1)^2+dVZ(1)^2);sqrt(dVX(2)^2+dVY(2)^2+dVZ(2)^2);sqrt(dVX(3)^2+dVY(3)^2+dVZ(3)^2);sqrt(dVX(4)^2+dVY(4)^2+dVZ(4)^2)];
Table = table(dX,dY,dZ,dVX,dVY,dVZ,dX_scalar,dVX_scalar,'RowNames',Integrators);

%T = table(Integrators,Init_diffs);

%% The differences
%difference_rkv89emb = abs(Gmat(:,1:5859) - orbit_rkv89_emb(:,1:5859));
%difference_ab8 = abs(Gmat - orbit_ab8);
%difference_rkv89 = abs(Gmat - orbit_rkv89);


%% Total Energy checks
if check_energy == true
    if RKV_89 == true && simpleRKV89 == true
        energy_rkv89 = zeros(1, length(et_vector));
        % First calculate the initial energies
        b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
        [init_total, init_kinetic, init_potential] = calculate_energy(b);
        Initial_energy = init_total;
        Initial_kinetic = init_kinetic;
        Initial_potential = init_potential;

        % Calculate for each step
        for epoch = 1:length(et_vector)
            % Create a structure for the satellite
            sat_at_this_time = create_sat_structure(orbit_rkv89(:,epoch));
            % Information about planets at a given epoch
            [earth, sun, moon, jupiter, venus, mars, saturn] = create_structure( planets_name_for_struct, et_vector(epoch), observer);
            bodies = [sat_at_this_time, earth, sun, moon, jupiter, venus, mars, saturn];
            [total, kinetic, potential] = calculate_energy(bodies);
            kin1 = kinetic - Initial_kinetic;
            pot1 = potential - Initial_potential;
            tot1 = total - Initial_energy;

            energy_rkv89(1,epoch) = abs(tot1);
        end
    end
    
    if ABM == true
        energy_abm = zeros(1, length(abm_et_vector));
        % First calculate the initial energies
        b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
        [init_total, init_kinetic, init_potential] = calculate_energy(b);
        Initial_energy = init_total;
        Initial_kinetic = init_kinetic;
        Initial_potential = init_potential;

        % Calculate for each step
        for epoch = 1:length(abm_et_vector)
            % Create a structure for the satellite
            sat_at_this_time = create_sat_structure(orbit_ab8(:,epoch));
            % Information about planets at a given epoch
            [earth, sun, moon, jupiter, venus, mars, saturn] = create_structure( planets_name_for_struct, et_vector(epoch), observer);
            bodies = [sat_at_this_time, earth, sun, moon, jupiter, venus, mars, saturn];
            [total, kinetic, potential] = calculate_energy(bodies);
            kin1 = kinetic - Initial_kinetic;
            pot1 = potential - Initial_potential;
            tot1 = total - Initial_energy;

            energy_abm(1,epoch) = abs(tot1);
        end
    end
    
    if RK45 == true
        energy_rk = zeros(1, length(et_vector));
        % First calculate the initial energies
        b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
        [init_total, init_kinetic, init_potential] = calculate_energy(b);
        Initial_energy = init_total;
        Initial_kinetic = init_kinetic;
        Initial_potential = init_potential;

        % Calculate for each step
        for epoch = 1:length(orbit.y)
            % Create a structure for the satellite
            sat_at_this_time = create_sat_structure(orbit.y(:,epoch));
            % Information about planets at a given epoch
            [earth, sun, moon, jupiter, venus, mars, saturn] = create_structure( planets_name_for_struct, et_vector(epoch), observer);
            bodies = [sat_at_this_time, earth, sun, moon, jupiter, venus, mars, saturn];
            [total, kinetic, potential] = calculate_energy(bodies);
            kin1 = kinetic - Initial_kinetic;
            pot1 = potential - Initial_potential;
            tot1 = total - Initial_energy;

            energy_rk(1,epoch) = abs(tot1);
        end
    end
    
    if PD78 == true
        energy_pd78 = zeros(1, length(et_vector));
        % First calculate the initial energies
        b = [sat, earth_init, sun_init, moon_init, jupiter_init, venus_init, mars_init, saturn_init];
        [init_total, init_kinetic, init_potential] = calculate_energy(b);
        Initial_energy = init_total;
        Initial_kinetic = init_kinetic;
        Initial_potential = init_potential;

        % Calculate for each step
        for epoch = 1:length(et_vector)
            % Create a structure for the satellite
            sat_at_this_time = create_sat_structure(orbit_ode87(:,epoch));
            % Information about planets at a given epoch
            [earth, sun, moon, jupiter, venus, mars, saturn] = create_structure( planets_name_for_struct, et_vector(epoch), observer);
            bodies = [sat_at_this_time, earth, sun, moon, jupiter, venus, mars, saturn];
            [total, kinetic, potential] = calculate_energy(bodies);
            kin1 = kinetic - Initial_kinetic;
            pot1 = potential - Initial_potential;
            tot1 = total - Initial_energy;

            energy_pd78(1,epoch) = abs(tot1);
        end
    end
end
%% Plotting
figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'b');% Reference
%plot3(Gmat(1,1:15000),Gmat(2,1:15000),Gmat(3,1:15000),'b');
if RK45 == true
    plot3(orbit.y(1,:),orbit.y(2,:),orbit.y(3,:),'r');% RK45
end
if ABM == true
    plot3(orbit_ab8(1,:),orbit_ab8(2,:),orbit_ab8(3,:),'g'); % ABM8
    difference_abm = abs(Gmat(:,1:length(orbit_ab8)) - orbit_ab8);
end
if RKV_89 == true
    if simpleRKV89 == true
    plot3(orbit_rkv89(1,:),orbit_rkv89(2,:),orbit_rkv89(3,:),'c'); % RKV89
    difference_rkv89 = abs(Gmat(:,1:length(orbit_rkv89)) - orbit_rkv89);
    end
    if embedded_estimation == true
    plot3(orbit_rkv89_emb(1,:),orbit_rkv89_emb(2,:),orbit_rkv89_emb(3,:),'m'); % RKV89 with real error estimate
    end
    %plot3(orbit_rkv89(1,:),orbit_rkv89(2,:),orbit_rkv89(3,:),'c');
end
if PD78 == true
    plot3(orbit_ode87(1,:),orbit_ode87(2,:),orbit_ode87(3,:),'y'); % RK87
    difference_pd78 = abs(Gmat(:,1:length(orbit_ode87)) - orbit_ode87);
end
% figure(2)
% grid on
% hold on
% plot(et_vector(1,1:5859),difference_rkv89emb(1,1:5859),et_vector(1,1:5859),difference_rkv89emb(2,1:5859),et_vector(1,1:5859),difference_rkv89emb(3,1:5859) );% Reference
if RKV_89 == true
    figure(3)
    grid on
    hold on
    plot(et_vector,difference_rkv89(1,:),et_vector,difference_rkv89(2,:),et_vector,difference_rkv89(3,:));% Reference

end

if PD78 == true
    figure(4)
    grid on
    hold on
    plot(et_vector,difference_pd78(1),et_vector,difference_pd78(2),et_vector,difference_pd78(3));% Reference
end

if ABM == true
    figure(5)
    grid on
    hold on
    plot(abm_et_vector,difference_abm(1,:),abm_et_vector,difference_abm(2,:),abm_et_vector,difference_abm(3,:));% Reference

end


%% Energy calculation figures
if check_energy == true
    if RKV_89 == true
        figure(8)
        grid on
        hold on
        plot(et_vector,energy_rkv89);
    end

    if ABM == true
        figure(9)
        grid on
        hold on
        plot(abm_et_vector,energy_abm);
    end

    if RK45 == true
        figure(10)
        grid on
        hold on
        plot(et_vector,energy_rk);
    end
    
    if PD78 == true
        figure(11)
        grid on
        hold on
        plot(et_vector,energy_pd78);
    end
end



%% Plots info
figure(3)
title('Reference vs Integration');
legend('Reference','RK45','ABM8', 'RKV89', 'RKV89 embedded');
xlabel('x');
ylabel('y');
grid on


