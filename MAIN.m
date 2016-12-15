% SCRIPT FOR THE MAIN RESULT
    % uses only one most effective integrator
    
    % 3 force models used

clc
clear all

METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

%% Settings

% Force model type
showsimple = false;
showsimplesrp = false;
showfull = false;
   

model = 'Simplified+SRP';

TEST = true;

%% Load kernel
cspice_furnsh ( METAKR );
planets_name_for_struct = {'EARTH','SUN','MOON','JUPITER','VENUS','MARS','SATURN';'EARTH','SUN',...
    '301','5','VENUS','4','6'};
observer = 'EARTH';

global G;
global L2frame;
global RKV_89_emb_check;
G = 6.673e-20;
L2frame = true;
RKV_89_emb_check = false;

%% Setting up some values and structures
% Satellite initial position w.r.t the Earth center
initial_state =  [-5.618445118318512e+005;  -1.023778587192635e+006;  -1.522315532439711e+005;...
    5.343825699573794e-001;  -2.686719669693540e-001;  -1.145921728828306e-001];

% initial_epoch = cspice_str2et ( 2030 MAY 22 00:03:25.693 ); USE IF ET EPOCH NOT KNOWN

initial_epoch = 958.910668311133e+006; % 22 May 2030
sat = create_sat_structure(initial_state);


%% Integration
if showsimplesrp
        simplesrp_orbit = [];
        simplesrp_epochs = [];
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
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        deltaVs = {dV1;dV2;dV3;dV4;dV5;dV6};
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 2;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            maneuver = deltaVs{maneuver_number};
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            phi0 = reshape(eye(6), 36, 1);
            init_state = [init_state(1:6); phi0];
            
            % RKV 89
            [epochs, y0state, orbit, y0state_E] = rkv89emb_maneuvers(@simplified_force_model_srp, init_t  , init_state);
            
            simplesrp_orbit = [simplesrp_orbit, orbit];
            simplesrp_epochs = [simplesrp_epochs, epochs];
            

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

if showfull
        full_orbit = [];
        full_epochs = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        % Maneuvers applied. Maneuvers for different integrators and models
        % are kept in the file MANEUVERS.TXT or can be calculated in
        % Calculate_Maneuvers.m script
        dV1 = [13.2530134946924e-003; -16.2338801932272e-003; 4.06482679813973e-003]; % CORRECT
        dV2 = [15.9090393682285e-003; 7.51017725084166e-003;-5.93871532821120e-003];
        %-5.93871532821120e-003]; % CORRECT less than origin
        %dV2 = [ 16.9057182084791e-003;10.7993500579884e-003;-3.55060805758670e-003]; % larger the origin
        dV3 = [26.5303022604513e-003;1.23533361008047e-003;4.00576011107553e-003];
        %dV3 = [26.5316539994299e-003; 1.23116499522077e-003; 4.00351131615303e-003];
       % dV3 = [25.6186952403759e-003; 954.286569193955e-006; 4.01796153186672e-003]; 
        %dV3 = [778.204109430765e-006;-6.23045235180619e-003;2.77553582883912e-003];
       % dV3 = [ 24.3633127693640e-003;3.42970257818512e-003;-2.06266375896980e-003];
        %dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
       % dV4 = [0;0;0];
        %dV4 = [137.011976423798e-003;70.3010595279474e-003;471.972514453013e-003];
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        
        
        %%%%% New try
        apply_simplified = false;
        if apply_simplified
        dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
        dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
            
            
        end

        
        
        
        %%%%%%

        
        deltaVs = {dV1;dV2;dV3;dV4;dV5;dV6};
        
        
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 2;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            maneuver = deltaVs{maneuver_number};
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            phi0 = reshape(eye(6), 36, 1);
            init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit, y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t  , init_state);

            full_orbit = [full_orbit, orbit];
            full_epochs = [full_epochs, epochs];
            
            
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

if showsimple
        simple_orbit = [];
        simple_epochs = [];
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
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        deltaVs = {dV1;dV2;dV3;dV4;dV5;dV6};
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 2;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            maneuver = deltaVs{maneuver_number};
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            phi0 = reshape(eye(6), 36, 1);
            init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit, y0state_E] = rkv89emb_maneuvers(@simplified_force_model, init_t , init_state);

            simple_orbit = [simple_orbit, orbit];
            simple_epochs = [simple_epochs, epochs];
            
            
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


%%
if TEST
        test_orbit = [];
        test_epochs = [];
        final_point = 10000.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;

        deltaVs = [];
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 6;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            %%%% Automatic maneuver calculation
            % these global values are passed into evalutaion function
            global R0;
            global V0;
            global start_time;
            
            R0 = init_state(1:3);
            V0 = init_state(4:6);
            start_time = init_t;
            
            deltaV = calculate_maneuver();
            %deltaVs(1:3, maneuver_number) = deltaV;
            deltaVs = [deltaVs; deltaV];
            %%%% 
            
            maneuver = deltaV;
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            phi0 = reshape(eye(6), 36, 1);
            init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit, y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t , init_state);

            test_orbit = [test_orbit, orbit];
            test_epochs = [test_epochs, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
        % Write cell array with maneuvers into a data file
       % dlmcell('rkv89emb_full_maneuvers_new.txt',deltaVs);
       save('rkv89emb_full_maneuvers_new.txt', 'deltaVs', '-ASCII');
        
        figure(3)
view(3)
grid on
hold on
plot3(0,0,0,'*r'); % nominal L2 point

    plot3(test_orbit(1,:),test_orbit(2,:),test_orbit(3,:),'b'); % orbit

end       

%% PLOTTING
figure(1)
view(3)
grid on
hold on
plot3(0,0,0,'*r'); % nominal L2 point
if showsimple
    plot3(simple_orbit(1,:),simple_orbit(2,:),simple_orbit(3,:),'b'); % orbit
end
if showsimplesrp
    plot3(simplesrp_orbit(1,:),simplesrp_orbit(2,:),simplesrp_orbit(3,:),'r'); % orbit
end
if showfull
    plot3(full_orbit(1,:),full_orbit(2,:),full_orbit(3,:),'g'); % orbit
end



%% Differences

sim_simsrp_x = abs(simple_orbit(1,:) - simplesrp_orbit(1,:);

%%

figure(2)
view(3)
grid on
hold on
plot3(0,0,0,'*r'); % nominal L2 point
if showsimple
    plot3(simple_orbit(1,:),simple_orbit(2,:),simple_orbit(3,:),'b'); % orbit
end
if showsimplesrp
    plot3(simplesrp_orbit(1,:),simplesrp_orbit(2,:),simplesrp_orbit(3,:),'r'); % orbit
end
if showfull
    plot3(full_orbit(1,:),full_orbit(2,:),full_orbit(3,:),'g'); % orbit
end



%% Plots info
figure(1)
title(['HALO orbit around L2 SEM. 3 Force Models']);
legend('Nominal L2 point','Simplified', 'Simplified + SRP', 'Full model');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

