clc
clear all

METAKR = 'planetsorbitskernels.txt';%'satelliteorbitkernels.txt';

%% Settings
RKV_89 = false;
RKV_89_emb = true;
ABM8 = false;
ODE113 = true;
ODE45 = false;
ODE87 = false;

global reverse_check;
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

global ode87_lastpiece;
ode87_lastpiece = false;

% Force model type
   % this is changed within the force_model function
   % change it here manually for displaying on the graphs
model = 'Full force model';

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

% RKV89 EMBEDDED
if RKV_89_emb
    
        orbit_RKV_89_emb = [];
        totalepochs_rkv89_emb = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        
        dV1 = [0.0132567757055320;-0.0162165135037194;0.00404055680235709]; 
        dV2 = [-0.000780848477937322;0.00185419489886210;-0.00724680203763881];
        dV3 = [0.00254415412213867;-0.00292144414696504;0.00770244866114741];
        dV4 = [-0.00162687225052243;-0.00312491482142916;-0.00808871229278551];
        dV5 = [-0.00292347011146131;-0.00338614109574762;0.00835218288400760];
        dV6 = [-0.00234439794209358;-0.00239804821203765;-0.00874185935583752];

       % dV1 = [13.2530134946924e-003; -16.2338801932272e-003; 4.06482679813973e-003]; % CORRECT
        %dV2 = [15.9090393682285e-003; 7.51017725084166e-003;-5.93871532821120e-003];
        
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
            %phi0 = reshape(eye(6), 36, 1);
            %init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit_rkv89_emb, rkv89emb_y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t , init_state);
            
            
            
            orbit_RKV_89_emb = [orbit_RKV_89_emb, orbit_rkv89_emb];
            totalepochs_rkv89_emb = [totalepochs_rkv89_emb, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rkv89emb_y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                 
               complete = true; 
            end
            
        end
        
       % point_to_compare1 = [-381206.986695521;-4.39933501183987e-07;267234.597932537;...
            %0.00388945312992179;0.397958033842957;0.0107094127665944];
         % point_to_compare1_rkv89emb =  [-381206.986695522;-3.10683390125632e-08;267234.597932545;...
             %   0.00388945312987548;0.397958033842940;0.0107094127665364];
             point_to_compare1_rkv89emb = [-381207.017827721;3.86062311008573e-08;267234.468573306;...
                    0.00388273530349280;0.397957765020769;0.0106802647964312];
            
        
end

% RKV89

if RKV_89
    
        orbit_RKV_89 = [];
        totalepochs_rkv89 = [];
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        dV1 = [0.0132567757055320;-0.0162165135037194;0.00404055680235709]; 
        dV2 = [-0.000780848477937322;0.00185419489886210;-0.00724680203763881];
        dV3 = [0.00254415412213867;-0.00292144414696504;0.00770244866114741];
        dV4 = [-0.00162687225052243;-0.00312491482142916;-0.00808871229278551];
        dV5 = [-0.00292347011146131;-0.00338614109574762;0.00835218288400760];
        dV6 = [-0.00234439794209358;-0.00239804821203765;-0.00874185935583752];

        %dV1 = [13.2530134946924e-003; -16.2338801932272e-003; 4.06482679813973e-003]; % CORRECT
        %dV2 = [15.9090393682285e-003; 7.51017725084166e-003;-5.93871532821120e-003];

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
            %phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit_rkv89, rkv89_y0state_E] = full_RKV89(@full_force_model, init_t, init_state);

            orbit_RKV_89 = [orbit_RKV_89, orbit_rkv89];
            totalepochs_rkv89 = [totalepochs_rkv89, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rkv89_y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
        %point_to_compare1_rkv89 =  [-381206.986695513;-9.33068804442883e-08;267234.597932549;...
         %   0.00388945312996990;0.397958033842970;0.0107094127665074];
        point_to_compare1_rkv89 = [-381207.017827715;-1.22963683679700e-08;267234.468573308;...
            0.00388273530357262;0.397957765020794;0.0106802647964059];
        
end

if ABM8
    
        orbit_ABM8 = [];
        totalepochs_abm8 = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        % Maneuvers applied. Maneuvers for different integrators and models
        % are kept in the file MANEUVERS.TXT or can be calculated in
        % Calculate_Maneuvers.m script
%         dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
%         dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [];
        dV6 = [];
        
        dV1 = [13.2530134946924e-003; -16.2338801932272e-003; 4.06482679813973e-003]; % CORRECT
        dV2 = [15.9090393682285e-003; 7.51017725084166e-003;-5.93871532821120e-003];

        
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
            %phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit_abm8, y0state_E] = ABM8(@full_force_model, init_t , init_state);

            orbit_ABM8 = [orbit_ABM8, orbit_abm8];
            totalepochs_abm8 = [totalepochs_abm8, epochs];
            
            
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

if ODE113
    
        orbit_ODE113 = [];
        totalepochs_ode113 = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        % Maneuvers applied. Maneuvers for different integrators and models
        % are kept in the file MANEUVERS.TXT or can be calculated in
%         % Calculate_Maneuvers.m script
%         dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
%         dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [];
        dV6 = [];
        
        dV1 = [0.0132624893374001;-0.0162513152365047;0.00404609981491367];
        dV2 = [-0.000780457538397596;0.00185525021798016;-0.00724844974057293];
        
        %deltaVs = {[0;0;0];dV1;dV2;dV3;dV4;dV5;dV6};
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
            %phi0 = reshape(eye(6), 36, 1);
            %init_state = [init_state(1:6); phi0];
            
            
%             options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
%             solution = ode113(@full_force_model,[init_t final_point],init_state,options);
%             epochs = solution.x;
%             orbit_ode113 = solution.y;
%             ode113_y0state_E = orbit_ode113(:,end); %solution.ye;
%             % Event handler does inner transformation to check y=0
%             % Now I need to transform the result to L2 frame
%             orbit_ode113 = EcenToL2frame( orbit_ode113, epochs );
            [epochs, y0state, orbit_ode113, ode113_y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t , init_state);


            orbit_ODE113 = [orbit_ODE113, orbit_ode113];
            totalepochs_ode113 = [totalepochs_ode113, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = ode113_y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
end

if ODE45
    
        orbit_ODE45 = [];
        totalepochs_ode45 = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        

        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [];
        dV6 = [];
        
        dV1 = [0.0132624893374001;-0.0162513152365047;0.00404609981491367];
        dV2 = [-0.000780457537688229;0.00185525021761572;-0.00724844974038747];
        
        deltaVs = {[0;0;0];dV1;dV2;dV3;dV4;dV5;dV6};
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        % Number of required integrations. One integration - approximately
        % 3 months or when y = 0. Half of the revolution
        n_integrations = 3;
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            maneuver = deltaVs{maneuver_number};
            init_state(1:6) = init_state(1:6) + [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
           % phi0 = reshape(eye(6), 36, 1);
            %init_state = [init_state(1:6); phi0];
            
            
             
            options = odeset('Events',@event_handler, 'MaxStep', 2700, 'InitialStep', 60);
            solution = ode45(@full_force_model,[init_t final_point],init_state,options);
            epochs = solution.x;
            orbit_ode45 = solution.y;
            ode45_y0state_E = orbit_ode45(:,end); %solution.ye;
            % Event handler does inner transformation to check y=0
            % Now I need to transform the result to L2 frame
            orbit_ode45 = EcenToL2frame( orbit_ode45, epochs );

            orbit_ODE45 = [orbit_ODE45, orbit_ode45];
            totalepochs_ode45 = [totalepochs_ode45, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = ode45_y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
end

if ODE87
    
        orbit_ODE87 = [];
        totalepochs_ode87 = [];
        final_point = 100.747744831378550e+08; % The integrator should not reach it.
        
        complete = false;
        init_t = initial_epoch;
        init_state = initial_state;
        
        % Maneuvers applied. Maneuvers for different integrators and models
        % are kept in the file MANEUVERS.TXT or can be calculated in
        % Calculate_Maneuvers.m script
%         dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
%         dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [];
        dV6 = [];
        
       dV1 = [0.0132567757114457;-0.0162165134933847;0.00404055677593769];
       dV2 =[-780.848492290756e-006;1.85419489025320e-003;-7.24680207511490e-003];
        
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
           % phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            
            [epochs, y0state, orbit_ode87, ode87_y0state_E] = full_ode87(@full_force_model, [init_t final_point] , init_state);

            orbit_ODE87 = [orbit_ODE87, orbit_ode87];
            totalepochs_ode87 = [totalepochs_ode87, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = ode87_y0state_E;
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
               complete = true; 
            end
            
        end
        
        
       % point_to_compare1_ode87 = [-381206.986695515;-7.76926754042506e-08;267234.597932548;...
            %0.00388945312995299;0.397958033842946;0.0107094127666117];
            point_to_compare1_ode87 = [-381207.017827733;-5.43659552931786e-08;267234.468573176;...
                0.00388273530147918;0.397957765026424;0.0106802647681097];
end


%% Checking accuracy of the integrators

if reverse_check == true
    %Reverse method
    if RKV_89_emb == true

        RKV_89_emb_check = true;
            
        reverse_orbit_RKV_89_emb = [];
        reverse_totalepochs_rkv89_emb = [];
        final_point = 938.910668311133e+006; % The integrator should not reach it.
        
        complete = false;
        init_t = totalepochs_rkv89_emb(end);
        init_state = rkv89emb_y0state_E;%orbit_RKV_89(1:6,end);
        % rkv89_y0state is the last point of the orbit but in EME!

        dV1 = [0.0132567757055320;-0.0162165135037194;0.00404055680235709]; 
        dV2 = [-0.000780848477937322;0.00185419489886210;-0.00724680203763881];
        dV3 = [0.00254415412213867;-0.00292144414696504;0.00770244866114741];
        dV4 = [-0.00162687225052243;-0.00312491482142916;-0.00808871229278551];
        dV5 = [-0.00292347011146131;-0.00338614109574762;0.00835218288400760];
        dV6 = [-0.00234439794209358;-0.00239804821203765;-0.00874185935583752];
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;

        n_integrations = 2;
        % 3 integrations means that one should start from maneuver 3
        deltaVs = {dV2;dV1;[0;0;0]};
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete

            %phi0 = reshape(eye(6), 36, 1);
            %init_state = [init_state(1:6); phi0];
             
            [epochs, y0state, rev_orbit_rkv89_emb, rev_rkv89emb_y0state_E] = full_rkv89emb_maneuvers(@full_force_model, init_t , init_state);

            maneuver = deltaVs{maneuver_number};
            rev_rkv89emb_y0state_E(1:6) = rev_rkv89emb_y0state_E(1:6) - [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
           
            reverse_orbit_RKV_89_emb = [reverse_orbit_RKV_89_emb, rev_orbit_rkv89_emb];
            reverse_totalepochs_rkv89_emb = [reverse_totalepochs_rkv89_emb, epochs];
            
            
            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rev_rkv89emb_y0state_E;      
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                point_to_compare2_rkv89emb = y0state;
                % Before this run a small piece from last y = 0 till state at t = t_initial
               start_time = 958.910668311133e+006;
               global rkv89emb_lastpiece;
               rkv89emb_lastpiece = true;
               [epochs, y0state, rev_orbit_rkv89emb, rev_rkv89emb_y0state_E] = full_rkv89emb_maneuvers(@full_force_model, [init_t start_time], init_state);
                
                reverse_orbit_RKV_89_emb = [reverse_orbit_RKV_89_emb, rev_orbit_rkv89emb];
                reverse_totalepochs_rkv89_emb = [reverse_totalepochs_rkv89_emb, epochs];
                
               complete = true; 
            end
            
        end
        rkv89emb_lastpiece = false; % set back to false

            %rkv89emb_conditions_difference = abs(fliplr(reverse_orbit_RKV_89_emb) - orbit_RKV_89_emb);
            rkv89emb_flp = fliplr(reverse_orbit_RKV_89_emb);
            %rkv89emb_initial_value_difference = abs(rkv89emb_flp(:,1) - orbit_RKV_89_emb(:,1));
            rkv89emb_initial_value_difference = abs(point_to_compare2_rkv89emb(1:6) - point_to_compare1_rkv89emb);
            disp('difference RKV_89_emb');
            disp(rkv89emb_initial_value_difference);

    end

     if RKV_89 == true

        RKV_89_check = true;
            
        reverse_orbit_RKV_89 = [];
        reverse_totalepochs_rkv89 = [];
        final_point = 938.910668311133e+006; % The integrator should not reach it.
        
        complete = false;
        % Init have to be in EME!
        init_t = totalepochs_rkv89(end);
        init_state = rkv89_y0state_E;%orbit_RKV_89(1:6,end);
        % rkv89_y0state is the last point of the orbit but in EME!

        dV1 = [0.0132567757055320;-0.0162165135037194;0.00404055680235709]; 
        dV2 = [-0.000780848477937322;0.00185419489886210;-0.00724680203763881];
        dV3 = [0.00254415412213867;-0.00292144414696504;0.00770244866114741];
        dV4 = [-0.00162687225052243;-0.00312491482142916;-0.00808871229278551];
        dV5 = [-0.00292347011146131;-0.00338614109574762;0.00835218288400760];
        dV6 = [-0.00234439794209358;-0.00239804821203765;-0.00874185935583752];
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        n_integrations = 2;
        % Adjust maneuvers manually
        % 3 integrations means that one should start from maneuver 3
        deltaVs = {dV2;dV1};
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
            %phi0 = reshape(eye(6), 36, 1);
            %init_state = [init_state(1:6); phi0];
            
            [epochs, y0state, rev_orbit_rkv89, rev_rkv89_y0state_E] = full_RKV89(@full_force_model, init_t , init_state);

            % Now subtract the maneuver
            
            maneuver = deltaVs{maneuver_number};
            rev_rkv89_y0state_E(1:6) = rev_rkv89_y0state_E(1:6) - [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            
            reverse_orbit_RKV_89 = [reverse_orbit_RKV_89, rev_orbit_rkv89];
            reverse_totalepochs_rkv89 = [reverse_totalepochs_rkv89, epochs];
            

            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rev_rkv89_y0state_E;     
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                point_to_compare2_rkv89 = y0state;
               % Before this run a small piece from last y = 0 till state at t = t_initial
               start_time = 958.910668311133e+006;
               global rkv89_lastpiece;
               rkv89_lastpiece = true;
               [epochs, y0state, rev_orbit_rkv89, rev_rkv89_y0state_E] = full_RKV89(@full_force_model, [init_t start_time], init_state);
                
                reverse_orbit_RKV_89 = [reverse_orbit_RKV_89, rev_orbit_rkv89];
                reverse_totalepochs_rkv89 = [reverse_totalepochs_rkv89, epochs];
                
               complete = true; 
            end
            
        end
            rkv89_lastpiece = false; % set back to false
            %rkv89_conditions_difference = abs(fliplr(reverse_orbit_RKV_89) - orbit_RKV_89);
            rkv89_flp = fliplr(reverse_orbit_RKV_89);
           % rkv89_initial_value_difference = abs(rkv89_flp(1:6,1) - orbit_RKV_89(1:6,1));
           rkv89_initial_value_difference = abs(point_to_compare2_rkv89(1:6) - point_to_compare1_rkv89);
            disp('difference RKV_89');
            disp(rkv89_initial_value_difference);

     end

    
     if ODE87 == true
         
        ODE87_check = true;
            
        reverse_orbit_ODE87 = [];
        reverse_totalepochs_ode87 = [];
        final_point = 938.910668311133e+006; % The integrator should not reach it.
        
        complete = false;
        % Init have to be in EME!
        init_t = totalepochs_ode87(end);
        init_state = ode87_y0state_E;%orbit_RKV_89(1:6,end);
        % rkv89_y0state is the last point of the orbit but in EME!
%         
%         dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
%         dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        
        dV1 = [0.0132567757114457;-0.0162165134933847;0.00404055677593769];
        dV2 =[-780.848492290756e-006;1.85419489025320e-003;-7.24680207511490e-003];
       
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        n_integrations = 2;
        % Adjust maneuvers manually
        % 3 integrations means that one should start from maneuver 3
        deltaVs = {dV2;dV1};
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
           % phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            [epochs, y0state, rev_orbit_ode87, rev_ode87_y0state_E] = full_ode87(@full_force_model, [init_t final_point], init_state);

            % Now subtract the maneuver
            
            maneuver = deltaVs{maneuver_number};
            rev_ode87_y0state_E(1:6) = rev_ode87_y0state_E(1:6) - [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            
            reverse_orbit_ODE87 = [reverse_orbit_ODE87, rev_orbit_ode87];
            reverse_totalepochs_ode87 = [reverse_totalepochs_ode87, epochs];
            

            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rev_ode87_y0state_E;     
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                point_to_compare2_ode87 = y0state;
                 % Before this run a small piece from last y = 0 till state at t = t_initial
               start_time = 958.910668311133e+006;
             

               [epochs, y0state, rev_orbit_ode87, rev_ode87_y0state_E] = full_ode87(@full_force_model, [init_t start_time], init_state);
                
                reverse_orbit_ODE87  = [reverse_orbit_ODE87, rev_orbit_ode87];
                reverse_totalepochs_ode87 = [reverse_totalepochs_ode87, epochs];
                
               complete = true; 
            end
            
            
            
        end
        
      

            %rkv89_conditions_difference = abs(fliplr(reverse_orbit_RKV_89) - orbit_RKV_89);
           	ode87_flp = fliplr(reverse_orbit_ODE87);
            %ode87_initial_value_difference = abs(ode87_flp(1:6,1) - orbit_ODE87(1:6,1));
            ode87_initial_value_difference = abs(point_to_compare2_ode87(1:6) - point_to_compare1_ode87);
            disp('difference ODE87');
            disp(ode87_initial_value_difference);

         
     end
     
      if ODE45 == true
         
        ODE45_check = true;
            
        reverse_orbit_ODE45 = [];
        reverse_totalepochs_ode45 = [];
        final_point = 938.910668311133e+006; % The integrator should not reach it.
        
        complete = false;
        % Init have to be in EME!
        init_t = totalepochs_ode45(end);
        init_state = ode45_y0state_E;%orbit_RKV_89(1:6,end);
        % rkv89_y0state is the last point of the orbit but in EME!

        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        
        dV1 = [0.0132624893374001;-0.0162513152365047;0.00404609981491367];
        dV2 = [-0.000780457537688229;0.00185525021761572;-0.00724844974038747];
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        n_integrations = 2;
        % Adjust maneuvers manually
        % 3 integrations means that one should start from maneuver 3
        deltaVs = {dV2;dV1;[0;0;0]};
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
          %  phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            options = odeset('Events',@event_handler, 'MaxStep', -2700, 'InitialStep', -60);
            solution = ode45(@full_force_model,[init_t final_point],init_state,options);
            epochs = solution.x;
            rev_orbit_ode45 = solution.y;
            rev_ode45_y0state_E = rev_orbit_ode45(:,end); %solution.ye;
            % Event handler does inner transformation to check y=0
            % Now I need to transform the result to L2 frame
            rev_orbit_ode45 = EcenToL2frame( rev_orbit_ode45, epochs );
            % Now subtract the maneuver
            
            maneuver = deltaVs{maneuver_number};
            rev_ode45_y0state_E(1:6) = rev_ode45_y0state_E(1:6) - [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            
            reverse_orbit_ODE45 = [reverse_orbit_ODE45, rev_orbit_ode45];
            reverse_totalepochs_ode45 = [reverse_totalepochs_ode45, epochs];
            

            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rev_ode45_y0state_E;     
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                
                 % Before this run a small piece from last y = 0 till state at t = t_initial
               start_time = 958.910668311133e+006;

            options = odeset('MaxStep', -2700, 'InitialStep', -60);
            solution = ode45(@full_force_model,[init_t start_time],init_state,options);
            epochs = solution.x;
            rev_orbit_ode45 = solution.y;
            rev_ode45_y0state_E = rev_orbit_ode45(:,end);
            rev_orbit_ode45 = EcenToL2frame( rev_orbit_ode45, epochs );
              
                reverse_orbit_ODE45  = [reverse_orbit_ODE45, rev_orbit_ode45];
                reverse_totalepochs_ode45 = [reverse_totalepochs_ode45, epochs];
                
               complete = true; 
            end
            
        end

            %rkv89_conditions_difference = abs(fliplr(reverse_orbit_RKV_89) - orbit_RKV_89);
           	ode45_flp = fliplr(reverse_orbit_ODE45);
            ode45_initial_value_difference = abs(ode45_flp(1:6,1) - orbit_ODE45(1:6,1));
            disp('difference ODE45');
            disp(ode45_initial_value_difference);

         
      end
     
      if ODE113 == true
         
        ODE113_check = true;
            
        reverse_orbit_ODE113 = [];
        reverse_totalepochs_ode113 = [];
        final_point = 938.910668311133e+006; % The integrator should not reach it.
        
        complete = false;
        % Init have to be in EME!
        init_t = totalepochs_ode113(end);
        init_state = ode113_y0state_E;%orbit_RKV_89(1:6,end);
        % rkv89_y0state is the last point of the orbit but in EME!
        
%         dV1 = [0.013256455593648; -0.016216516507728; 0.004041602572279]; % 3 months!
%         dV2 = [-7.803777280688135e-04; 0.001854569833090;-0.007247538179753]; 
        dV3 = [0.002544242144491; -0.002921527856874; 0.007703415162441];
        dV4 = [-0.001625936670348; -0.003125208256016; -0.008088501084076];
        dV5 = [-0.002918536114165;-0.003384664726700;0.008333531253574];
        dV6 = [-0.002355831016669;-0.002402859984804;-0.008729136657509];
        
        dV1 = [0.0132624893374001;-0.0162513152365047;0.00404609981491367];
        dV2 = [-0.000780457538397596;0.00185525021798016;-0.00724844974057293];
        
        % Shows the consecutive number of the maneuver applied
        maneuver_number = 1;
        
        n_integrations = 2;
        % Adjust maneuvers manually
        % 3 integrations means that one should start from maneuver 3
        deltaVs = {dV2;dV1;[0;0;0];};
        
        % Keep track on integration number. Don't change!
        n = 1;
        
        while ~complete
            
           % phi0 = reshape(eye(6), 36, 1);
           % init_state = [init_state(1:6); phi0];
            
            options = odeset('Events',@event_handler, 'MaxStep', -2700, 'InitialStep', -60);
            solution = ode113(@full_force_model,[init_t final_point],init_state,options);
            epochs = solution.x;
            rev_orbit_ode113 = solution.y;
            rev_ode113_y0state_E = rev_orbit_ode113(:,end); %solution.ye;
            % Event handler does inner transformation to check y=0
            % Now I need to transform the result to L2 frame
            rev_orbit_ode113 = EcenToL2frame( rev_orbit_ode113, epochs );
            % Now subtract the maneuver
            
            maneuver = deltaVs{maneuver_number};
            rev_ode113_y0state_E(1:6) = rev_ode113_y0state_E(1:6) - [0;0;0;maneuver(1);maneuver(2);maneuver(3)];
            
            reverse_orbit_ODE113 = [reverse_orbit_ODE113, rev_orbit_ode113];
            reverse_totalepochs_ode113 = [reverse_totalepochs_ode113, epochs];
            

            % Change the values for the next orbit part integration
            init_t = epochs(end);
            init_state = rev_ode113_y0state_E;     
            
            maneuver_number = maneuver_number + 1;
            
            n = n+1;
            if n > n_integrations 
                
                 % Before this run a small piece from last y = 0 till state at t = t_initial
               start_time = 958.910668311133e+006;

            options = odeset('MaxStep', -2700, 'InitialStep', -60);
            solution = ode113(@full_force_model,[init_t start_time],init_state,options);
            epochs = solution.x;
            rev_orbit_ode113 = solution.y;
            rev_ode113_y0state_E = rev_orbit_ode113(:,end);
            rev_orbit_ode113 = EcenToL2frame( rev_orbit_ode113, epochs );
              
                reverse_orbit_ODE113  = [reverse_orbit_ODE113, rev_orbit_ode113];
                reverse_totalepochs_ode113 = [reverse_totalepochs_ode45, epochs];
                
               complete = true; 
            end
            
        end

            %rkv89_conditions_difference = abs(fliplr(reverse_orbit_RKV_89) - orbit_RKV_89);
           	ode113_flp = fliplr(reverse_orbit_ODE113);
            ode113_initial_value_difference = abs(ode113_flp(1:6,1) - orbit_ODE113(1:6,1));
            disp('difference ODE113');
            disp(ode113_initial_value_difference);

         
     end
     
     
    %% Create a table with results
    abm_initial_value_difference = zeros(6,1);

    Integrators = {'RKV89';'ABM';'RK45';'PD78';'RKV89_embedded';'ODE113'};
    
    Init_diffs = [rkv89_initial_value_difference;abm_initial_value_difference;ode45_initial_value_difference;ode87_initial_value_difference;rkv89emb_initial_value_difference;ode113_initial_value_difference];
    %x_diffs = [rkv89_initial_value_difference;abm_initial_value_difference;rk_initial_value_difference;pd78_initial_value_difference ];
    dX = [rkv89_initial_value_difference(1);abm_initial_value_difference(1);ode45_initial_value_difference(1);ode87_initial_value_difference(1);rkv89emb_initial_value_difference(1);ode113_initial_value_difference(1)];
    dY = [rkv89_initial_value_difference(2);abm_initial_value_difference(2);ode45_initial_value_difference(2);ode87_initial_value_difference(2);rkv89emb_initial_value_difference(2);ode113_initial_value_difference(2)];
    dZ = [rkv89_initial_value_difference(3);abm_initial_value_difference(3);ode45_initial_value_difference(3);ode87_initial_value_difference(3);rkv89emb_initial_value_difference(3);ode113_initial_value_difference(3)];
    dVX = [rkv89_initial_value_difference(4);abm_initial_value_difference(4);ode45_initial_value_difference(4);ode87_initial_value_difference(4);rkv89emb_initial_value_difference(4);ode113_initial_value_difference(4)];
    dVY = [rkv89_initial_value_difference(5);abm_initial_value_difference(5);ode45_initial_value_difference(5);ode87_initial_value_difference(5);rkv89emb_initial_value_difference(5);ode113_initial_value_difference(5)];
    dVZ = [rkv89_initial_value_difference(6);abm_initial_value_difference(6);ode45_initial_value_difference(6);ode87_initial_value_difference(6);rkv89emb_initial_value_difference(6);ode113_initial_value_difference(6)];
    dX_scalar = [sqrt(dX(1)^2+dY(1)^2+dZ(1)^2);sqrt(dX(2)^2+dY(2)^2+dZ(2)^2);sqrt(dX(3)^2+dY(3)^2+dZ(3)^2);sqrt(dX(4)^2+dY(4)^2+dZ(4)^2);sqrt(dX(5)^2+dY(5)^2+dZ(5)^2);sqrt(dX(6)^2+dY(6)^2+dZ(6)^2)];
    dVX_scalar = [sqrt(dVX(1)^2+dVY(1)^2+dVZ(1)^2);sqrt(dVX(2)^2+dVY(2)^2+dVZ(2)^2);sqrt(dVX(3)^2+dVY(3)^2+dVZ(3)^2);sqrt(dVX(4)^2+dVY(4)^2+dVZ(4)^2);sqrt(dVX(5)^2+dVY(5)^2+dVZ(5)^2);sqrt(dVX(6)^2+dVY(6)^2+dVZ(6)^2)];
    Table = table(dX,dY,dZ,dVX,dVY,dVZ,dX_scalar,dVX_scalar,'RowNames',Integrators);
end

%% Plotting

% figure(2)
% plot3(reverse_orbit_RKV_89_emb(1,:),reverse_orbit_RKV_89_emb(2,:),reverse_orbit_RKV_89_emb(3,:),'g'); % orbit


% figure(1)
% view(3)
% grid on
% hold on
% plot3(0,0,0,'*r'); % nominal L2 point
% if RKV_89
%    plot3(orbit_RKV_89(1,:),orbit_RKV_89(2,:),orbit_RKV_89(3,:),'m'); % orbit
%    plot3(reverse_orbit_RKV_89(1,:),reverse_orbit_RKV_89(2,:),reverse_orbit_RKV_89(3,:),'b'); % orbit
% end


% figure(2)
% view(3)
% grid on
% hold on
% plot3(0,0,0,'*r'); % nominal L2 point
% if RKV_89_emb
%     plot3(orbit_RKV_89_emb(1,:),orbit_RKV_89_emb(2,:),orbit_RKV_89_emb(3,:),'m'); % orbit
%     plot3(reverse_orbit_RKV_89_emb(1,:),reverse_orbit_RKV_89_emb(2,:),reverse_orbit_RKV_89_emb(3,:),'b'); % orbit
% 
% end
% 
% figure(3)
% view(3)
% grid on
% hold on
% plot3(0,0,0,'*r'); % nominal L2 point
% if ODE87 
%     plot3(orbit_ODE87(1,:),orbit_ODE87(2,:),orbit_ODE87(3,:),'m');
%     plot3(reverse_orbit_ODE87(1,:),reverse_orbit_ODE87(2,:),reverse_orbit_ODE87(3,:),'b'); % orbit
% 
% end
% 
% 
% figure(4)
% view(3)
% grid on
% hold on
% plot3(0,0,0,'*r'); % nominal L2 point
% if ODE45
%     plot3(orbit_ODE45(1,:),orbit_ODE45(2,:),orbit_ODE45(3,:),'m');
%     plot3(reverse_orbit_ODE45(1,:),reverse_orbit_ODE45(2,:),reverse_orbit_ODE45(3,:),'b'); % orbit
% 
% end
% 
% figure(5)
% view(3)
% grid on
% hold on
% plot3(0,0,0,'*r'); % nominal L2 point
% if ODE113
%     plot3(orbit_ODE113(1,:),orbit_ODE113(2,:),orbit_ODE113(3,:),'m');
%     plot3(reverse_orbit_ODE113(1,:),reverse_orbit_ODE113(2,:),reverse_orbit_ODE113(3,:),'b'); % orbit
% 
% end


figure(6)
view(3)
grid on
hold on
plot3(0,0,0,'*r'); % nominal L2 point
%plot3(orbit_RKV_89(1,:),orbit_RKV_89(2,:),orbit_RKV_89(3,:),'b'); % orbit
plot3(orbit_RKV_89_emb(1,:),orbit_RKV_89_emb(2,:),orbit_RKV_89_emb(3,:),'b'); % orbit
%plot3(orbit_ODE87(1,:),orbit_ODE87(2,:),orbit_ODE87(3,:),'g');
%plot3(orbit_ODE45(1,:),orbit_ODE45(2,:),orbit_ODE45(3,:),'y'); % orbit
plot3(orbit_ODE113(1,:),orbit_ODE113(2,:),orbit_ODE113(3,:),'m');

difference = orbit_RKV_89_emb(:,:) - orbit_ODE113(:,1:length(orbit_RKV_89_emb));
difX = difference(1,:);
difY = difference(2,:);
difZ = difference(3,:);
difVX = difference(4,:);
difVY = difference(5,:);
difVZ = difference(6,:);

epochs = totalepochs_rkv89_emb;

% figure(7)
% subplot(2,1,1)
% grid on
% hold on
% plot(epochs,difX,'r');
% plot(epochs,difY,'g');
% plot(epochs,difZ,'b');
% subplot(2,1,2)
% grid on
% hold on
% plot(epochs,difVX,'r');
% plot(epochs,difVY,'g');
% plot(epochs,difVZ,'b');

figure(7)
grid on
hold on
plot(epochs,difX,'r');
plot(epochs,difY,'g');
plot(epochs,difZ,'b');
figure(8)
grid on
hold on
plot(epochs,difVX,'r');
plot(epochs,difVY,'g');
plot(epochs,difVZ,'b');

%% Plots info
figure(1)
title(['HALO orbit around L2 SEM. RKV89', model]);
legend('Nominal L2 point','HALO orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(2)
title(['HALO orbit around L2 SEM. RKV89emb', model]);
legend('Nominal L2 point','HALO orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(3)
title(['HALO orbit around L2 SEM. ODE87', model]);
legend('Nominal L2 point','HALO orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(4)
title(['HALO orbit around L2 SEM. ODE45', model]);
legend('Nominal L2 point','HALO orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

figure(5)
title(['HALO orbit around L2 SEM. ODE113', model]);
legend('Nominal L2 point','HALO orbit');
xlabel('x');
ylabel('y');
zlabel('z');
grid on

% figure(6)
% title('HALO orbit around SEM L2. Full force model');
% legend('Nominal L2 point','RKV8','RKV89 Embedded','PD78','ODE45','ODE113');
% set(gca,'fontsize',13)
% xlabel('x, km');
% ylabel('y, km');
% zlabel('z, km');
% grid on



figure(6)
title('HALO orbit around SEM L2. Full force model');
legend('Nominal L2 point','RKV89 Normal Maneuvers','RKV89 ODE113 Maneuvers','Location','northeast');
set(gca,'fontsize',13)
xlabel('x, km');
ylabel('y, km');
zlabel('z, km');
grid on


% figure(7)
% subplot(2,1,1)
% title('Difference in position components.');
% legend('X', 'Y', 'Z');
% set(gca,'fontsize',16)
% ylabel('Difference, km');
% xlabel('Epoch, sec');
% subplot(2,1,2)
% title('Differences in velocity components.');
% legend('VX', 'VY', 'VZ');
% set(gca,'fontsize',16)
% ylabel('Difference, km/s');
% xlabel('Epoch, sec');

figure(7)
title('Difference in position components.');
legend('X', 'Y', 'Z');
set(gca,'fontsize',20)
ylabel('Difference, km');
xlabel('Epoch, sec');
figure(8)
title('Differences in velocity components.');
legend('VX', 'VY', 'VZ');
set(gca,'fontsize',20)
ylabel('Difference, km/s');
xlabel('Epoch, sec');

