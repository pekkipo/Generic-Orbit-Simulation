function [ desired_t_for_maneuver, state_at_desired_t ] = find_T_foryzero( initials, init_state, ytol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Check the interpolation code
        % using binary search tree
        found = false;
        int_step = 0.1;
        initials = initials(1):int_step:initials(length(initials));%epochs(5872):int_step:epochs(5873);
        %init_state = orbit_rkv89_emb(:,5872);
        %ytol = 0.000001;
        desired_t_for_maneuver = 0;
        state_at_desired_t = zeros(6,1);
        yvalue = 0; % Desired value of y-component of the sat in L2centered frame
        while ~found
            %options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep', 0.001,'InitialStep',0.001);
            %options = odeset('MaxStep', 1,'InitialStep',0.1);
            [ti, oiE] = ode45(@force_model, initials, init_state);
            ti = ti';
            oiE = oiE';
            oiE = [oiE;ti];
            % now in this oi array I have to check second row to find the
            % closest to 0 +- tolerance
            
            % Convert to L2centered
            xform = cspice_sxform('J2000','L2CENTERED', ti);
            for g = 1:length(oiE)
                oi(1:6,g) = xform(:,:,g)*oiE(1:6,g);
            end
           
            
            center_epoch = floor(length(oi)/2); % integer epoch
            disp(center_epoch);
            center_state = oi(1:6,center_epoch);
            disp(center_state);
            ycenter = oi(2,center_epoch);
            center_t = oiE(7,center_epoch);
            disp(ycenter);
            disp(center_t);
            
            if ycenter > yvalue 
                initials = [initials(1) center_t]; 
                disp('bigger');
            end
            
            if ycenter < yvalue 
                initials = [center_t initials(length(initials))]; 
                init_state = center_state;
                disp('smaller');
            end
            
             left_border = yvalue - ytol;
             right_border = yvalue + ytol;
             
             if ycenter <= right_border && ycenter >= left_border
                 %index = find(abs(oi(2,:))<ytol); 
                 [closest_value, N] = min((abs(oi(2,:))));
                 desired_t_for_maneuver = oiE(7,N);
                 %state_at_desired_t = oi(1:6,N); % If I wanted L2frame
                 state_at_desired_t = oiE(1:6,N);
                 found = true;       
             end
            
        end
end


