function y0state = rkv89emb_maneuvers(f, t_range, y, numb, maneuvers_implementation)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    y0state = [];

    %%%% Maneuvers info
    t_tolerance = 1e-9;
    possible_t_for_maneuver1 = 9.747581737552394e+08; 
    
    %possible_t_for_maneuver2 = 9.902502994711838e+08;
    possible_t_for_maneuver2 = 9.823832611683108e+08;
    maneuvers = [possible_t_for_maneuver1, possible_t_for_maneuver2];
    % use vector of possible ts for manevuers
    % when this value is reached - run checking function
    % when the t* and state* are found stop the integrator
    stop = false;
    %21 Nov 2030 09:06:53.955
    % first is somewhere between epoch Nr. 5872 = 9.747554737552394e+08 sec and 
    % end Nr. 5873 = 9.747581737552394e+08
    %%%%
    
    global L2frame;

    t = t_range(1); %Initial epoch
    tfinal = t_range(length(t_range));
        
    % Preallocation
    output_state = [];
    epoch = [];
    last_point_in_E = [];
    
%     phi0 = [1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1];
%     y = [y;phi0];
    
    % Set the first point
    output_state(:,1) = y;
    epoch(1) = t;
  
    % intermediary Earth frame state
    E_output_state(:,1) = y;
    
    % Settings
    global checkrkv89_emb; % If true then the reverse check will be run
    max_attempts = 50;
    tolerance = 1e-11; % 1e-13;
    minimumStep = 1e-13;
    maximumStep = 2700;
    sigma = 0.9;
    
    if ~checkrkv89_emb
        initial_step = 60;
    else
        initial_step = -60;
    end
    
    stepSize = initial_step;
%
%while ~stop

% FOR NOW
checkrkv89_emb = false;

    
if ~checkrkv89_emb
    %% Normal Integration
    
   
    
    while t < tfinal && ~stop
            if (t + stepSize) > tfinal
                stepSize = tfinal - t;
            end

            goodStepTaken = false;
            currentAttempts = 0;

            while not(goodStepTaken)

                % Check stepSize
                if (abs(stepSize) < minimumStep) %&& not(finalStep))
                   if stepSize > 0.0 
                    stepSize = minimumStep;
                    else
                    stepSize = -minimumStep;
                   end
                end

                 if (abs(stepSize) > maximumStep)
                    if stepSize > 0.0 
                    stepSize = maximumStep;
                 else
                    stepSize = -maximumStep;
                    end
                 end
                 
                 % Check if with this step size we don't jump over the
                 % maneuver time
%                  m = maneuvers(1);
%                  for k = 1:length(maneuvers)
%                      potential_m = m;
%                      if t+stepSize > maneuvers(k)
%                         potential_m = maneuvers(k); 
%                      end
%                      
%                      if potential_m >= m
%                         m = potential_m;
%                      end
%   
%                  end
                 
                 %
                 %if ((t+stepSize) >= possible_t_for_maneuver1 && (t+stepSize) < possible_t_for_maneuver2)  || ((t+stepSize) >= possible_t_for_maneuver2)
                
             if maneuvers_implementation
                 if numb == 2
                    possible_t_for_maneuver1 = possible_t_for_maneuver2;
                 end
                 
                 if (t+stepSize) >= possible_t_for_maneuver1  
                   ytol = 1e-6;
                   
                   % here should add the check if sign of y changes from t
                   % to t + stepsize. Otherwise, should take one of the
                   % previous steps as t
                   
                   % I check only for OVERSHOOTING
                   % Undershooting is not considered. TODO
                   
                    [er, possible_new_state] = RungeKutta89_2(f,y,t,stepSize);
                    xform = cspice_sxform('J2000','L2CENTERED', t+stepSize);
                    possible_new_stateL2 = xform*possible_new_state(1:6);

%                     if maneuvers_implementation
%                         phi = reshape(possible_new_state(7:end), 6, 6);
%                         phi = xform*phi*xform^(-1);
%                         phi = reshape(phi, 36,1);
%                         possible_new_stateL2 = [possible_new_stateL2; phi];
%                     end
                    
                    % determine the direction, from + to - or the other way
                    % round
                    syms from_plus_to_minus;
                    if output_state(2,1) > 0
                        from_plus_to_minus = true;
                    else
                        from_plus_to_minus = false;
                    end
                    % if both less than zero - move previous t to the right
%                    if L2state(2) > 0 && possible_new_stateL2 > 0
%                          ind = epoch(length(epoch));
%                          found = false;
%                        while ~found
%                             check_y = output_state(ind);
%                             if check_y < 0 
%                                 found = false;
%                                 ind = ind - 1;
%                             else
%                                 found = true;
%                             end
%                        end
%                        % When found, set this state to right and it will be
%                        % fed into next checking function
%                        state = output_state(ind);
%                        t = epoch(ind);
%                    end
                   
                   % if both less than zero - move previous t to the left
                   if from_plus_to_minus == true
                       y1 = L2state(2,1);
                       y2 = possible_new_stateL2(2,1);
                       if (y1 < 0) && (y2 < 0)
                           ind = length(epoch);
                           found = false;
                           while ~found
                                check_y = output_state(2,ind);
                                if check_y < 0 
                                    found = false;
                                    ind = ind - 1;
                                else
                                    found = true;
                                end
                           end
                           % When found, set this state to right and it will be
                           % fed into next checking function
                           state = E_output_state(:,ind);
                           t = epoch(ind);
                       end
                       
                   else
                            y1 = L2state(2,1);
                            y2 = possible_new_stateL2(2,1);
                            if (y1 > 0) && (y2 > 0)
                               ind = length(epoch);
                               found = false;
                               while ~found
                                    check_y = output_state(2,ind);
                                    if check_y > 0 
                                        found = false;
                                        ind = ind - 1;
                                    else
                                        found = true;
                                    end
                               end
                               % When found, set this state to right and it will be
                               % fed into next checking function
                               state = E_output_state(:,ind);
                               t = epoch(ind);
                           end

                   end
                   % L2state - wrong, give state as I need in it in
                   % Earth-centerd. Convert within the func
                   % t+stepSize - should use possible t man instead
                   [desired_t_for_maneuver, state_at_desired_t, state_at_desired_t_E ] = find_T_foryzero( [t possible_t_for_maneuver1], state, ytol);                  
                   output_state = [output_state, state_at_desired_t];
                   epoch = [epoch, desired_t_for_maneuver];
                   last_point_in_E = state_at_desired_t_E;
                   y0state = state_at_desired_t;
                   stop = true;
                   break; 
                 end
                 
            end
                 
                 % Calculate raw state and the error (column vector)
                 [errh, state] = RungeKutta89_2(f,y,t,stepSize);
                 % Find max error component
                 error = maxerror(errh(1:6), state(1:6), y(1:6)); %1:6 because need to ignore I 36x1 matrix
                 currentStep = stepSize;
                
                 % Check the error
                 if (error ~= 0.0)

                     if (error > tolerance)    

                         stepSize = sigma * stepSize * ((tolerance/error)^(1/8));
                            if (abs(stepSize) < minimumStep)
                               if stepSize < 0.0 
                                    stepSize = -minimumStep;
                               else
                                    stepSize = minimumStep;
                               end
                            end

                            if (abs(stepSize) > maximumStep)
                               if stepSize > 0.0 
                                    stepSize = maximumStep;
                               else
                                    stepSize = -maximumStep;
                               end
                            end
                        currentAttempts = currentAttempts+1;
                     else 
                        
                           y = state;
                           t = t+currentStep;
                         
                           stepSize = sigma * stepSize * ((tolerance/error)^(1/9)); 
                           currentAttempts = 0;

                           if (abs(stepSize) > maximumStep)
                               if stepSize > 0.0 
                                    stepSize = maximumStep;
                               else
                                    stepSize = -maximumStep;
                               end
                           end

                           if (abs(stepSize) < minimumStep) %&& not(finalStep))
                                if stepSize > 0.0 
                                    stepSize = minimumStep;
                                else
                                    stepSize = -minimumStep;
                                end
                           end

                           goodStepTaken = true;
                    end
               %%%%

                 else  
                     y = state;
                     t = t+currentStep;
                     currentAttempts = 0;
                     goodStepTaken = true;
                end

                 if (currentAttempts >= max_attempts)
                       disp('Bad step');
                       y = state;
                       t = t+currentStep;
                       goodStepTaken = true; % actually no, but I have to leave the while loop
                 end

            end
            
             % Here add the stuff about maneuvers
            
        if ~stop     
            % Convert state into L2-centered frame if needed
            if L2frame
                xform = cspice_sxform('J2000','L2CENTERED', t);
                L2state = xform*state(1:6);
                if maneuvers_implementation
%                     phi = reshape(state(7:end), 6, 6);
%                     phi = xform*phi*xform^(-1);
%                     phi = reshape(phi, 36,1);
                    phi = state(7:end);
                    L2state = [L2state; phi];
                end
                output_state = [output_state, L2state];   
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
                
            else  % Earth-centered frame
                output_state = [output_state, state];% , - column ; - row
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
            end
            
            epoch = [epoch, t];
        end
            
            
            % Here goodstep is taken, and solution is ready
            

    end
    
    % Convert also the first point to L2
       if L2frame
       xform = cspice_sxform('J2000','L2CENTERED', t_range(1));
       output_state(1:6,1) = xform*output_state(1:6,1);
       end
    
end

end