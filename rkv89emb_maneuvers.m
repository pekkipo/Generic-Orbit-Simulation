function y0state = rkv89emb_maneuvers(f, t_range, y, numb, maneuvers_implementation)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    y0state = [];

    %%%% Maneuvers info
    t_tolerance = 1e-9;
    possible_t_for_maneuver1 = 9.747581737552394e+08; 
    possible_t_for_maneuver2 = 9.832987163183606e+08;%9.911611163179582e+08;%9.830395163183095e+08;%9.902502994711838e+08;%9.823832611683108e+08;
                                % FEB
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
                
                % Subract coordinates of L2!
                L2point = cspice_spkezr('392', t, 'J2000', 'NONE', '399');
                conv_state = state;
                conv_state(1:6) = state(1:6) - L2point;
                
                xform = cspice_sxform('J2000','L2CENTERED', t);
                L2state = xform*conv_state(1:6);
                if maneuvers_implementation
                    phi = reshape(state(7:end), 6, 6);
                    phi = xform*phi*xform^(-1);
                    phi = reshape(phi, 36,1);
                    %phi = state(7:end);
                    L2state = [L2state; phi];
                end
                output_state = [output_state, L2state];   
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
                epoch = [epoch, t];
                
                % Now do the checking

             if size(output_state,2) > 10 % skip first points
                if ~isequal(sign(output_state(2,end-1)), sign(L2state(2,1)))
                   
                    ytol = 1e-5;
                    
                   [desired_t_for_maneuver, state_at_desired_t, state_at_desired_t_E ] = find_T_foryzero( [epoch(end-1) epoch(end)], E_output_state(:,end-1), ytol);                  
                   output_state = [output_state, state_at_desired_t];
                   epoch = [epoch, desired_t_for_maneuver];
                   last_point_in_E = state_at_desired_t_E;
                   y0state = state_at_desired_t;
                   stop = true;
                   break;
                    
                    
                end    
                
            end 
                
                
            else  % Earth-centered frame
                output_state = [output_state, state];% , - column ; - row
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
                epoch = [epoch, t];
            end
            
            
        end
            
            
            % Here goodstep is taken, and solution is ready
            

    end
    
   % Convert also the first point to L2
       if L2frame
           L2point = cspice_spkezr('392', t_range(1), 'J2000', 'NONE', '399');
           conv_out1 = output_state(:,1);
           conv_out1(1:6) = output_state(1:6,1) - L2point;
           
       xform = cspice_sxform('J2000','L2CENTERED', t_range(1));
       output_state(1:6,1) = xform*conv_out1(1:6,1);
       end
    
end

end