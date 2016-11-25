function [epoch, output_state] = rkv89emb(f, t_range, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %%%% Maneuvers info
    t_tolerance = 1e-9;
    possible_t_for_maneuver1 = 9.747581737552394e+08; 
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
        
    % Set the first point
    output_state(:,1) = y;
    epoch(1) = t;
    
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
                 if (t+stepSize) >= possible_t_for_maneuver1  
                   ytol = 1e-8;
                   [desired_t_for_maneuver, state_at_desired_t] = find_T_foryzero( [t t+stepSize], L2state, ytol);                  
                   output_state = [output_state, state_at_desired_t];
                   epoch = [epoch, desired_t_for_maneuver];
                   stop = true;
                   break; 
                 end
                 
                 % Calculate raw state and the error (column vector)
                 [errh, state] = RungeKutta89_2(f,y,t,stepSize);
                 % Find max error component
                 error = maxerror(errh, state, y); 
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
                L2state = xform*state;    
                output_state = [output_state, L2state];   
            else  % Earth-centered frame
                output_state = [output_state, state];% , - column ; - row
            end
            
            epoch = [epoch, t];
        end
            
            
            % Here goodstep is taken, and solution is ready
            

    end
    
    % Convert also the first point to L2
       if L2frame
       xform = cspice_sxform('J2000','L2CENTERED', t_range(1));
       output_state(:,1) = xform*output_state(:,1);
       end
    
end
%% Reverse check
if checkrkv89_emb == true
    while t > tfinal
            if (t + stepSize) < tfinal
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
                 error = maxerror(errh, state, y); 
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
            
            % Convert state into L2-centered frame if needed
            if L2frame
                xform = cspice_sxform('J2000','L2CENTERED', t);
                L2state = xform*state;
                output_state = [output_state, L2state];
                
                % here check function that will determine if y = 0 and t <
                % t_potential_maneuver
                % if t > t_potential_maneuver -> increase the step and make
                % goodStepTaken = false to run the loop again but with
                % smaller t!
                
                
            else  % Earth-centered frame
                output_state = [output_state, state];% , - column ; - row
            end
           
            epoch = [epoch, t];
            
            % Here goodstep is taken, and solution is ready
            
            
    end
    
        % Convert also the first point to L2
           if L2frame
           xform = cspice_sxform('J2000','L2CENTERED', t_range(1));
           output_state(:,1) = xform*output_state(:,1);
           end


    end

    %end
end