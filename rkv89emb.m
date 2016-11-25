function [epoch, output_state] = rkv89emb(f, t_range, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

if ~checkrkv89_emb
    %% Normal Integration
    while t < tfinal
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