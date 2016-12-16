function [epoch, y0state, output_state, last_point_in_E] = simple_rkv89emb_maneuvers(f, t_range, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    y0state = zeros(42,1);


    stop = false;

    global L2frame;
    global rkv89emb_lastpiece;

    t = t_range(1); %Initial epoch
    tfinal = t_range(length(t_range));
        
    % Preallocation
    output_state = [];
    epoch = [];
    last_point_in_E = [];

    % Set the first point
    output_state(:,1) = y;
    epoch(1) = t;
  
    % intermediary Earth frame state
    E_output_state(:,1) = y;
    
 
    
    % Settings
    max_attempts = 50;
    tolerance = 1e-11; % 1e-13;
    minimumStep = 1e-13;
    maximumStep = 2700;
    sigma = 0.9;
    
    global RKV_89_emb_check;
    
   
    
    if ~RKV_89_emb_check
        initial_step = 60;
    else
        initial_step = -60;
    end
    
     if rkv89emb_lastpiece
        %initial_step = -0.1;
        %initial_step = -0.1;
     end
    
    
    stepSize = initial_step;
%
%while ~stop



    %% Normal Integration
    
   
    
    while ~stop
%             if (t + stepSize) > tfinal
%                 stepSize = tfinal - t;
%             end

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
                
                    phi = reshape(state(7:end), 6, 6);
                    phi = xform*phi*xform^(-1);
                    phi = reshape(phi, 36,1);
                    %phi = state(7:end);
                    L2state = [L2state; phi];
                
                output_state = [output_state, L2state];   
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
                epoch = [epoch, t];
                
                % Now do the checking
                
                skip = 10;
%                 if RKV_89_emb_check 
%                    skip = 10; 
%                 end
               % skip = 3500; % Skip first y=0 crossing
             if size(output_state,2) > skip % skip first points
                 
                % New condition! Operate it manually
                %intermediate_maneuver = true;
%                 if intermediate_maneuver 
%                    if t > t_range(1) + 3.942e+6 % approx 1.5 month
%                       stop = true;
%                       break
%                    end
%                  
%                 end
                 
              %if ~intermediate_maneuver   
                 
                if ~isequal(sign(output_state(2,end-1)), sign(L2state(2,1)))  %&& t >= t_range+7.884e+6;%&& (L2state(1,1) < 0)
                   %7.845e+6  
                   ytol = 1e-5;
                    
                   [desired_t_for_maneuver, state_at_desired_t, state_at_desired_t_E ] = find_T_foryzero( [epoch(end-1) epoch(end)], E_output_state(:,end-1), ytol);                  
                   %output_state = [output_state, state_at_desired_t];
                   output_state(:,end) = state_at_desired_t;
                   epoch(end) = desired_t_for_maneuver;
                   %epoch = [epoch, desired_t_for_maneuver];
                   last_point_in_E = state_at_desired_t_E;
                   y0state = state_at_desired_t;
                   
                   stop = true;
                   break;
                    
                    
                end  
             %end
                
            end 
                
                
            else  % Earth-centered frame
                output_state = [output_state, state];% , - column ; - row
                E_output_state = [E_output_state, state]; 
                last_point_in_E = state;
                epoch = [epoch, t];
                
            end
            
            
        end
          if RKV_89_emb_check  
                if rkv89emb_lastpiece 
                  if isempty(t_range(2))
                     disp('Provide second epoch!'); 
                  end
                  tfinal = t_range(2); %must be provided as a second argument

                      if (t + stepSize) < tfinal 
                             stepSize = tfinal - t; 


                          [errhh, state] = RungeKutta89_2(f,y,t,stepSize);

                          if L2frame

                            % Subract coordinates of L2!
                            L2point = cspice_spkezr('392', t, 'J2000', 'NONE', '399');
                            conv_state = state;
                            conv_state(1:6) = state(1:6) - L2point;

                            xform = cspice_sxform('J2000','L2CENTERED', t);
                            L2state = xform*conv_state(1:6);

                                phi = reshape(state(7:end), 6, 6);
                                phi = xform*phi*xform^(-1);
                                phi = reshape(phi, 36,1);
                                L2state = [L2state; phi];

                            output_state = [output_state, L2state];   
                            E_output_state = [E_output_state, state]; 
                            last_point_in_E = state;
                            epoch = [epoch, t];
                          end

                           stop = true;
                           break
                           % break out of the loop!
                      end
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
