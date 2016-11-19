function [output_state, stepTaken, final] = rkv(f, t,y, stepSize, tfinal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

max_attempts = 50;
tolerance = 1e-13; % 1e-13;

minimumStep = 1e-13;
maximumStep = 2700;
goodStepTaken = false; % 1 if good step was taken
currentAttempts = 0;
sigma = 0.9;
final = false;

temporary_raw_state = y;

  if (t + stepSize) > tfinal
      finalStep = true;
      final = true;
  else
      finalStep = false;
      
  end


%   if (abs(stepSize) < minimumStep) %&& not(finalStep))
%         if stepSize > 0.0 
%             stepSize = minimumStep;
%         else
%             stepSize = -minimumStep;
%         end
%   end
%   
%   if (abs(stepSize) > maximumStep)
%        if stepSize > 0.0 
%             stepSize = maximumStep;
%        else
%             stepSize = -maximumStep;
%        end
%   end
  
  % Start checking
  while not(goodStepTaken)
      
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
      
        disp(stepSize);
        [errh, state] = RungeKutta89(f,y,t,stepSize);
        if currentAttempts == 0
            temporary_raw_state = state;
        end
        error = maxerror(errh, state, y); % 3rd can be temporary_raw or y, don't know for sure yet
        %error = abs(max(errh));
        stepTaken = stepSize;
        
        if (error ~= 0.0)
           %%%%
           if (error > tolerance)
               
                stepSize = sigma * stepSize * (tolerance/error)^(1/7);

                if (abs(stepSize) < minimumStep)
                    if stepSize < 0.0 
                        stepSize = -minimumStep;
                    else
                        stepSize = minimumStep;
                    end

                    currentAttempts = currentAttempts+1;
                end 
           else % if error is okay, within the tolerance boundaries
               stepSize = sigma * stepSize * (tolerance/error)^(1/8);
               output_state = state;
               currentAttempts = 0;
               stepTaken = stepSize;
               goodStepTaken = true;
           end
           %%%%
           
        else  % 0.0 means no error control; in that case leave step alone
           %memcpy(outState, candidateState, dimension*sizeof(Real));
           output_state = state; 
           currentAttempts = 0;
           goodStepTaken = true;
        end
        
        if (currentAttempts >= max_attempts)
           %return false;
           goodStepTaken = true; % actually no, but I have to leave the while loop
           stepTaken = stepSize;
        end
      
  end


end