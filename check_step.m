function step = check_step(stepSize, finalStep)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
minimumStep = 1e-13;
maximumStep = 2700;

 if (abs(stepSize) > maximumStep)
     if stepSize > 0.0 
        step = maximumStep;
     else
        step = -maximumStep;
     end
  end
               
   if (abs(stepSize) < minimumStep) && not(finalStep)
      if stepSize > 0.0 
         step = minimumStep;
      else
         step = -minimumStep;
      end
   end
               
end

