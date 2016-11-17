function [output_solution, newstep] = Embedded_Verner89_new(f, t, y, h, tmax, tolerance)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   err_exponent = 1.0 / 8.0;


   % Verify that the step size is positive and that the upper endpoint
   % of integration is greater than the initial endpoint.              

   if (tmax < t) || (h <= 0.0)
      % flag = -2;
       output_solution = 0;
       disp('ERROR!')
       return
   end
   
   % If the upper endpoint of the independent variable agrees with the 
   % initial value of the independent variable.  Set the value of the  
   % dependent variable and return success.                            

   %hnext = h;
   newstep = h;
   output_solution = y;
   
   if (tmax == t)
      % flag = 0;
       return
   end
  
    ATTEMPTS = 12;

    sigma = 0.9;
    decPower = 1/7;
    minimumstep = 2.6964e+03;
    
    temp_y = y;
    while ( t+h < tmax) 
        for i = 0:1:ATTEMPTS
            disp(i);
            
             if ( t + h + 0.5 * h > tmax )
                 h = 0.5 * (tmax - t);
             end
            
             [errvect, err,  solution] = RungeKutta89(f, temp_y, t, h); % Solution and estimate
             error = abs(max(errvect));
             
             if error > tolerance
                 h = sigma * h * (tolerance/error)^decPower;
                 
                 if (fabs(h) < minimumstep)
                    h = minimumstep;
                    if i < ATTEMPTS
                        continue
                    else
                        disp('no more attemts');
                        newstep = h;
                        temp_y = solution; 
                        output_solution = temp_y;
                        return
                    end
                 end
                 
             else 
                 h = sigma * h * (tolerance/error)^err_exponent;
                 newstep = h;
                 temp_y = solution; 
                 output_solution = temp_y;
                 return
             end
             
        end
        
        
    end
    
   % flag = 0;
   % newstep = h;
   % output_solution = solution;

end

