function [flag, output_solution, newstep] = Embedded_Verner89(f, t, y, h, tmax, hnext, tolerance)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

   err_exponent = 1.0 / 8.0;
   %temp_y = zeros(0,0);
%    err;
%    yy;
%    i;
   last_interval = 0;
   %scale = 1;
   scale = [];

   % Verify that the step size is positive and that the upper endpoint
   % of integration is greater than the initial endpoint.              

   if (tmax < t) || (h <= 0.0)
       flag = -2;
       output_solution = 0;
       disp('ERROR!')
       return
   end
   
   % If the upper endpoint of the independent variable agrees with the 
   % initial value of the independent variable.  Set the value of the  
   % dependent variable and return success.                            

   hnext = h;
   solution = y;
   
   if (tmax == t)
       flag = 0;
       return
   end

   % Redefine the error tolerance to an error tolerance per unit   
   % length of the integration interval.                            

   tolerance = tolerance / (tmax - t);

    % Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to 
    % maintain an error less than tolerance * (xmax-x) using an     
    % initial step size of h and initial value: y = y[0]    
    ATTEMPTS = 12;
    MIN_SCALE_FACTOR = 0.125;
    MAX_SCALE_FACTOR = 4.0;
    
    temp_y = y;
    while ( ~last_interval ) 
         for i=0:1:ATTEMPTS
             last_interval = 0;
             if( t + h >= tmax )
                 h = tmax - t;
                 last_interval = 1; 
             elseif ( t + h + 0.5 * h > tmax )
                 h = 0.5 * (tmax - t);
             end
             [err, solution] = RungeKutta89(f, temp_y, t, h);
             err = max(abs(err));
             if (err == 0.0) 
                 scale = MAX_SCALE_FACTOR;
                 break 
             end
             if max(temp_y) == 0.0
                 yy = tolerance;
             else
                 yy = max(abs(temp_y));
             end
             scale = 0.8 * (tolerance * max(yy) /  err)^err_exponent; % max to reduce it to 1 dimension!
             % disp(scale);
             if scale < MIN_SCALE_FACTOR 
                scale = MIN_SCALE_FACTOR;
             else
                scale = scale;
             end
             
             if scale < MAX_SCALE_FACTOR 
                 scale = scale;
             else 
                 scale = MAX_SCALE_FACTOR;
             end
             
             if (err < (tolerance * yy))
                 break
             end
             
             h = h*scale;
         end
         
      if (i >= ATTEMPTS) 
          hnext = h;
          flag = -1;
          newstep = hnext;
          return
      end
      
      temp_y = solution;  
      
      t = t + h;
    end
    
   hnext = h;
   newstep = hnext;
   output_solution = solution;
   flag = 0;
   
   return
    
end

