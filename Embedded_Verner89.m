function [flag, output_solution, newstep] = Embedded_Verner89(f, t, y, h, tmax, tolerance)
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

   %hnext = h;
   newstep = h;
   output_solution = y;
   
   if (tmax == t)
       flag = 0;
       return
   end

   % Redefine the error tolerance to an error tolerance per unit   
   % length of the integration interval.                            

   tolerance = tolerance / (tmax - t);%10^(-13);%tolerance / (tmax - t);

    % Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to 
    % maintain an error less than tolerance * (xmax-x) using an     
    % initial step size of h and initial value: y = y[0]    
    ATTEMPTS = 12;
    MIN_SCALE_FACTOR = 0.125;
    MAX_SCALE_FACTOR = 4.0;
    
    temp_y = y;
    %while ( ~last_interval ) 
    while ( last_interval ~= 1 )
         for i=0:1:ATTEMPTS
             %i = i+1;
             last_interval = 0;
             if( t + h >= tmax )
                 h = tmax - t;
                 last_interval = 1; 
             elseif ( t + h + 0.5 * h > tmax )
                 h = 0.5 * (tmax - t);
             end
             [errvect,  solution] = RungeKutta89(f, temp_y, t, h);
             err = abs(max(errvect));
             if (err == 0.0) 
                 scale = MAX_SCALE_FACTOR;
                 disp('err = 0 break');
                 break 
             end
             % make scalar position
             temp_y_scalar = [temp_y(1);temp_y(2);temp_y(3)];%sqrt(temp_y(1)^2 + temp_y(2)^2 + temp_y(3)^2);
             temp_y_scalar = max(temp_y_scalar);
             if temp_y_scalar == 0.0
                 yy = tolerance;
             else
                 yy = abs(temp_y_scalar);
             end
             
             %% My code 
%              maxval = 2700;
%              minval = 10^(-13);
%              
%              if yy < minval
%              yy = minval;
%              else
%              yy = yy;
%              end
%              
%              if yy < maxval 
%                  yy = yy;
%              else 
%                  yy = maxval;
%              end
             %%
             
             scale = 0.8 * (tolerance * yy /  err)^err_exponent; 
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
                 %disp('err < tolerance');
                 break
             end
             
             h = h*scale;
         end
         % Not sure but think I should put it here
        % h = h*scale;
         
      if (i >= ATTEMPTS) 
          %hnext = h;
          flag = -1;
          newstep = h;
          return
      end
      
      temp_y = solution;  
      
      t = t + h;
    end
    
   %hnext = h;
   newstep = h;
   output_solution = solution;
   flag = 0;
%    % my stuff
%    
   minstep = (2.6964e+03);
  % maxstep = 2700;
   %86400/1000;%(2.6964e+03)/5;%30;
   
   if newstep < minstep
       newstep = minstep;
   end
   
   return
    
end

