function [y, tour] = RKV89(f,t,y0,n)

%h = step;
y(:,1) = y0;
s6 = sqrt(6);
tour(1) = t(1);
global epochs_numbers;
global maneuvers;
tolerance = 1e-13; % not sure

%n_et = 9120-3244; % 3244 - number of epoch before HALO orbit starts
%maneuver1 = [0.003696355989169846;-0.004709746685339394;0.01461216953990576]; 
% my new maneuver
%maneuver1 = [0.0016;-0.0013;0.0011]; % my recalculated maneuver



%corrected_step = t(2)-t(1);

    for i = 1 : n-1
     
        h = t(i+1) - t(i);
       % h = corrected_step;
        
        
   %% Value estimation  
   % For instance: a5 = 2 + 2*s6/15
   
   a2 = 1.0 / 12.0;
   a3 = 1.0 / 9.0;
   a4 = 1.0 / 6.0;
   a5 = 2.0 * (1.0 + s6) / 15.0;
   a6 = (6.0 + s6) / 15.0;
   a7 = (6.0 - s6) / 15.0;
   a8 = 2.0 / 3.0;
   a9 = 1.0 / 2.0;
   a10 = 1.0 / 3.0;
   a11 = 1.0 / 4.0;
   a12 = 4.0 / 3.0;
   a13 = 5.0 / 6.0;
   a15 = 1.0 / 6.0;
   
   
   % For instance:
   % h * ( b51*k1 + b53*k3 + b54*k4) ) = ..
   % ..= h/375*((4+94s6)*k1 - (282+252s6)*k3 + (328+208s6)*k4))
   % b51 = (4.0 + 94.0 * s6) / 375.0;
  
   
   b31 = 1.0 / 27.0;
   b32 = 2.0 / 27.0;
   b41 = 1.0 / 24.0;
   b43 = 3.0 / 24.0;
   b51 = (4.0 + 94.0 * s6) / 375.0;
   b53 = -(282.0 + 252.0 * s6) / 375.0;
   b54 = (328.0 + 208.0 * s6) / 375.0;
   b61 = (9.0 - s6) / 150.0;
   b64 = (312.0 + 32.0 * s6) / 1425.0;
   b65 = (69.0 + 29.0 * s6) / 570.0;
   b71 = (927.0 - 347.0 * s6) / 1250.0;
   b74 = (-16248.0 + 7328.0 * s6) / 9375.0;
   b75 = (-489.0 + 179.0 * s6) / 3750.0;
   b76 = (14268.0 - 5798.0 * s6) / 9375.0;
   b81 = 4.0 / 54.0;
   b86 = (16.0 - s6) / 54.0;
   b87 = (16.0 + s6) / 54.0;
   b91 = 38.0 / 512.0;
   b96 = (118.0 - 23.0 * s6) / 512.0;
   b97 = (118.0 + 23.0 * s6) / 512.0;
   b98 = - 18.0 / 512.0;
   b10_1 = 11.0 / 144.0;
   b10_6 = (266.0 - s6) / 864.0;
   b10_7 = (266.0 + s6) / 864.0;
   b10_8 = - 1.0 / 16.0;
   b10_9 = - 8.0 / 27.0;
   b11_1 = (5034.0 - 271.0 * s6) / 61440.0;
   b11_7 = (7859.0 - 1626.0 * s6) / 10240.0;
   b11_8 = (-2232.0 + 813.0 * s6) / 20480.0;
   b11_9 = (-594.0  + 271.0 * s6) / 960.0;
   b11_10 = (657.0 - 813.0 * s6) / 5120.0;
   b12_1 = (5996.0 - 3794.0 * s6) / 405.0;
   b12_6 = (-4342.0 - 338.0 * s6) / 9.0;
   b12_7 = (154922.0 - 40458.0 * s6) / 135.0;
   b12_8 = (-4176.0 + 3794.0 * s6) / 45.0;
   b12_9 = (-340864.0 + 242816.0 * s6) / 405.0;
   b12_10 = (26304.0 - 15176.0 * s6) / 45.0;
   b12_11 = -26624.0 / 81.0;
   b13_1 = (3793.0 + 2168.0 * s6) / 103680.0;
   b13_6 = (4042.0 + 2263.0 * s6) / 13824.0;
   b13_7 = (-231278.0 + 40717.0 * s6) / 69120.0;
   b13_8 = (7947.0 - 2168.0 * s6) / 11520.0;
   b13_9 = (1048.0 - 542.0 * s6) / 405.0;
   b13_10 = (-1383.0 + 542.0 * s6) / 720.0;
   b13_11 = 2624.0 / 1053.0;
   b13_12 = 3.0 / 1664.0;
   b14_1 = -137.0 / 1296.0;
   b14_6 = (5642.0 - 337.0 * s6) / 864.0;
   b14_7 = (5642.0 + 337.0 * s6) / 864.0;
   b14_8 = -299.0 / 48.0;
   b14_9 = 184.0 / 81.0;
   b14_10 = -44.0 / 9.0;
   b14_11 = -5120.0 / 1053.0;
   b14_12 = -11.0 / 468.0;
   b14_13 = 16.0 / 9.0;
   b15_1 = (33617.0 - 2168.0 * s6) / 518400.0;
   b15_6 = (-3846.0 + 31.0 * s6) / 13824.0;
   b15_7 = (155338.0 - 52807.0 * s6) / 345600.0;
   b15_8 = (-12537.0 + 2168.0 * s6) / 57600.0;
   b15_9 = (92.0 + 542.0 * s6) / 2025.0;
   b15_10 = (-1797.0 - 542.0 * s6) / 3600.0;
   b15_11 = 320.0 / 567.0;
   b15_12 = -1.0 / 1920.0;
   b15_13 = 4.0 / 105.0;
   b16_1 = (-36487.0 - 30352.0 * s6) / 279600.0;
   b16_6 = (-29666.0 - 4499.0 * s6) / 7456.0;
   b16_7 = (2779182.0 - 615973.0 * s6) / 186400.0;
   b16_8 = (-94329.0 + 91056.0 * s6) / 93200.0;
   b16_9 = (-232192.0 + 121408.0 * s6) / 17475.0;
   b16_10 = (101226.0 - 22764.0 * s6) / 5825.0;
   b16_11 = - 169984.0 / 9087.0;
   b16_12 = - 87.0 / 30290.0;
   b16_13 =  492.0 / 1165.0;
   b16_15 =  1260.0 / 233.0;
   
   

   k1 = f(t(i),y(:,i));                                                   
   k2 = f(t(i)+a2*h, y(:,i) + h*k1/12);                                  
   k3 = f(t(i)+a3*h, y(:,i)+ h*(b31*k1 + b32*k2));                           
   k4 = f(t(i)+a4*h, y(:,i) + h * ( b41*k1 + b43*k3));
   k5 = f(t(i)+a5*h, y(:,i) + h * (b51*k1 + b53*k3 + b54*k4));
   k6 = f(t(i)+a6*h, y(:,i) + h * (b61*k1 + b64*k4 + b65*k5));
   k7 = f(t(i)+a7*h, y(:,i) + h * (b71*k1 + b74*k4 + b75*k5 + b76*k6));
   k8 = f(t(i)+a8*h, y(:,i) + h * (b81*k1 + b86*k6 + b87*k7));
   k9 = f(t(i)+a9*h, y(:,i) + h * (b91*k1 + b96*k6 + b97*k7 + b98*k8));
   k10 = f(t(i)+a10*h, y(:,i) + h * (b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8 + b10_9*k9));
   k11 = f(t(i)+a11*h, y(:,i) + h * (b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10));
   k12 = f(t(i)+a12*h, y(:,i) + h * (b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8 + b12_9*k9 + b12_10 * k10 + b12_11 * k11));
   k13 = f(t(i)+a13*h, y(:,i) + h * (b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12));
   k14 = f(t(i)+h, y(:,i) + h * (b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8 + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13));
   k15 = f(t(i)+a15*h, y(:,i) + h * (b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8 + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13));
   k16 = f(t(i)+h, y(:,i) + h * (b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8 + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13 + b16_15*k15));
     
   
    
   % Convenient notation 
   
   c1 = 103.0 / 1680.0;
   c8 = -27.0 / 140.0;
   c9 = 76.0 / 105.0;
   c10 = -201.0 / 280.0;
   c11 = 1024.0 / 1365.0;
   c12 = 3.0 / 7280.0;
   c13 = 12.0 / 35.0;
   c14 = 9.0 / 280.0;
     
   % Turn this:  y(:,i+1) = y(:,i) +  h * ( 103/1680 * k1 - 27/140 * k8 + 76/105 * k9    
   %               - 201/280 * k10 + 1024/1365 * k11 + 3/7280 k12 + 12/35 k13  + 9/280 k14)
   % into this:
   
   
   
   % Calculate the value
   y(:,i+1) = y(:,i) +  h * (c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11 + c12 * k12 + c13 * k13 + c14 * k14);
   
        % Add maneuver if this is required epoch - 1. Output state for 9120
        % for instance. That means that i is the epoch 9119
        for k = 1:length(epochs_numbers)
            if i == epochs_numbers(k) - 1
                % If this epoch is one of the epoch presented in maneuvers
                % array - add dV to its components
                applied_maneuver = maneuvers{k};
                y(4,i+1) = y(4,i+1) + applied_maneuver(1);
                y(5,i+1) = y(5,i+1) + applied_maneuver(2);
                y(6,i+1) = y(6,i+1) + applied_maneuver(3);
            end
        end
        
        % If this is a reverse integation to check precision - subtract
      % maneuvers
      if checkrkv89 == true
           for k = 1:length(epochs_numbers)
                if i == epochs_numbers(k) - 1
                    % If this epoch is one of the epoch presented in maneuvers
                    % array - add dV to its components
                    applied_maneuver = maneuvers{k};
                    y(4,i+1) = y(4,i+1) - applied_maneuver(1);
                    y(5,i+1) = y(5,i+1) - applied_maneuver(2);
                    y(6,i+1) = y(6,i+1) - applied_maneuver(3);
                end
          end
      end
        
        
%      if i == n_et-1
%         y(4,i+1) = y(4,i+1) + maneuver1(1);
%         y(5,i+1) = y(5,i+1) + maneuver1(2);
%         y(6,i+1) = y(6,i+1) + maneuver1(3);
%     end
    
%    %% Error estimation
%    
% %    The error is estimated to be
% %    err = - h*( 1911 k1 - 34398 k8 + 61152 k9 - 114660 k10 + 114688 k11
% %    + 63 k12 + 13104 k13 + 3510 k14 - 39312 k15 - 6058 k16 / 109200
% %    The step size h is then scaled by the scale factor
% %    scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/8
% %    The scale factor is further constrained 0.125 < scale < 4.0.
% %    The new step size is h := scale * h
% %    
%     
%     e1 = -1911.0 / 109200.0;
%     e8 = 34398.0 / 109200.0;
%     e9 = -61152.0 / 109200.0;
%     e10 = 114660.0 / 109200.0;
%     e11 = -114688.0 / 109200.0;
%     e12 = -63.0 / 109200.0;
%     e13 = -13104.0 / 109200.0;
%     e14 = -3510.0 / 109200.0;
%     e15 = 39312.0 / 109200.0;
%     e16 = 6058.0 / 109200.0;
% 
%     error = h* (e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 + e12*k12 + e13*k13 + e14*k14 + e15*k15 + e16*k16);
%     error = [error(1); error(2); error(3)];
%     error = max(abs(error)); 
%     % max - because I need to take only value in column vector
%     
%     MIN_SCALE_FACTOR = 0.125;
%     MAX_SCALE_FACTOR = 10.0;%/4.0;
%     scale = [];
%     
%     
%     if (error == 0.0) 
%            scale = MAX_SCALE_FACTOR;
%            
%     else 
%             temp_y = [solution(1);solution(2);solution(3)];
%             temp_y = max(temp_y);
%              if temp_y == 0.0
%                  yy = tolerance;
%              else
%                  yy = abs(temp_y);
%              end
%              scale = 0.8 * (tolerance * yy /  error)^(1/8);
%              
%              if scale < MIN_SCALE_FACTOR 
%                 scale = MIN_SCALE_FACTOR;
%              else
%                 scale = scale;
%              end
%              
%              if scale < MAX_SCALE_FACTOR 
%                  scale = scale;
%              else 
%                  scale = MAX_SCALE_FACTOR;
%              end
%             
%     end
%      
%      corrected_step = corrected_step*scale;
%        
%        % h = h_new; % now the step for further integrations is changed to new value

    
    tour(i+1) = t(i) + h;
    
    end

end