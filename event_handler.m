 function [value,isterminal,direction] = event_handler(t,y)
      % Event function -- y0 shared with the outer function.
      
      
      % First transform to L2
       L2point = cspice_spkezr('392', t, 'J2000', 'NONE', '399');
       conv_state = y;
       conv_state(1:6) = y(1:6) - L2point;
                
       xform = cspice_sxform('J2000','L2CENTERED', t);
       L2state = xform*conv_state(1:6);
%        phi = reshape(y(7:end), 6, 6);
%        phi = xform*phi*xform^(-1);
%        phi = reshape(phi, 36,1);
       %L2state = [L2state; phi];
                
       y = L2state;
    
      % y = y(2) - L2point(2);
      
     
      
      value = y(2); % second row which is Y component
      isterminal = 1;         % stop at local minimum
      direction  = 0;         % [local minimum, local maximum]
   end


