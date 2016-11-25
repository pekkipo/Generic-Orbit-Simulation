function orbitL2 = EcenToL2frame( orbit, et_vector )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    orbitL2 = zeros(6, length(orbit));
    % Convert to L2 frame
    xform = cspice_sxform('J2000','L2CENTERED', et_vector);
    
    for g = 1:length(orbit)
        orbitL2(:,g) = xform(:,:,g)*orbit(:,g);
    end


end
