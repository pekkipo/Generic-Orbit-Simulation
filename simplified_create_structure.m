function [Earth, Sun, Moon] = simplified_create_structure( bodies, t, observer)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%global G;
% local G
G=6.67e-20; %km

    for pl=1:length(bodies) 
    
    field1 = 'name';  value1 = bodies{1,pl};
    field2 = 'x';  value2 = [];
    field3 = 'y';  value3 = [];
    field4 = 'z';  value4 = [];
    field5 = 'vx';  value5 = [];
    field6 = 'vy';  value6 = [];
    field7 = 'vz';  value7 = [];
    field8 = 'mass'; value8 = [];
    field9 = 'GM'; value9 = [];
    field10 = 'coords'; value10 = [];
        if pl == 1
        Earth = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Earth.GM = cspice_bodvrd( bodies{2,pl}, 'GM', 1 );
        Earth.mass = (Earth.GM)/G;
        Earth_coords = cspice_spkezr ( bodies{2,pl}, t, 'J2000', 'NONE', observer );
        Earth.x = Earth_coords(1);
        Earth.y = Earth_coords(2);
        Earth.z = Earth_coords(3);
        Earth.vx = Earth_coords(4);
        Earth.vy = Earth_coords(5);
        Earth.vz = Earth_coords(6);
        Earth.coords = Earth_coords(1:3);
        elseif pl == 2
        Sun = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Sun.GM = cspice_bodvrd( bodies{2,pl}, 'GM', 1 );
        Sun.mass = (Sun.GM)/G;  
        Sun_coords = cspice_spkezr ( bodies{2,pl}, t, 'J2000', 'NONE', observer );
        Sun.x = Sun_coords(1);
        Sun.y = Sun_coords(2);
        Sun.z = Sun_coords(3);
        Sun.vx = Sun_coords(4);
        Sun.vy = Sun_coords(5);
        Sun.vz = Sun_coords(6);
        Sun.coords = Sun_coords(1:3);
        elseif pl == 3
        Moon = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
        Moon.GM = cspice_bodvrd( bodies{2,pl}, 'GM', 1 );
        Moon.mass = (Moon.GM)/G;
        Moon_coords = cspice_spkezr ( bodies{2,pl}, t, 'J2000', 'NONE', observer );
        Moon.x = Moon_coords(1);
        Moon.y = Moon_coords(2);
        Moon.z = Moon_coords(3);
        Moon.vx = Moon_coords(4);
        Moon.vy = Moon_coords(5);
        Moon.vz = Moon_coords(6);
        Moon.coords = Moon_coords(1:3);
        end
    end

end