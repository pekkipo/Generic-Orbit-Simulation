% Checking GMAT data
clc
clear all
close all

METAKR = 'planetsorbitskernels.txt';
cspice_furnsh ( METAKR );

%% Define local variables
load('irassihalotime.mat', 'Date');
load('irassihalogmat.mat', 'Gmat');

et_vector = zeros(1,length(Date));
   for d=1:length(Date)
        utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
        et_vector(d) = cspice_str2et (utcdate);
   end

Gmat = EcenToL2frame( Gmat, et_vector );

%% Energy
figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'r')% 
