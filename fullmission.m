
METAKR = 'planetsorbitskernels.txt';
cspice_furnsh ( METAKR );

%% Define local variables
load('IRASSIFullMissionDate.mat', 'Date');
load('IRASSIFullMission.mat', 'Gmat');

et_vector = zeros(1,length(Date));
   for d=1:length(Date)
        utcdate = datestr((datetime(Date(d,:),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'UTC')), 'yyyy mmm dd HH:MM:SS.FFF');
        et_vector(d) = cspice_str2et (utcdate);
   end

gmE = Gmat;
Gmat = EcenToL2frame( Gmat(:,1:14866), et_vector(1:14866) );
l2 = gmE(1:3) - Gmat(1:3);


figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,1:14866),Gmat(2,1:14866),Gmat(3,1:14866),'b');% Reference
%plot3(Gmat(1),Gmat(2),Gmat(3),'*r');