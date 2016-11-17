load('IRASSIFullMissionDate.mat', 'Date');
load('IRASSIFullMission.mat', 'Gmat');

figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'b');% Reference