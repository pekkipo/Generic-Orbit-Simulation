% Checking GMAT data
clc
clear all
close all

%% Define local variables
load('irassihalotime.mat', 'Date')
load('irassihalogmat.mat', 'Gmat')


%% Energy
figure(1)
view(3)
grid on
hold on
plot3(Gmat(1,:),Gmat(2,:),Gmat(3,:),'r')% 
