%% Script for example 1: Bilinear academic example -- LMI approach
% Inputs: 
%   - none 
%
% Outputs: 
%   - none
%
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"

clear;clc;clearvars;
% Set MATLAB figure style
set(groot, 'defaultAxesColorOrder', get(gca,'colororder')); % Default color order
set(groot,'defaultAxesFontSize', 14); % Set font size
set(groot,'defaultLineLineWidth', 1.5); % Set line width

addpath("fcn\")
format long;

%% System dynamics
sys.A = [1,1;0,1];
sys.B0 = [1;1];
sys.Aux = [1,0;0,1];
sys.dynamics = @(x,u) sys.A*x + sys.B0*u + sys.Aux*kron(u,x);
%
[sys.n_x,sys.n_u] = size(sys.B0); % state and input dimension

%% data properties
param.sigma_x = 1; 
param.stdNoise = 0.1; 
param.delta = 0.05;
param.umax = 2; param.umin = -2;

%% Compute RoA for different combinations of bounds and Rz
Rz = 0.1; 
param.boundType='dataIndividual'; param.T_samples=360;
[data1,~,~,~,~,compTime1] = indirectDataDrivenControl_LMI(sys,param,Rz);
param.boundType='dataEllipsoidal'; param.T_samples=360;
[data2,~,~,~,~,compTime2] = indirectDataDrivenControl_LMI(sys,param,Rz); 
param.boundType='dataEllipsoidal'; param.T_samples=33; 
[data3,~,~,~,~,compTime3] = indirectDataDrivenControl_LMI(sys,param,Rz);
fprintf('Average computation time for Rz=%.1f: %.4f\n\n',Rz,(compTime1+compTime2+compTime3)/2)

Rz = 0.6; 
param.boundType='dataIndividual'; param.T_samples=2263;
[data4,~,~,~,~,compTime4] = indirectDataDrivenControl_LMI(sys,param,Rz);
param.boundType='dataEllipsoidal'; param.T_samples=2263;
[data5,~,~,~,~,compTime5] = indirectDataDrivenControl_LMI(sys,param,Rz);
param.boundType='dataEllipsoidal'; param.T_samples=213; 
[data6,~,~,~,~,compTime6] = indirectDataDrivenControl_LMI(sys,param,Rz);
fprintf('Average computation time for Rz=%.1f: %.4f\n\n',Rz,(compTime4+compTime5+compTime6)/2)

Rz = 0.9; 
param.boundType='dataIndividual'; param.T_samples=34668;
[data7,~,~,~,~,compTime7] = indirectDataDrivenControl_LMI(sys,param,Rz);
param.boundType='dataEllipsoidal'; param.T_samples=34668;
[data8,~,~,~,~,compTime8] = indirectDataDrivenControl_LMI(sys,param,Rz);
param.boundType='dataEllipsoidal'; param.T_samples=3999;
[data9,~,~,~,~,compTime9] = indirectDataDrivenControl_LMI(sys,param,Rz);
fprintf('Average computation time for Rz=%.1f: %.4f\n',Rz,(compTime7+compTime8+compTime9)/2)

%% Plotting
figure(1);clf;hold all;grid on;
xlabel('$x_1$', Interpreter='latex');
ylabel('$x_2$', Interpreter='latex');

% Rz = 0.10
p1=plot(data1(1,:),data1(2,:),'b-');
plot(data2(1,:),data2(2,:),'b.');
plot(data3(1,:),data3(2,:),'b--');

% Rz = 0.60
p2=plot(data4(1,:),data4(2,:),'r-');
plot(data5(1,:),data5(2,:),'r.');
plot(data6(1,:),data6(2,:),'r--');

% Rz = 0.90
p3=plot(data7(1,:),data7(2,:),'k-');
plot(data8(1,:),data8(2,:),'k.');
plot(data9(1,:),data9(2,:),'k--');

legend([p1,p2,p3],{'$c=0.1$','$c=0.6$','$c=0.9$'},Interpreter='latex')