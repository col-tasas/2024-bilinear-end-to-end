%% Script for example 1: Bilinear academic example -- SOS approach
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
%
degrees.alpha = 1;
degrees.beta = degrees.alpha;

%% Compute controller for different bounds
param.boundType='dataIndividual'; param.T_samples=5145;
[u1,~,~,compTime1] = indirectDataDrivenControl_SOS(sys,param,degrees);
param.boundType='dataEllipsoidal'; param.T_samples=959;
[u2,~,~,compTime2] = indirectDataDrivenControl_SOS(sys,param,degrees); 
param.boundType='dataEllipsoidal'; param.T_samples=5145; 
[u3,~,~,compTime3] = indirectDataDrivenControl_SOS(sys,param,degrees);
fprintf('Average computation time: %.4f\n\n',(compTime1+compTime2+compTime3)/2)