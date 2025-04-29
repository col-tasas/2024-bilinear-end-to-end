%% Script for example 3: Bilinear CSTR cooling -- SOS approach
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
sys.A = [1.425, 0.1; -0.625, 0.8];
sys.B0 = [-0.25;0];
sys.Aux = [-0.1,0; 0,0];
sys.dynamics = @(x,u) sys.A*x + sys.B0*u + sys.Aux*kron(u,x);
%
[sys.n_x,sys.n_u] = size(sys.B0); % state and input dimension

%% data properties
param.sigma_x = 1; 
param.stdNoise = 0.1; 
param.delta = 0.05;
param.umax = 2; param.umin = -2;
%
degrees.alpha = 4;
degrees.beta = degrees.alpha;

%% Compute controller for different bounds
param.boundType='dataIndividual'; param.T_samples=7234;
[u1,~,sys_hat1] = indirectDataDrivenControl_SOS(sys,param,degrees);
param.boundType='dataEllipsoidal'; param.T_samples=969;
[u2,~,sys_hat2] = indirectDataDrivenControl_SOS(sys,param,degrees);