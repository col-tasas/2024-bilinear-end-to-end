%% Script for example 2: Bilinear CSTR cooling
% Inputs: 
%   - none 
%
% Outputs: 
%   - none
%
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"

clear;clc;clearvars;

addpath("fcn\")
format long;

%% System dynamics
sys.A = [1.425, 0.1; -0.625, 0.8];
sys.B0 = [-0.025;0];
sys.Aux = [-0.1,0; 0,0];
sys.dynamics = @(x,u) sys.A*x + sys.B0*u + sys.Aux*kron(u,x);
%
[sys.n_x,sys.n_u] = size(sys.B0); % state and input dimension

%% data properties
param.sigma_x = [1;1]; 
param.stdNoise = 0.1; 
param.delta = 0.05;
param.umax = 2; param.umin = -2;
param.Qz = 'optimize';

%% Compute RoA for different combinations of bounds and Rz
Rz = 1e-4; 
param.boundType='dataIndividual'; param.T_samples=8430;
estimateRoA(Rz,sys,param);
param.boundType='dataEllipsoidal'; param.T_samples=1650;
estimateRoA(Rz,sys,param);

Rz = 1e-3; 
param.boundType='dataIndividual'; param.T_samples=8767; 
estimateRoA(Rz,sys,param);
param.boundType='dataEllipsoidal'; param.T_samples=1797;
estimateRoA(Rz,sys,param);

Rz = 1e-2; 
param.boundType='dataIndividual'; param.T_samples=11403; 
estimateRoA(Rz,sys,param);
param.boundType='dataEllipsoidal'; param.T_samples=2497;  
estimateRoA(Rz,sys,param);