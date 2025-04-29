%% Script for example 2: Nonlinear inverted pendulum -- SOS approach
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
m = 1; l = 1; b = 0.5; g = 9.81; Ts=0.01;
sys.dynamics = @(z,u) [z(1,:) + Ts*z(2,:); 
                       z(2,:) + Ts*(g/l*sin(z(1,:)) - b/m/l^2*z(2,:) + 1/m/l^2*u)];
%
sys.n_x = 2; sys.n_u = 1;

%% Koopman lifting
sys.lifting = @(z) [z; sin(z(1))];

%% data properties
param.sigma_x = 1; 
param.stdNoise = 1e-3; 
param.delta = 0.05;
param.umax = 2; param.umin = -2;
%
degrees.alpha = 2;
degrees.beta = degrees.alpha;
%
param.overapproxQDelta = true;

%% Compute controller
param.boundType='dataEllipsoidal'; param.T_samples = 5000; degrees.alpha=3;
degrees.beta = degrees.alpha;
[u,Pinv,sys] = indirectDataDrivenControl_SOS(sys,param,degrees);

%% Plotting
fprintf('Plotting...')
figure(1);clf;hold all;grid on;
xlabel('$x_1$', Interpreter='latex');
ylabel('$x_2$', Interpreter='latex');

% Plot closed-loop trajectories
param.Tsim = 2/Ts;
x0s = [-pi  ,-pi  ,-pi  ,-pi  ,-2.00,-1.00, 0.00, 1.00, 2.00, 2.685, pi  , pi  , pi  , pi  , 2.00, 1.00, 0.00,-1.00,-2.00,-2.675;
        0.00, 1.00, 2.00, pi  , pi  , pi  , pi  , pi  , pi  , pi   , 0.00,-1.00,-2.00,-pi  ,-pi  ,-pi  ,-pi  ,-pi  ,-pi  ,-pi  ];
for i = 1:size(x0s,2)
    xsim = NaN(sys.n_z,param.Tsim+1);
    xsim(:,1) = x0s(1:sys.n_z,i);
    usim = NaN(sys.n_u,param.Tsim);
    for t = 1:param.Tsim
        usim(:,t) = u(sys.lifting(xsim(:,t)));
        xsim(:,t+1) = sys.dynamics(xsim(:,t),usim(:,t));
    end
    plot(xsim(1,:),xsim(2,:),'k')
    drawnow
end
axis equal
xlim([-pi,pi])
ylim([-pi,pi])
fprintf('Done.\n')