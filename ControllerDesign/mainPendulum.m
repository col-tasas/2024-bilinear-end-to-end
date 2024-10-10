%% Script for example 3: Nonlinear inverted pendulum 
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
% Set MATLAB figure style
set(groot, 'defaultAxesColorOrder', get(gca,'colororder')); % Default color order
set(groot,'defaultAxesFontSize', 14); % Set font size
set(groot,'defaultLineLineWidth', 1.5); % Set line width

addpath("fcn\")
format long;

%% System dynamics
m = 1; l = 1; b = 0.5; g = 9.81; Ts=0.1;
sys.dynamics = @(z,u) [z(1,:) + Ts*z(2,:); 
                       z(2,:) + Ts*(g/l*sin(z(1,:)) - b/m/l^2*z(2,:) + 1/m/l^2*u)];
%
sys.n_x = 2; sys.n_u = 1;

%% Koopman lifting
sys.lifting = @(z) [z; sin(z(1))];

%% data properties
param.sigma_x = [1;1]; 
param.stdNoise = 1e-3; 
param.delta = 0.05;
param.umax = 2; param.umin = -2;
param.Qz = 'eye';
param.overapproxQDelta = true;

%% Compute RoA and controller gains
Rz = 11;  
param.boundType='dataEllipsoidal'; param.T_samples = 50000;
[~,K,Kw,Pinv,sys] = estimateRoA(Rz,sys,param);

%% Plotting
fprintf('Plotting...')
figure(1);clf;hold all;grid on;
xlabel('$x_1$', Interpreter='latex');
ylabel('$x_2$', Interpreter='latex');

% Sample data from RoA in original state space
x_dist = 0.01; y_dist = x_dist;
xplotmax = 3.3; xplotmin = -xplotmax;
yplotmax = 3.3; yplotmin = -yplotmax;
xx = xplotmin:x_dist:xplotmax;
yy = yplotmin:y_dist:yplotmax;
[XX,YY] = meshgrid(xx,yy);
xSamples = [horzcat(XX(:))';horzcat(YY(:))'];
xRoA = zeros(sys.n_z,1);
for i = 1:size(xSamples,2)
    if sys.lifting(xSamples(:,i))'*Pinv*sys.lifting(xSamples(:,i)) <= 1
        xRoA = [xRoA,xSamples(:,i)];
    end
end
xRoABoundary = boundary(xRoA(1,:)',xRoA(2,:)');

% Plot RoA in original state space
plot(xRoA(1,xRoABoundary),xRoA(2,xRoABoundary),'k')

% Plot closed-loop trajectories
u = @(z) (eye(sys.n_u)-Kw*kron(eye(sys.n_u),sys.lifting(z))) \ (K*sys.lifting(z));
param.Tsim = 100/Ts;
xRoAinit = [ 2.00,-2.00, 1.00,-1.00,-pi  , pi  ,-0.16, 0.16,-1.00, 1.00,-2.00, 2.00;
             1.50,-1.50, 2.20,-2.20, 0.00, 0.00, 3.05,-3.05, 2.60,-2.60, 2.05,-2.05];
for i = 1:size(xRoAinit,2)
    xsim = NaN(sys.n_z,param.Tsim+1);
    xsim(:,1) = xRoAinit(1:sys.n_z,i);
    usim = NaN(sys.n_u,param.Tsim);
    for t = 1:param.Tsim
        usim(:,t) = u(xsim(:,t));
        xsim(:,t+1) = sys.dynamics(xsim(:,t),usim(:,t));
    end
    plot(xsim(1,:),xsim(2,:),'r')
end
xlim([xplotmin,xplotmax])
ylim([yplotmin,yplotmax])
fprintf('Done.\n')
