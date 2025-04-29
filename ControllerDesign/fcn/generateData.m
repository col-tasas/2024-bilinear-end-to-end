function [X, Xp, Y, Yp, sys] = generateData(sys, T_samples, sigma_z, sigma_w)
%% Function for data generation 
% Inputs: 
%   - sys: Bilinear system description
%   - T_samples: (sys.n_u+1 x 1) number of samples generated for each of the
%   linear simulations
%   - sigma_z, sigma_w (sys.n_x x 1), variance of each state/noise entry 
%
% Outputs: 
%   - X: data matrix (states) for xp = Ax + w
%   - Xp: data matrix (xp) for xp = Ax + w 
%   - Y: data matrix (states and input 1) for xp = B_i x + [B_0]_i + w
%   - Yp: data matrix (xp and input 1) for xp = B_i x + [B_0]_i + w
%
%
% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"

rng('default')

sys.n_z = sys.n_x;          % original state space dimension 
% Check whether lifted bilinear system representation
if ~isfield(sys,'lifting')
    sys.lifting = @(x) x;
else 
    sys.n_x = size(sys.lifting(zeros(sys.n_z,1)),1); % lifted state space dimension 
end

% Init stuff
X = zeros(T_samples(1), sys.n_x);
Xp = zeros(T_samples(1), sys.n_x);

Y = cell(sys.n_u,1);
Yp = cell(sys.n_u,1);

for i = 1:1:sys.n_u
    Y{i} =zeros(T_samples(i), sys.n_x + 1);
    Yp{i} = zeros(T_samples(i), sys.n_x);
end


% Data for x^+ = Ax + w
u0 = zeros(sys.n_u,1);
for t = 1:1:T_samples(1)
    z = normrnd(zeros(sys.n_z, 1), sigma_z);
    z_p = sys.dynamics(z,u0);
    w = normrnd(zeros(sys.n_x, 1), sigma_w);
    x = sys.lifting(z);
    x_p = sys.lifting(z_p) + w;
    X(t, :) = x;
    Xp(t,:) = x_p;
end

% Data for x^+ = B_i x + [B_0]_i + w
for i = 1:1:sys.n_u
    ui=u0; ui(i,1)=1;
    for t = 1:1:T_samples(i+1)
        z = normrnd(zeros(sys.n_z, 1), sigma_z);
        z_p = sys.dynamics(z,ui);
        w = normrnd(zeros(sys.n_x, 1), sigma_w);
        x = sys.lifting(z);
        x_p = sys.lifting(z_p) + w;
        Y{i}(t, :) = [x; 1];
        Yp{i}(t, :) = x_p;
    end
end
end