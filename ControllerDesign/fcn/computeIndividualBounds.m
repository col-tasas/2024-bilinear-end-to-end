function [eps] = computeIndividualBounds(sys, X, Y, T_samples, delta, stdNoise)
%% Function to compute LS estimates and bounds
% Inputs: 
%   - sys: Bilinear system description
%   - X: data matrix (states) for xp = Ax + w
%   - Y: data matrix (states and input 1) for xp = B_i x + [B_0]_i + w
%   - T_samples: (sys.n_u+1 x 1) number of samples generated for each of the
%   linear simulations
%   - delta: desired failure probability
%   - stdNoise: standard deviation of the noise (stdNoise^2 * I)
%
% Outputs: 
%   - eps: Data dependent errors 

% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"


%% Check burn-in time 
T_burnin = zeros(sys.n_u+1, 1);
T_burnin(1,1) = 0.5 * log( 2 * sys.n_u * 9^(2*sys.n_x)/delta); % 64 * (3 + 2 * sqrt(2)) * log(8* sys.n_u * 9^(sys.n_x) / delta);
for i = 1:1:sys.n_u 
    T_burnin(i+1,1) = 0.5 * log( 2 * sys.n_u * 9^(2*sys.n_x)/delta); % 64 * (3 + 2 * sqrt(2)) * log(8* sys.n_u * 9^(sys.n_x) / delta);
end 
 
if ~all(T_samples > T_burnin)
        error("Burn-in time not satisfied")
end

%% Individual bounds 
eps = zeros(sys.n_u + 1, 1);

% Identification of A
lambda_min = min(eig(X' * X));

eps(1) = 4 * stdNoise * (sqrt(2 * T_samples(1) * log(2 * sys.n_u * 9 ^ (sys.n_x) / delta)))/lambda_min; 

% Identification of B_i, [B_0]_i
for i=1:1:sys.n_u
    lambda_min = min(eig(Y{i}(:,:)' * Y{i}(:,:)));
    eps(i+1) = 4 /3 * sqrt(10) * stdNoise * (sqrt(2 * T_samples(i+1) * log(2 * sys.n_u * 9 ^ (sys.n_x) / delta)))/lambda_min;
end 


end 