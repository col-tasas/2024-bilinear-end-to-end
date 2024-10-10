function [epsEll] = computeEllipsoidalBounds(sys, X, Y, T_samples, delta, stdNoise)
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
%   - epsEll: Ellipsoidal Data dependent errors (Complete RHS)
%
% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"


%% Check burn-in time
% Ellipsoidal
T_burnin_ell = zeros(sys.n_u+1, 1);

T_burnin_ell(1,1) = sys.n_x;            % 64 * (3 + 2 * sqrt(2)) * log(8* sys.n_u * 9^(sys.n_x) / delta);
for i = 1:1:sys.n_u
    T_burnin_ell(i+1,1) = sys.n_x + 1;  % 64 * (3 + 2 * sqrt(2)) * log(8* sys.n_u * 9^(sys.n_x) / delta);
end

if ~all(T_samples > T_burnin_ell)
    error("Burn-in time not satisfied")
end

%% Elipsoidal bounds
epsEll = cell(sys.n_u+1, 1);

% Identification of A
M = (X' * X);
epsEll{1} = stdNoise^2 * (2 * sqrt(sys.n_x) + sqrt(2 * log(1/delta)))^2 * (eye(size(M)) / M);

% Identification of B_i, [B_0]_i
for i=1:1:sys.n_u
    M = (Y{i}(:,:)' * Y{i}(:,:));
    epsEll{i+1} = stdNoise^2 * (sqrt(sys.n_x+1)+ sqrt(sys.n_x) + sqrt(2 * log(1/delta)))^2 * (eye(size(M)) / M);
end

end