function sys_hat = computeOLS(sys, X, Xp, Y, Yp)
%% Function to compute OLS from data
% Inputs:
%   - sys: system
%   - X: data matrix (states) for xp = Ax + w
%   - Xp: data matrix (xp) for xp = Ax + w
%   - Y: data matrix (states and input 1) for xp = B_i x + [B_0]_i + w
%   - Yp: data matrix (xp and input 1) for xp = B_i x + [B_0]_i + w
%
% Output:
%   - sys_hat: Estimate
%
% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"


%% Estimates (using CL solution)
A_hat = ((X' * X) \ (X' * Xp))';
Ai_hat = zeros(sys.n_x, sys.n_x*sys.n_u);
B0_hat = zeros(sys.n_x, sys.n_u);
for i = 1:1:sys.n_u
    theta_hat = ((Y{i}(:,:)' * Y{i}(:,:)) \ (Y{i}(:,:)' * Yp{i}(:,:)))';
    Ai_hat(:, (i-1)*sys.n_x+1 : i*sys.n_x)  = theta_hat(:, 1:sys.n_x) - A_hat;
    B0_hat(:, i) = theta_hat(:, sys.n_x+1:end);
end

sys_hat = sys;
sys_hat.A = A_hat;
sys_hat.B0 = B0_hat;
sys_hat.Aux = Ai_hat;

end