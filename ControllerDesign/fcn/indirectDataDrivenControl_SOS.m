function [ufunc,Pinv,sys_hat,compTime] = indirectDataDrivenControl_SOS(sys,param,degrees)
%% Function to compute the region of attraction for a given Rz, sample size and bound type
% Inputs: 
%   - sys: Bilinear system description
%   - param: Parameters defining the area of interest and LS bounds
%   - degrees: degrees for SOS controller design
%
% Outputs:
%   - ufunc: SOS-based control law
%   - Pinv: (n_x x n_x)-dimensional Lyapunov matrix 
%   - sys_hat: Bilinear system description of identified system
%   - compTime: computation time of the controller design
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"

    T_samples=[param.T_samples;param.T_samples];     
    
    [X, Xp, Y, Yp, sys] = generateData(sys, T_samples, param.sigma_x, param.stdNoise);
    epsEll = computeEllipsoidalBounds(sys, X, Y, T_samples, param.delta, param.stdNoise);
    epsIndiv = computeIndividualBounds(sys, X, Y, T_samples, param.delta, param.stdNoise);
    sys_hat = computeOLS(sys, X, Xp, Y, Yp);
    
    param.epsA = epsIndiv(1);
    param.epsB0 = epsIndiv(2:end);
    param.epsB = epsIndiv(2:end)./sqrt(min(param.sigma_x));
    %
    param.EpsA = epsEll{1};
    param.EpsB = cell(sys.n_u,1);
    for j = 1:sys.n_u
        param.EpsB{j,1} = epsEll{j+1};
    end
    
    % matrix characterizing the residual bound
    param = computeProportionalBounds(param,sys);
    if isfield(param,'overapproxQDelta') && param.overapproxQDelta
        param.QDelta = max(eig(param.QDelta))*eye(sys.n_x+sys.n_u);
    end
    
    %% Controller design
    eps.P = 1e-6;
    eps.F = 1e-6;
    eps.Lambda = 1e-7;
    eps.tau = 1e-7;
    eps.nu = 1e-7;
    eps.rho = 0;
    %
    x = sdpvar(sys.n_x,1);
    [~,c,v] = polynomial(x,2*degrees.alpha,2*degrees.alpha);
    ud = 1 + ones(size(c))'*v;
    %
    [Kn,Pinv,compTime] = controllerDesign_SOS(sys_hat,eps,param,ud,degrees,x,'[]');

    Knfunc = @(xvar) replace(Kn,x,xvar);
    udfunc = @(xvar) replace(ud,x,xvar);
    ufunc = @(x) 1/(udfunc(x))*Knfunc(x)*x;
    
    %% Computation times
    fprintf('Computation time of controller design: %fs\n',compTime)
end