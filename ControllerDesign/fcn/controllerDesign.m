function [K,Kw,Pinv,sys,compTime] = controllerDesign(sys,eps,param)
%% Function to design a (robust) controller for bilinear systems
% Inputs: 
%   - sys: Bilinear system description
%   - eps: small margin to ensure positive definiteness
%   - param: Parameters defining the area of interest and LS bounds
%
% Outputs: 
%   - K: (n_u x n_u) controller gain (numerator)
%   - Kw: (n_x x n_x) controller gain (denominator)
%   - Pinv: (n_x x n_x)-dimensional Lyapunov matrix 
%   - sys: Bilinear system description
%   - compTime: computation time of the controller design
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"

fprintf('Design a stabilizing controller...')
% Optimization: Decision variables
lmis = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(sys.n_x,sys.n_x,'symmetric');
    lmis = [lmis,P >= eps.P*eye(sys.n_x)];
    % State-feedback gain
    L = sdpvar(sys.n_u,sys.n_x,'full');
    % Full-information feedback gain
    Lw = sdpvar(sys.n_u,sys.n_u*sys.n_x,'full');
%     Lw = zeros(sys.n_u,sys.n_u*sys.n_x);

    % multiplier bilinearity
    Piinv = inv([sys.Pi.Qz,sys.Pi.Sz;sys.Pi.Sz',sys.Pi.Rz]);
    tQz = Piinv(1:sys.n_x,1:sys.n_x);
    tSz = Piinv(1:sys.n_x,sys.n_x+1:end);
    tRz = Piinv(sys.n_x+1:end,sys.n_x+1:end);
    Lambda = sdpvar(sys.n_u,sys.n_u,'symmetric');
    lmis = [lmis,Lambda - eps.Lambda*eye(sys.n_u) >= 0];
    
    if min(eig(param.QDelta)) < 0
        error('The matrix defining the proportional bound needs to be positive!')
    elseif param.QDelta == 0
        tau = 0;
        F13 = []; F23 = []; F33 = []; F34 = []; F35 = [];
    else
        % multiplier residual
        tau = sdpvar(1);
        lmis = [lmis,tau >= eps.tau];
        
        F13 = zeros(sys.n_x,sys.n_x+sys.n_u);
        F23 = -kron(eye(sys.n_u),tSz')*[zeros(sys.n_x*sys.n_u,sys.n_x),Lw'];
        F33 = tau*inv(param.QDelta);
        F34 = [P;L];
        F35 = -[zeros(sys.n_x,sys.n_x*sys.n_u);Lw];
    end
    
    F11 = P - tau*eye(sys.n_x);
    F12 = -sys.Aux*kron(Lambda,tSz) - sys.B0*Lw*kron(eye(sys.n_u),tSz);
    F14 = sys.A*P + sys.B0*L;
    F15 = sys.Aux*kron(Lambda,eye(sys.n_x)) + sys.B0*Lw;
    F22 = kron(Lambda,tRz) - Lw*kron(eye(sys.n_u),tSz) - (Lw*kron(eye(sys.n_u),tSz))';
    F24 = L;
    F25 = Lw;
    F44 = P;
    F45 = zeros(sys.n_x,sys.n_x*sys.n_u);
    F55 = -kron(Lambda,inv(tQz));

    F = [F11 ,F12 ,F13 ,F14 ,F15 ;
         F12',F22 ,F23 ,F24 ,F25 ;
         F13',F23',F33 ,F34 ,F35 ;
         F14',F24',F34',F44 ,F45 ;
         F15',F25',F35',F45',F55];
    
    lmis = [lmis,F - eps.F*eye(size(F)) >= 0];
    
    % invariance of safe operating region
    nu = sdpvar(1);
    lmis = [lmis,nu >= eps.nu];
    FI11 = nu*tRz - 1;
    FI12 = -nu*tSz';
    FI22 = nu*tQz + P;
    FI = [ FI11 ,FI12 ;
           FI12',FI22];
    lmis = [lmis,FI <= 0];
    
    % Solve optimization
    cost = -trace(P);
    opt = optimize(lmis,cost,sdpsettings('solver','mosek'));

    compTime = opt.solvertime;
    switch opt.problem 
        case -1
            error('Controller design was not successful for Rz=%.2g and T_samples=%i: %s: ILL_POSED',sys.Pi.Rz,param.T_samples,opt.info)
        case 0 
            fprintf('Controller design completed. %s: PRIMAL_AND_DUAL_FEASIBLE\n', opt.info)
            % Store obtained decision variables
            P = double(P);
            Pinv = P \ eye(sys.n_x);
            Lambda = value(Lambda);
        
            K = double(L) / P;
            Kw = double(Lw)*kron(inv(Lambda),eye(sys.n_x));
        case 1
            error('Controller design was not successful for Rz=%.2g and T_samples=%i: %s: DUAL_INFEASIBLE\n',sys.Pi.Rz,param.T_samples,opt.info)
        case 2
            error('Controller design was not successful for Rz=%.2g and T_samples=%i: %s: PRIMAL_INFEASIBLE\n',sys.Pi.Rz,param.T_samples,opt.info)
        case 4
            error('Controller design was not successful for Rz=%.2g and T_samples=%i: %s: UNKNOWN\n',sys.Pi.Rz,param.T_samples,opt.info)
        otherwise
            error('Controller design was not successful for Rz=%.2g and T_samples=%i: %s\n',sys.Pi.Rz,param.T_samples,opt.info)
    end
end