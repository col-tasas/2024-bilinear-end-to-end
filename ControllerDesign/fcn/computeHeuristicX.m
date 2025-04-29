function QzHeuristics = computeHeuristicX(sys,eps,param)
%% Function to compute an heuristically optimized region X as outer-approximation of the region of attraction
% Inputs: 
%   - sys: Bilinear system description
%   - eps: small margin to ensure positive definiteness
%   - param: Parameters defining the area of interest and LS bounds
%
% Outputs: 
%   - QzHeuristics: (n_x x n_x) upper left entry of matrix defining X
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"



% Optimization: Decision variables
lmis = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(sys.n_x,sys.n_x,'symmetric');
    lmis = [lmis,P >= eps.P*eye(sys.n_x)];
    % State-feedback gain
    L = sdpvar(sys.n_u,sys.n_x,'full');
    % Full-information feedback gain
    Lw = sdpvar(sys.n_u,sys.n_u*sys.n_x,'full');

    % multiplier bilinearity
    sys.Pi.Qz = -eye(sys.n_x);
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
        % multiplier Koopman approximation error
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
    
    % Solve optimization
    opt = optimize(lmis,[],sdpsettings('solver','mosek'));

    if opt.problem ~= 0
        warning('Pre-computing Qz not successful! Using Qz=-I instead.')
        QzHeuristics = -eye(sys.n_x);
    else
        Pinv = double(P) \ eye(sys.n_x);
        QzHeuristics = -Pinv/norm(Pinv);
    end
end