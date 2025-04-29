function [Kn,Pinv,compTime] = controllerDesign_SOS(sys,eps,param,ud,degrees,z,cost)
%% Function to design a controller for a discrete-time bilinear system
% Inputs: 
%   - sys: Bilinear system description
%   - eps: small margin to ensure positive definiteness
%   - param: Parameters defining the area of interest and LS bounds
%   - ud: scalar polynomial denominator
%   - degrees: vector defining alpha and beta for the polynomial degrees
%   - z: substitute for the lifted state
%   - cost: optional, defined a objective for the optimization
%
% Outputs: 
%   - Kn: (n_u x n_x) polynomial matrix (numerator)
%   - Pinv: (n_x x n_x)-dimensional Lyapunov matrix 
%   - compTime: computation time of the controller design
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"
if nargin <= 6
    cost = '[]';
end
% Optimization: Decision variables
lmis = []; params = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(sys.n_x,sys.n_x,'symmetric');
    lmis = [lmis,P-eps.P*eye(size(P))>=0];
    % Controller matrix
    [Ln,Ln_coeff] = polymatrix(z,2*degrees.alpha-1,[sys.n_u,sys.n_x],'full');
    params = [params;vertcat(Ln_coeff{:})];

    % Decay factor rho
    rho = sdpvar(1);
    params = [params;rho];
    lmis = [lmis;rho>=eps.rho];

    % controller denominator needs to be SOS: check if ud is SOS
    [~,info_ud,~,Q_ud] = evalc("solvesos(sos(ud))");
    if info_ud.problem ~= 0 || min(eig(Q_ud{:})) < 0
        error("Chosen denominator is NOT SOS!")
    end
    
    if min(eig(param.QDelta)) < 0
        error('The matrix defining the proportional bound needs to be positive!')
    elseif param.QDelta == 0
        tau = 0;
        F12 = []; F22 = []; F23 = [];
    else
        % multiplier residual
        [tau,tau_coeff] = polynomial(z,2*degrees.beta); % needs to be even degree
        params = [params;tau_coeff];
        lmis = [lmis,tau_coeff(1)>=eps.tau];

        F12 = zeros(sys.n_x,sys.n_x + sys.n_u);
        F22 = tau*inv(param.QDelta);
        F23 = [P*ud;Ln];
    end
    F11 = ud*P - rho*ud*eye(sys.n_x) - tau*eye(sys.n_x);
    F13 = ud*sys.A*P + sys.B0*Ln + sys.Aux*kron(Ln,z);
    F33 = ud*P;
    
    F = [F11 ,F12 ,F13;
         F12',F22 ,F23;
         F13',F23',F33];
    
    [~,cost] = evalc(cost);
    lmis = [lmis,sos(F)];
    [~,opt] = evalc("solvesos(lmis,cost,sdpsettings('solver','mosek','sos.clean',1e-12),params)");

    compTime = opt.solvertime;
    switch opt.problem 
        case -1
            error('Controller design was not successful for T_samples=%i: %s: ILL_POSED',param.T_samples,opt.info)
        case 0 
            fprintf('Controller design completed for T_samples=%i: %s: PRIMAL_AND_DUAL_FEASIBLE\n',param.T_samples,opt.info)
            % Store obtained decision variables
            Pinv = double(P) \ eye(sys.n_x);
            Ln = replace(Ln,vertcat(Ln_coeff{:}),value(vertcat(Ln_coeff{:})));
            Kn = Ln / double(P);
        case 1
            error('Controller design was not successful for T_samples=%i: %s: DUAL_INFEASIBLE\n',param.T_samples,opt.info)
        case 2
            error('Controller design was not successful for T_samples=%i: %s: PRIMAL_INFEASIBLE\n',param.T_samples,opt.info)
        case 4
            error('Controller design was not successful for T_samples=%i: %s: UNKNOWN\n',param.T_samples,opt.info)
        otherwise
            error('Controller design was not successful for T_samples=%i: %s\n',param.T_samples,opt.info)
    end
end