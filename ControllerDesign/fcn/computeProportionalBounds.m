function [param] = computeProportionalBounds(param,sys)
%% Function to compute proportional bounds for control based on LS bounds
% Inputs: 
%   - param: Parameters defining the area of interest and LS bounds
%   - sys: Bilinear system description
%
% Outputs: 
%   - param: Updated parameters to contain proportional bound
%   characterization
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2024/10/01"

if isfield(param,'boundType')
    switch param.boundType
        case {'sampleComplexity','dataIndividual'}
            cx = max(abs(1-sys.n_u*param.umin),abs(1-sys.n_u*param.umax))*param.epsA ...
                + sum(max(abs(param.umin),abs(param.umax))*param.epsB);
            cu = sqrt(sum(param.epsB0.^2));
            param.QDelta = blkdiag(cx^2*eye(sys.n_x),cu^2*eye(sys.n_u));
        case 'dataEllipsoidal'
            tildeEpsB11 = [];
            tildeEpsB12 = [];
            tildeEpsB21 = [];
            tildeEpsB22 = [];
            for j=1:sys.n_u
                tildeEpsB11 = blkdiag(tildeEpsB11,param.EpsB{j,1}(1:end-1,1:end-1));
                tildeEpsB12 = blkdiag(tildeEpsB12,param.EpsB{j,1}(1:end-1,end));
                tildeEpsB21 = blkdiag(tildeEpsB21,param.EpsB{j,1}(end,1:end-1));
                tildeEpsB22 = blkdiag(tildeEpsB22,param.EpsB{j,1}(end,end));
            end
            tildeEpsB = [tildeEpsB11,tildeEpsB12;tildeEpsB21,tildeEpsB22];
            hatEpsB ... 
                = blkdiag(kron(max(abs(param.umin),abs(param.umax))*ones(sys.n_u,1),eye(sys.n_x)),eye(sys.n_u))' *...
                  tildeEpsB *...  
                  blkdiag(kron(max(abs(param.umin),abs(param.umax))*ones(sys.n_u,1),eye(sys.n_x)),eye(sys.n_u));
            param.QDelta ...
                = blkdiag(...
                    (sys.n_u+1)*max(abs(1-sys.n_u*param.umin),abs(1-sys.n_u*param.umax))^2*param.EpsA,...
                    zeros(sys.n_u)...
                ) + (sys.n_u+1)*hatEpsB;
        otherwise
            param.QDelta = 0; 
    end
else
    param.QDelta = 0;
end
end