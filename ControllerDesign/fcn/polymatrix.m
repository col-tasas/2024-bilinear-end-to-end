function [Q_poly,Q_coeff] = polymatrix(x,d,dim,flag)    
%% Function to construct a polynominal matrix based on the Yalmip toolbox
% Inputs: 
%   - x: variable in which the matrix is polynomial
%   - d: polynomial degree
%   - dim: dimension of the polynomial matrix
%   - flag: flag if 'full' or 'symmetric' matrix
%
% Outputs: 
%   - Q_poly: dim-dimensional polynomial matrix
%   - Q_coeff: dim-dimensional polynomial matrix

%   - compTime: computation time of the controller design
%
% __author__ = "Robin Straesser"
% __contact__ = "robin.straesser@ist.uni-stuttgart.de"
% __date__ = "2025/04/24"
    if dim(1) ~= dim(2) && isequal(flag,'symmetric')
        error('Symmetric matrix must be square!')
    end
    Q_poly = sdpvar(dim(1),dim(2));
    Q_coeff = cell(dim(1),dim(2));
    for j = 1:dim(2)
        switch flag 
            case 'symmetric'
                for i = 1:j
                    [Q_poly(i,j),Q_coeff{i,j}] = polynomial(x,d);
                    if i ~= j
                        Q_poly(j,i) = Q_poly(i,j);
                        Q_coeff{j,i} = Q_coeff{i,j};
                    end
                end
            case 'full'
                for i = 1:dim(1)
                    [Q_poly(i,j), Q_coeff{i,j}] = polynomial(x,d); %Polynomiale Matrix
                end
        end
    end
end