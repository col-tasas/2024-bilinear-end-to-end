%% Matlab file to plot the comparison of identification errror and bounds

% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "24/10/01"

clear all 
close all 
clc 

% Load data
delta = 0.05;
name = "dataErrornx25Tmax250000nMC100delta5";
data = load("../data/" + name + ".mat");

% Set MATLAB figure style
set(groot, 'defaultAxesColorOrder', get(gca,'colororder')); % Default color order
set(groot,'defaultAxesFontSize', 14); % Set font size
set(groot,'defaultLineLineWidth', 1.5); % Set line width

data.errorData = data.errorData';
[nMC, n_evals] = size(data.errorBoundData);

normalize = 0;
ylog = 1;
conf = 3;


t_eval = data.t_eval;
t_sqrt = sqrt(t_eval);

error_m = mean(data.errorData, 1);

errorDataNorm = zeros(size(data.errorData));
errorBoundDataNorm = zeros(size(data.errorBoundData));

errorBoundNorm = data.errorBound .* sqrt(t_eval);
for j = 1:nMC
    errorDataNorm(j, :) = data.errorData(j,:) .* sqrt(t_eval);
    errorBoundDataNorm(j, :) = data.errorBoundData(j, :) .* sqrt(t_eval);
end

fig = figure(1);

for i = 1:2
    if i == 1
        t_eval = data.t_eval;
        errorData = data.errorData;
        errorBound = data.errorBound;
        errorBoundData = data.errorBoundData;
    elseif i == 2
        t_eval = data.t_eval;
        errorData = errorDataNorm;
        errorBound = errorBoundNorm;
        errorBoundData = errorBoundDataNorm;
    end

    subplot(2,1, i)
    ax = gca;
    hold(ax, "on")

    % Plot error bound
    plot(ax, t_eval, errorBound, 'DisplayName', 'A priori bound', Color='black');


    % Compute mean and std of error bound data
    errorBoundData_m = mean(errorBoundData, 1);
    errorBoundData_std = std(errorBoundData, 0, 1);

    % Plot Sigma band of error in Monte Carlo
    t_conf = [t_eval t_eval(end:-1:1)];
    errorBoundData_conf = [errorBoundData_m - conf * errorBoundData_std, errorBoundData_m(end:-1:1) + conf * errorBoundData_std(end:-1:1)];
    p1 = fill(ax, t_conf, errorBoundData_conf, 'red', HandleVisibility='off');
    p1.FaceColor = [1, 0.8, 0.8];
    p1.EdgeColor = 'none';


    % Plot mean error bound data
    plot(ax, t_eval, errorBoundData_m, 'DisplayName', 'Data-based bound (mean)', Color='red');

    % Compute the mean and std of the error
    error_m = mean(errorData, 1);
    error_std = std(errorData, 0, 1);

    % Plot Sigma band of error in Monte Carlo
    t_conf = [t_eval t_eval(end:-1:1)];
    error_conf = [error_m - conf * error_std, error_m(end:-1:1) + conf * error_std(end:-1:1)];
    p1 = fill(ax, t_conf, error_conf, 'blue', HandleVisibility='off');
    p1.FaceColor = [0.8, 0.8, 1];
    p1.EdgeColor = 'none';

    % Plot mean error
    plot(ax, t_eval, error_m, 'DisplayName', 'Monte Carlo Sim (mean)', Color='blue');

    hold(ax, 'off');

    % title(ax, sprintf('Error in B_1 (Bound with confidence \\delta = %.2f)', delta), 'Interpreter', 'latex');

    if ylog
        set(ax, 'YScale', 'log');
    end

    % set(ax, 'XScale', 'log');
    

    xlabel(ax, '$T$', Interpreter='latex');
    grid on 
    if i ==1
        ylabel(ax, '$\Vert \hat B_1- B_1\Vert_2$', Interpreter='latex');
        legend(ax, 'show', 'Location', 'northeast');
    elseif i==2
         ylabel(ax, '$\sqrt{T}\Vert \hat B_1- B_1\Vert_2$', Interpreter='latex');
        
    end
end 


