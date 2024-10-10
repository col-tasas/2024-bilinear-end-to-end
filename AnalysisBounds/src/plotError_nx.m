%% Matlab file to plot the comparison of identification errror and bounds over problem dimension n_x

% __author__ = "Nicolas Chatzikiriakos"
% __contact__ = "nicolas.chatzikiriakos@ist.uni-stuttgart.de"
% __date__ = "24/10/01"

clear all
close all
clc


% Load data
delta = 0.05;
conf = 3;
names = ["dataErrornx1Tmax250000nMC100delta5", ... 
    "dataErrornx2Tmax250000nMC100delta5", ...
    "dataErrornx3Tmax250000nMC100delta5", ...
    "dataErrornx4Tmax250000nMC100delta5", ...
    "dataErrornx5Tmax250000nMC100delta5", ...
    "dataErrornx7Tmax250000nMC100delta5", ...
    "dataErrornx10Tmax250000nMC100delta5", ...
    "dataErrornx15Tmax250000nMC100delta5",...
    "dataErrornx20Tmax250000nMC100delta5", ...
    "dataErrornx25Tmax250000nMC100delta5", ...
    "dataErrornx30Tmax250000nMC100delta5"];

stateDims = [1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30];
L = length(names);

error_mean = zeros(1, L);
error_std = zeros(1, L);

errorBound = zeros(1, L);

errorBoundData_mean = zeros(1, L);
errorBoundData_std =  zeros(1, L);

fig = figure(1);


for j=1:1:L

    data = load("../data/" + names(j) + ".mat");
    data.errorData = data.errorData';

    error_mean(j) = mean(data.errorData(:, end));
    error_std(j) = std(data.errorData(:, end));

    errorBound(j) = data.errorBound(end);

    errorBoundData_mean(j) =  mean(data.errorBoundData(:, end));
    errorBoundData_std(j) =   std(data.errorBoundData(:, end));

end
subplot(1,2,1)
hold on 
ax = gca;
% Plot error bound
plot(ax, stateDims, errorBound, 'DisplayName', 'A priori bound', Color='black', Marker='X');


% Plot Sigma band of error in Monte Carlo
stateDims_conf = [stateDims stateDims(end:-1:1)];
errorBoundData_conf = [errorBoundData_mean - conf * errorBoundData_std, errorBoundData_mean(end:-1:1) + conf * errorBoundData_std(end:-1:1)];
p1 = fill(ax, stateDims_conf, errorBoundData_conf, 'red', HandleVisibility='off');
p1.FaceColor = [1, 0.8, 0.8];
p1.EdgeColor = 'none';


% Plot mean error bound data
plot(ax, stateDims, errorBoundData_mean, 'DisplayName', 'Data-based bound (mean)', Color='red', Marker='X');


% Plot Sigma band of error in Monte Carlo
stateDims_conf = [stateDims stateDims(end:-1:1)];
error_conf = [error_mean - conf * error_std, error_mean(end:-1:1) + conf * error_std(end:-1:1)];
p1 = fill(ax, stateDims_conf, error_conf, 'blue', HandleVisibility='off');
p1.FaceColor = [0.8, 0.8, 1];
p1.EdgeColor = 'none';

% Plot mean error
plot(ax, stateDims, error_mean, 'DisplayName', 'Monte Carlo Sim (mean)', Color='blue', Marker='x');


grid on
xlabel(ax, '$n_x$', Interpreter='latex');
ylabel(ax, '$\Vert \hat B_1- B_1\Vert_2$', Interpreter='latex');
legend(ax, 'show', 'Location', 'northeast');


% figure()
subplot(1,2,2)
ax2 = gca ;
plot(ax2, stateDims,errorBoundData_mean./error_mean, 'DisplayName', 'Data-Based $\varepsilon_{B_1}$', Marker="x" )
hold on
plot(ax2, stateDims, errorBound./error_mean, 'DisplayName', 'A priori $\varepsilon_{B_1}$', Marker='X')
% yline(ax2, 1)


grid on 
xlabel(ax2, '$n_x$', Interpreter='latex');
ylabel(ax2, '$\frac{\varepsilon_{B_1}}{\Vert \hat B_1- B_1\Vert_2}$', Interpreter='latex');
legend(ax2, 'show', 'Location', 'northeast', Interpreter='latex');
ylim(ax2, [ 2.5, 40])
xlim(ax2, [1, 30])

