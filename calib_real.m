clearvars
close all
addpath('common')

[settings] = readConfig('config/config_real.json');

currentDir = pwd; % Store current directory
%%
data=load('data/magcalMains.mat');
allan_dev = load("Allan/allan_dev.mat");
real_datasets = struct('settings', [], 'data', []);
settings_arr = repmat(settings, 1, 30);
for i = 1 : 30
    settings.magNoiseDensity = allan_dev.allan.mag(i);
    settings_arr(i) = settings;
    real_datasets(i).settings = settings;
    real_datasets(i).data.y_acc = data.y_acc';
    real_datasets(i).data.y_gyro = data.y_gyro';
    real_datasets(i).data.y_mag = data.y_mag(:, (i-1)*3+1:i*3)';
end
%%
poolobj = gcp; 


%% JointMAP

cd("JointMAP");
disp("==================Proposed===================")
tic
calibres_JointMAP = arrayfun(@calibrate_real, real_datasets, settings_arr);
time_proposed = toc;

cd(currentDir); % Return to current directory

%% EKF
cd("EKF");
disp("==================Wu et al.========================")
tic
calibres_EKF = arrayfun(@calibrate_real,real_datasets, settings_arr);
time_wu = toc;
disp("==============================================")
cd(currentDir); % Return to current directory


%% ML

cd("ML");
disp("==================Kok et al.========================")
tic
calibres_ML = arrayfun(@calibrate_real, real_datasets, settings_arr);
time_kok=toc;
disp("==============================================")
cd(currentDir);

%% save results to file
% if Folder results does not exist create it
if ~exist('results', 'dir')
    mkdir('results');
end

theta_all = zeros(12, 30);
for i = 1 : 30
    theta_all(1:9, i) = calibres_JointMAP(i).D(:);
    theta_all(10:12, i) = calibres_JointMAP(i).mag_bias;
end
save("results/calibJointMAP", 'theta_all');

for i = 1 : 30
    theta_all(1:9, i) = calibres_EKF(i).D(:);
    theta_all(10:12, i) = calibres_EKF(i).mag_bias;
end
save("results/calibEKF", 'theta_all');

for i = 1 : 30
    theta_all(1:9, i) = calibres_ML(i).D(:);
    theta_all(10:12, i) = calibres_ML(i).mag_bias;
end
save("results/calibML", 'theta_all');




%% plot results
% Define colors for consistency
bar_colors = [0 0.4470 0.7410;       % Wu et al. (Blue)
              0.8500 0.3250 0.0980;  % Kok et al. (Orange)
              0.4660 0.6740 0.1880]; % Proposed (Green)


figure;
hold on;

xpos = [0.6, 2, 3.4];

bar(xpos(1), time_kok, 'FaceColor', bar_colors(2, :)); % Assigning individual colors
bar(xpos(2), time_proposed, 'FaceColor', bar_colors(3, :));
bar(xpos(3), time_wu, 'FaceColor', bar_colors(1, :));

hold off;
set(gca, 'YScale', 'log'); % Log scale due to large differences
xticks(1:3);
xticks(xpos);
xticklabels({'Kok et al. (20 cores)','Proposed', 'Wu et al.'});

ylabel('Computation Time (s)');
%title('Computation Time Comparison');
ylim([1, max([time_kok]) * 1.1]);
grid on;
set(gca, 'Box', 'on', 'LineWidth', 1); % 'Box' 'on' creates solid border around the plot
ax=gca;

if ~exist('computation_results', 'dir')
    mkdir('computation_results');
end
save('computation_results/processing_time', 'time_kok', 'time_proposed', 'time_wu');