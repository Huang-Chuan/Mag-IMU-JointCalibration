clearvars
close all
addpath('common')
%% Load data and settings
[settings] = readConfig('config/config.json');
currentDir = pwd; % Store current directory
%% Generate simulation data
cd('simulate')
simulation_datasets = simulate(settings);
num_datasets = length(simulation_datasets);
cd(currentDir); % Return to current directory
%%
poolobj = gcp; 
%%
downsample_rates = [1, 2, 4, 8];

% Initialize arrays to store computation times
ekf_times = zeros(size(downsample_rates));
jointmap_times = zeros(size(downsample_rates));
ML_times = zeros(size(downsample_rates));


%%
% Loop through each downsample rate
for i = 1:length(downsample_rates)
    % Update settings
    settings.downsample_rate = downsample_rates(i);
    [low_pass_simulation_datasets] = lowpass_data(simulation_datasets, settings);
    % EKF computation
    cd("EKF");
    disp(['==================EKF (Downsample Rate = ', num2str(settings.downsample_rate), ')========================']);
    tic;
    calibres_EKF = arrayfun(@calibrate_real, low_pass_simulation_datasets, repmat(settings, num_datasets, 1));
    ekf_times(i) = toc;
    disp(['Time: ', num2str(ekf_times(i)), ' seconds']);
    cd(currentDir);
    statEKF(i) = compute_stat(calibres_EKF, [low_pass_simulation_datasets.calib_paras]);
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
    
    % JointMAP computation
    cd("JointMAP");
    disp(['==================Joint MAP (Downsample Rate = ', num2str(settings.downsample_rate), ')===================']);
    tic;
    calibres_JointMAP = arrayfun(@calibrate_real, low_pass_simulation_datasets, repmat(settings, num_datasets, 1));
    jointmap_times(i) = toc;
    disp(['Time: ', num2str(jointmap_times(i)), ' seconds']);
    cd(currentDir);
    statMAP(i) = compute_stat(calibres_JointMAP, [low_pass_simulation_datasets.calib_paras]);
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");


    % ML computation
    cd("ML");
    disp(['==================ML (Downsample Rate = ', num2str(settings.downsample_rate), ')===================']);
    tic;
    calibres_ML = arrayfun(@calibrate_real, low_pass_simulation_datasets, repmat(settings, num_datasets, 1));
    ML_times(i) = toc;
    disp(['Time: ', num2str(ML_times(i)), ' seconds']);
    cd(currentDir);
    statML(i) = compute_stat(calibres_ML, [low_pass_simulation_datasets.calib_paras]);
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

end

%% save results
if ~exist('computation_results', 'dir')
    mkdir('computation_results');
end
save('computation_results/results_ratio.mat', "statML", "statEKF", "statMAP", "ekf_times", "ML_times", "jointmap_times","downsample_rates");
visualize_stat(ekf_times, statEKF, jointmap_times, statMAP, ML_times, statML, downsample_rates, 'ratio');
