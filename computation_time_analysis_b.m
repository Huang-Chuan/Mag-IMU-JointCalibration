clearvars
close all
addpath('common')
%% Load data and settings
[settings] = readConfig('config/config.json');
currentDir = pwd; % Store current directory

sample_freq = [10,20,40,80];
%sample_freq = [10];
% Initialize arrays to store computation times
ekf_times = zeros(size(sample_freq));
jointmap_times = zeros(size(sample_freq));
ML_times = zeros(size(sample_freq));


%%
poolobj = gcp; % If no pool, do not create new one.

%% Loop through each downsample rate
for i = 1:length(sample_freq)
    % Update settings
    settings.downsample_rate = 1;
    settings.freq = sample_freq(i);
    settings.dt = 1/settings.freq;
    % Generate simulation data
    cd('simulate')
    simulation_datasets = simulate(settings);
    num_datasets = length(simulation_datasets);
    cd(currentDir); % Return to current directory

    % EKF computation
    cd("EKF");
    disp(['==================EKF (Sample frequency = ', num2str(settings.freq), ')========================']);
    tic;
    calibres_EKF = arrayfun(@calibrate_real, simulation_datasets, repmat(settings, num_datasets, 1));
    ekf_times(i) = toc;
    disp(['Time: ', num2str(ekf_times(i)), ' seconds']);
    cd(currentDir);
    statEKF(i) = compute_stat(calibres_EKF, [simulation_datasets.calib_paras]);
    %compute_ANEES(calibres_EKF, [simulation_datasets.calib_paras])
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
    
    %JointMAP computation
    cd("JointMAP");
    disp(['==================Joint MAP (Sample frequency = ', num2str(settings.freq), ')===================']);
    tic;
    calibres_JointMAP = arrayfun(@calibrate_real, simulation_datasets, repmat(settings, num_datasets, 1));
    jointmap_times(i) = toc;
    disp(['Time: ', num2str(jointmap_times(i)), ' seconds']);
    cd(currentDir);
    statMAP(i) = compute_stat(calibres_JointMAP, [simulation_datasets.calib_paras]);
    %compute_ANEES(calibres_JointMAP, [simulation_datasets.calib_paras])
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
    % 
    % 
    % % Joint ML computation
    cd("ML");
    disp(['==================ML (Sample frequency = ', num2str(settings.freq), ')===================']);
    tic;
    calibres_ML = arrayfun(@calibrate_real, simulation_datasets, repmat(settings, num_datasets, 1));
    ML_times(i) = toc;
    disp(['Time: ', num2str(ML_times(i)), ' seconds']);
    cd(currentDir);
    statML(i) = compute_stat(calibres_ML, [simulation_datasets.calib_paras]);
    disp("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

end
%%
if ~exist('computation_results', 'dir')
    mkdir('computation_results');
end
save('computation_results/results_freq.mat', "statML", "statEKF", "statMAP", "ekf_times", "ML_times", "jointmap_times","sample_freq");
visualize_stat(ekf_times, statEKF, jointmap_times, statMAP, ML_times, statML, sample_freq, 'freq');
