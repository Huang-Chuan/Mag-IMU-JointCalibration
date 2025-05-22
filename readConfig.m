% Copyright 2023, Chuan Huang
    
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at

%   http://www.apache.org/licenses/LICENSE-2.0

% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
function [settings] = readConfig(configFilePath)

    try
        configData = jsondecode(fileread(configFilePath));
    catch
        error('Error reading or parsing the JSON configuration file.');
    end
    
    settings.freq = configData.Simulation.frequency;
    settings.dt = 1 / settings.freq;
    settings.nominal_gyro_bias = configData.Simulation.nominal_gyro_bias;
    settings.accNoiseDensity = configData.Simulation.accNoiseDensity;
    settings.magNoiseDensity = configData.Simulation.magNoiseDensity;
    settings.gyroNoiseDensity = configData.Simulation.gyroNoiseDensity;
    settings.g = configData.Simulation.gravity;
    settings.N = configData.Simulation.num_datasets;
    settings.downsample_rate = configData.downsample_rate;
    settings.num_padding = configData.num_padding_samples;
    settings.compute_covariance = configData.compute_covariance;




    settings.minAcc = configData.minAcc;
    settings.maxAcc = configData.maxAcc;
    settings.dip_angle = configData.dip_angle;
    settings.intrisic_calib_samples = configData.intrisic_calib_samples;
    settings.alignment_interval = configData.alignment_interval;


    settings.ML.init_acc_bias = configData.ML.init_acc_bias;
    settings.ML.init_gyro_bias = configData.ML.init_gyro_bias;
    settings.ML.init_pose = configData.ML.init_pose;
    settings.ML.init_pose_uncertainty = configData.ML.init_pose_uncertainty;

    settings.EKF.m0_noise_std = configData.EKF.m0_noise_std;
    settings.EKF.g0_noise_std = configData.EKF.g0_noise_std;
    settings.EKF.init_gyro_bias = configData.EKF.init_gyro_bias;
    settings.EKF.init_gyro_bias_std = configData.EKF.init_gyro_bias_std;
    settings.EKF.init_acc_bias_std = configData.EKF.init_acc_bias_std;
    settings.EKF.init_S_std = configData.EKF.init_S_std;
    settings.EKF.init_h_std = configData.EKF.init_h_std;
    settings.EKF.init_m0_std = configData.EKF.init_m0_std;
    settings.EKF.init_g0_std = configData.EKF.init_g0_std;
    settings.EKF.gyro_bias_noise_std = configData.EKF.gyro_bias_noise_std;
    settings.EKF.acc_bias_noise_std = configData.EKF.acc_bias_noise_std;
    
    if(configData.EKF.ANIS == 0)
        settings.EKF.compute_ANIS = false;
    else
        settings.EKF.compute_ANIS = true;
    end




end
