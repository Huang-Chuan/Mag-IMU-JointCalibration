function [simulation_data] = simulate(settings)
    rng(1);
    N = settings.N; % the number of simulation datasets

    simulation_data = repmat(struct('calib_paras', [], 'settings', [], 'data', []), N, 1);
    for i = 1 : N
        % generate the calibration parameters
        simulation_data(i).calib_paras = generate_calib_paras(settings);
        % generate the settings
        simulation_data(i).settings = generate_settings(settings);
        % generate the data
        [simulation_data(i).data.y_acc, simulation_data(i).data.y_mag, simulation_data(i).data.y_gyro] = generate_data(simulation_data(i).calib_paras, simulation_data(i).settings);
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Generate calibration parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calib_paras = generate_calib_paras(settings)
    calib_paras.acc_bias = (2 * rand(3, 1) - 1) * 0.5;                      % accelerometer bias uniformly distributed ~[-0.5, 0.5] m/s^2
    calib_paras.mag_bias = (2 * rand(3, 1) - 1) * 2;                        % magnetometer bias uniformly distributed ~[-2, 2] uT
    calib_paras.D_scale =  diag([1, 1, 1] + (2 * rand(1, 3) - 1) * 0.1);    % scale matrix, diagonal matrix with diagonal elements uniformly distributed ~[0.9, 1.1]
    skewAngle = (2 * rand(1, 3) - 1) * deg2rad(10);                         % skew angle, uniformly distributed ~[-10, 10] degrees
    calib_paras.D_skew = [1, 0, 0;                                          % skew-symmetric matrix
                         sin(skewAngle(1)), cos(skewAngle(1)), 0;
                         -sin(skewAngle(2)), cos(skewAngle(2)) * sin(skewAngle(3)), cos(skewAngle(2)) * cos(skewAngle(3))];
    misalignAngle = (2 * rand(1, 3) - 1) * deg2rad(5);                     % misalignment angle, uniformly distributed ~[-5, 5] degrees
    calib_paras.D_rot = rotx(misalignAngle(1)) * roty(misalignAngle(2)) * rotz(misalignAngle(3)); % rotation matrix
    calib_paras.D = 27*calib_paras.D_scale * calib_paras.D_skew * calib_paras.D_rot; % the calibration matrix                    
    calib_paras.dip_angle = deg2rad(settings.dip_angle) + (2 * rand() - 1) * deg2rad(5);           % dip angle, uniformly distributed ~[69.5, 79.5] degrees
    calib_paras.gyro_bias = settings.nominal_gyro_bias + (2 * rand(3, 1) - 1) * deg2rad(0.1); % gyro bias, uniformly distributed ~[-0.1, 0.1] degree/s;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Generate settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulationSettings = generate_settings(settings)

%           accNoiseDensity: accelerometer noise density (m/s^2/sqrt(Hz))
%           gyroNoiseDensity: gyroscope noise density (rad/s/sqrt(Hz))
%           magNoiseDensity: magnetometer noise density (u T/sqrt(Hz))

    simulationSettings.g = -settings.g;
    simulationSettings.freq = settings.freq;
    simulationSettings.dt = 1/settings.freq;

    simulationSettings.accNoiseDensity = settings.accNoiseDensity;
    simulationSettings.magNoiseDensity = settings.magNoiseDensity ;
    simulationSettings.gyroNoiseDensity = settings.gyroNoiseDensity; 

    simulationSettings.sigma_acc = settings.accNoiseDensity * sqrt(settings.freq);
    simulationSettings.sigma_mag = settings.magNoiseDensity * sqrt(settings.freq);
    simulationSettings.sigma_gyro = settings.gyroNoiseDensity * sqrt(settings.freq);

    simulationSettings.rotationTrajectory = generate_rotation_trajectory(settings);
    simulationSettings.K = size(simulationSettings.rotationTrajectory, 3);
    simulationSettings.verbose = false;
end


function rotationTrajectory = generate_rotation_trajectory(settings)
    rotationTrajectory = zeros(3, 3, 1000);
    rotationTrajectory(:, :, 1) = eye(3);
    idx = 1;
    R = rotationTrajectory(:, :, 1);
    dt = settings.dt;
    v_mtx = 0.1*[1 0 0; 0 1 0; 0 0 1;1 1 0; 1 0 1; 0 1 1];
    for i = 1 : size(v_mtx, 1)
        v = v_mtx(i, :)';
        t = 0;
        while(t < (2 * pi)/norm(v) )
            R = R * SO3_exp((v + 0.05 * randn(3, 1)) * dt);
            t = t + dt; 
            idx = idx + 1;
            rotationTrajectory(:, :, idx) = R;
        end
    end
end


function Rx = rotx(theta)
    Rx = [1, 0, 0; 
            0, cos(theta), -sin(theta); 
            0, sin(theta), cos(theta)];
end

function Ry = roty(theta)
    Ry = [cos(theta), 0, sin(theta); 
            0, 1, 0; 
            -sin(theta), 0, cos(theta)];
end

function Rz = rotz(theta)
    Rz = [cos(theta), -sin(theta), 0; 
            sin(theta), cos(theta), 0; 
            0, 0, 1];
end