function [calibres, covs, innovations, Ss, ANIS] = jointCalibration_wu(IMU, MAG, setting)
    % Copyright (C) 2024 by Chuan Huang
    % Input:
    %       imu: N x 6 matrix of IMU readings (accelerometer and gyroscope)
    %       mag: N x 3 matrix of magnetometer readings


    downsample_rate = setting.downsample_rate;


    covs = zeros(27, length(IMU));
    innovations = repmat(struct('acc',[], 'mag', []), length(IMU), 1);
    Ss = repmat(struct('acc',[], 'mag', []), length(IMU), 1);
    ANIS = 0;


    dt = 1 / setting.freq;
    [state, P] = init_state(IMU(1, 1:3), MAG(1, :), setting.EKF);
    
    % process noise
    Q = diag([(setting.gyroNoiseDensity)^2 * dt* ones(3, 1); ...
              (setting.EKF.acc_bias_noise_std )^2 * dt *  ones(3, 1); ...
              (setting.EKF.gyro_bias_noise_std)^2 * dt * ones(3, 1); ...
               zeros(12, 1); ...
               (setting.EKF.m0_noise_std)^2 * dt * ones(3, 1); ...
               (setting.EKF.g0_noise_std)^2 * dt * ones(3, 1)]);



    % measurement noise
    R_acc = setting.accNoiseDensity^2 * (setting.freq / downsample_rate) * eye(3);
    R_mag = setting.magNoiseDensity^2 * (setting.freq / downsample_rate) * eye(3);
    R = blkdiag(R_mag, R_acc);

    if setting.EKF.compute_ANIS
        update_counts = 0;
    end

    for i = 1 : size(IMU, 1)
        acc = IMU(i, 1: 3).';
       gyro = IMU(i, 4 : 6).';
        mag = MAG(i, :).';
       

        innovations(i).acc = acc - ( -state.ori' * state.g0 + state.acc_bias);
        innovations(i).mag = mag - (state.S * state.ori' * state.m0 + state.h); 
        
        if(mod(i, downsample_rate) == 0)
            H = [getHm(state); getHg(state)];
            residual = [innovations(i).mag;  innovations(i).acc];
            [state, P, S] = update(state, P, H, R, residual);
            if setting.EKF.compute_ANIS
                ANIS = ANIS + residual' * (S \ residual);
                update_counts = update_counts + 1;  
            end
            
            Ss(i).mag = diag(S(1:3, 1:3));
            Ss(i).acc = diag(S(4:6, 4:6));

        end
        % predict
        [state, P] = predict(state, P, gyro, dt, Q);
        covs(:, i) = diag(P);


    end

    if setting.EKF.compute_ANIS
        ANIS = ANIS / update_counts;
    end

    % Extract the covariance matrix for the calibration parameters
    calibres.uncertainty = P([7:9, 4:6, 10:18, 19:21], [7:9, 4:6, 10:18, 19:21]);



    calibres.S = state.S;
    calibres.h = state.h;
    calibres.m0 = state.m0;
    calibres.g0 = state.g0;
    calibres.gyro_bias = state.gyro_bias;
    calibres.acc_bias = state.acc_bias;
end


function Hg = getHg(state)
    Hg = zeros(3, 27);
    Hg(:, 1:3) =  -state.ori' * vect2skew(state.g0);
    Hg(:, 7:9) = eye(3);
    Hg(:, end - 2: end) = -state.ori';
end

function Hm = getHm(state)
    Hm = zeros(3, 27);
    Hm(:, 1:3) = state.S * state.ori' * vect2skew(state.m0);
    Hm(:, 10:10+9-1) = kron((state.ori' * state.m0)', eye(3));
    Hm(:, 19:21) = eye(3);
    Hm(:, 22:24) = state.S * state.ori';
end




function [state, P] = init_state(acc, mag, setting)
    state.ori = eye(3);
    state.gyro_bias = setting.init_gyro_bias;
    state.acc_bias = zeros(3, 1);
    state.S = eye(3);
    state.h = zeros(3, 1);
    state.m0 = mag(1, :).';
    state.g0 = -acc(1, 1:3).';

    P = diag([zeros(3, 1); ...
             setting.init_gyro_bias_std^2 * ones(3, 1); ...
              setting.init_acc_bias_std^2 * ones(3, 1); ...
              setting.init_S_std^2 * ones(9, 1); ...
              setting.init_h_std^2 * ones(3, 1); ...
              setting.init_m0_std^2 * ones(3, 1); ...
              setting.init_g0_std^2 * ones(3, 1)]);

end

function [state, P] = predict(state, P, gyro, dt, Q)
    ori = state.ori;
    
    
    rotvec = (gyro - state.gyro_bias) * dt;
    state.ori = state.ori * SO3_exp(rotvec);
    state.gyro_bias = state.gyro_bias;
    state.acc_bias = state.acc_bias;
    state.S = state.S;
    state.h = state.h;
    state.m0 = state.m0;
    state.g0 = state.g0;

    % update covariance
    F = eye(27);
    F(1:3, 4:6) = -ori*dt;

    P = F * P * F' + Q;

end

function [state, P, S] = update(state, P, H, R, delta_z)

    S = H * P * H' + R;
    K = (P * H') / S;   
    delta = K * delta_z;


    state.ori = SO3_exp(delta(1:3)) * state.ori;
    state.gyro_bias = state.gyro_bias + delta(4:6);
    state.acc_bias = state.acc_bias + delta(7:9);
    state.S = state.S + reshape(delta(10:10+9-1), 3, 3);
    state.h = state.h + delta(19:21);
    state.m0 = state.m0 + delta(22:24);
    state.g0 = state.g0 + delta(25:27);

    % update covariance
    P = (eye(size(P)) - K * H) * P;
end 

