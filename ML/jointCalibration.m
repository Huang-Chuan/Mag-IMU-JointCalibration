function [theta] = jointCalibration(imu, mag, x0, P0, setting)
% Copyright (C) 2024 by Chuan Huang
% Input:
%       imu: N x 6 matrix of IMU readings (accelerometer and gyroscope)
%       mag: N x 3 matrix of magnetometer readings
%        x0: 4 x 1 initial rotation w.r.t. the navigation frame
%        P0: 3 x 3 initial covariance matrix
%       setting: struct with the following fields
%           g0: 3 x 1 gravity vector
%           freq: sampling frequency
%           accNoiseDensity: accelerometer noise density (m/s^2/sqrt(Hz))
%           gyroNoiseDensity: gyroscope noise density (rad/s/sqrt(Hz))
%           magNoiseDensity: magnetometer noise density (uT/sqrt(Hz))
%           minAcc: minimum accelerometer reading norm
%           maxAcc: maximum accelerometer reading norm
%           ml_samples: number of samples for maximum likelihood estimation
% Output:
%       theta: 9 x 1 + 3 + 2N vector

    Ra =  setting.accNoiseDensity^2 * setting.freq / setting.downsample_rate * eye(3);
    Rw =  setting.gyroNoiseDensity^2 * setting.freq * eye(3); 
    Rm =  setting.magNoiseDensity^2 * setting.freq / setting.downsample_rate * eye(3);
    
    %init_theta = computemagpara(imu, mag, setting, Ra, Rw);
    init_theta = [computemagpara(imu, mag, setting, Ra, Rw); setting.ML.init_acc_bias; setting.ML.init_gyro_bias];

    % options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
    %                        'UseParallel',true, 'PlotFcn', 'optimplotfval');
    options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton', ...
                           'UseParallel',true);


    obj = @(theta) calneglkh(imu, mag, theta, x0, P0, Ra, Rw, Rm, setting);
    theta = fminunc(obj, init_theta, options);


end

function [initTheta] = computemagpara(imu, mag, setting, Ra, Rw)
    % First step: Ellipsoid fitting
    % [D_tilde, b] = ellipsoid_fitting(mag.');
    % mag_m = (D_tilde \ (mag.' - b)).';

    % [A, b, expmfs] = magcal(mag);
    % D_tilde = expmfs * inv(A.');
    % mag_m = ((mag - b) * A ) / expmfs;
    % b = b.';

    [D_tilde, b] = ml_estimate(mag.', setting);
    mag_m = (D_tilde \ (mag.' - b)).';

    Rnb = computepose(imu, setting, Ra, Rw);
    % set up nonlinear least square optimization
    M.RD = rotationsfactory(3, 1);
    M.vz = euclideanfactory(1, 1);
    problem.M = productmanifold(M);
    problem.cost = @(x) cost(x, mag_m.', Rnb);
    problem = manoptAD(problem);
    options.maxiter = 30;
    options.verbosity = 0;
    [sol, ~, ~] = trustregions(problem, [], options);
    D_h = D_tilde * sol.RD;
    dipAngle_h = atan2(-sol.vz, sqrt(1-sol.vz^2));
    %disp(rad2deg(dipAngle_h))
    initTheta = [D_h(:); b; dipAngle_h];

end

function f = cost(x, mag_m, rotmtx_h)
    RD= x.RD;
    v = x.vz;
    tmp = pagemtimes(pagemtimes(rotmtx_h, RD.'), reshape(mag_m, [size(mag_m, 1), 1, size(mag_m, 2)]));
    f = 1/2 * sum((v - tmp(3, :, :)).^2, 'all');
end



function [Rnb] = computepose(imu, setting, Ra, Rw)
    g0 = setting.g0;
    dt = 1 / setting.freq;
    % Initialize pose using accelerometer reading
    yaw   = 0;                                                  % since the yaw is not observable, we set it to 0
    pitch = atan2(-imu(1, 1), sqrt(imu(1, 2)^2 + imu(1, 3)^2));
    roll  = atan2(imu(1, 2), imu(1, 3));
    
    x0 = euler2quaternion([yaw pitch roll]);
    P0 = (setting.ML.init_pose_uncertainty)^2 * eye(3);                                 % initial covariance matrix
    
    Rnb = zeros(3, 3, length(imu));
    Rnb(:, :, 1) = q2r(x0);
    x = x0;
    P = P0;
    for i = 1 : size(imu, 1)
        acc = imu(i, 1 : 3).' - setting.ML.init_acc_bias;
       gyro = imu(i, 4 : 6).' - setting.ML.init_gyro_bias;
        % only update if acc is within resonable range
        accAmplitude = vecnorm(acc);
        if accAmplitude < setting.maxAcc && accAmplitude > setting.minAcc
            H = getH(x, g0);
            S = H * P * H' + Ra;
            SI = inv(S);
            K = P * H' * SI;   
            delta_z = acc - q2r(x).' * g0;
            delta_x = K * delta_z;
            x = q_R([1; 1/2*delta_x]) *  x;
            x = x / norm(x); 
            P = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K * Ra * K';
            % reset error state 
            G = eye(3) - vect2skew(1/2*delta_x);
            P = G * P * G';
        end
        % predict
        rotvec = gyro * dt;
        rot12m = SO3_exp(rotvec);
        x = q_R(rotvec2quat(rotvec)) *  x;

        
        
        F = rot12m';
        Fi = eye(3) * dt;
        P = F * P * F' + Fi * Rw * Fi';

        if(i ~= length(imu)) 
            Rnb(:, :, i + 1) = q2r(x); 
        end

    end
end


function H = getH(x, g0)
    [Q0, Q1, Q2, Q3] = dQqdq(x);
    Hx  = [Q0.' * g0, Q1.' * g0, Q2.' * g0, Q3.' * g0];
    HDeltax = Sq(x);
    H = Hx * HDeltax;
end

