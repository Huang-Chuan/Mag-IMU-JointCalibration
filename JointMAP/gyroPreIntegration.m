function [deltaR_ij, noise_cov, dDeltaR_ij_dGyro_bias] = gyroPreIntegration(gyro_data, gyro_bias, gyro_noise, L, dt)
% gyroPreIntegration: Pre-integrate the gyroscope data to get the relative
% rotation between the poses
% INPUTS 
%   gyro_data: 3xN matrix of gyroscope data
%   gyro_bias: 3x1 vector of gyroscope bias
%   gyro_noise: 3x3 covariance matrix of gyroscope noise
%   L: number of data points to integrate
%   dt: time step
% OUTPUTS
%   deltaR_ij: 3x3x(floor(N/L)) matrix of relative rotations
%   noise_cov: 3x3x(floor(N/L)) matrix of noise covariance

% Reference: On-Manifold Preintegration for Real-Time Visual-Inertial Odometry, Christian Forster et al. IEEE Transactions on Robotics, 2017

N = size(gyro_data, 2);

M = floor(N/L); % number of relative rotations
deltaR_ij = repmat(eye(3, 3), [1, 1, M]);
noise_cov = zeros(3, 3, M);
dDeltaR_ij_dGyro_bias = zeros(3, 3, M);

i = 1;
for kk = 1 : M
    j = i + L;
    for k = j - 1 : -1 : i 
        J =  rightSO3Jacobian((gyro_data(:, k) - gyro_bias)*dt);
        G = (deltaR_ij(:,:,kk)' * J * dt); 
        dDeltaR_ij_dGyro_bias(:, :, kk) = dDeltaR_ij_dGyro_bias(:, :, kk) - G;
        noise_cov(:, :, kk) = noise_cov(:, :, kk)  + G * gyro_noise * G';
        deltaR_ij(:, :, kk) = SO3_exp((gyro_data(:, k) - gyro_bias) * dt) * deltaR_ij(:, :, kk) ;
    end
    i = j;
end

end