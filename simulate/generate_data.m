function [y_acc, y_mag, y_gyro] = generate_data(calib_paras, settings)
% generate data for accelerometer, magnetometer and gyroscope
% calib_paras: calibration parameters
%              acc_bias: accelerometer bias
%              mag_bias: magnetometer bias
%                     D: magnetometer distortion matrix
%             dip_angle: dip angle
%              gyro_bias: gyroscope bias
% settings: settings for the data generation
%              g: gravity
%              K: number of data points
%              rotationTrajectory: rotation trajectory
%              dt: time interval
%              sigma_acc: accelerometer noise
%              sigma_mag: magnetometer noise
%              sigma_gyro: gyroscope noise

% y_acc: accelerometer data
% y_mag: magnetometer data
% y_gyro: gyroscope data

acc_bias = calib_paras.acc_bias;
mag_bias = calib_paras.mag_bias;
D = calib_paras.D;
dip_angle = calib_paras.dip_angle;
gyro_bias = calib_paras.gyro_bias;


g = settings.g;
K = settings.K;
R = settings.rotationTrajectory;
dt = settings.dt;
sigma_acc = settings.sigma_acc;
sigma_mag = settings.sigma_mag;
sigma_gyro = settings.sigma_gyro;


y_acc = zeros(3, K);
y_mag = zeros(3, K);
y_gyro = zeros(3, K - 1);
for i = 1 : K
    y_acc(:, i) =  R(:,:,i)' * -g + acc_bias + sigma_acc * randn(3, 1);
    y_mag(:, i) =  D * R(:,:,i)' * [0 cos(dip_angle) -sin(dip_angle)]' + mag_bias + sigma_mag * randn(3, 1);
    if i < K
        y_gyro(:, i) = log_SO3(R(:,:,i)' * R(:,:,i+1)) / dt + gyro_bias + sigma_gyro * randn(3, 1);
    else
        y_gyro(:, i) = [0;0;0];
    end
end

end