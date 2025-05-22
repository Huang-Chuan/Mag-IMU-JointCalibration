function [D_hat, o_hat] = ml_estimate(y, setting)
% Copyright 2024 by Chuan Huang
% Input: 
%       y: 3 x N matrix of magnetometer readings
% Output:
%       D_hat: 3 x 3 matrix, o_hat: 3 x 1 vector

    % downsample
    N = setting.intrisic_calib_samples;
    n_step = floor((size(y, 2) - 1) / (N - 1));
    y = y(:, 1 : n_step : 1 + n_step * (N - 1) );


    theta_0 = zeros(12 + 2 * N, 1);
    theta_0(1) = norm(y(:, 1));
    theta_0(5) = norm(y(:, 1));
    theta_0(9) = norm(y(:, 1));
    
    % compute initial angles
    yaw = atan2(y(1,:), vecnorm(y(2:3, :)));
    roll = atan2(-y(3,:), y(2,:));
    theta_0(13:2:end) = yaw;
    theta_0(14:2:end) = roll;
    
    % problem
    problem.objective = @(theta) costFunctionMag(theta, y, N);
    problem.x0 = theta_0;
    problem.solver = 'lsqnonlin';
    problem.options = optimoptions('lsqnonlin', 'Display', 'none', 'SpecifyObjectiveGradient', true);
    theta_est = lsqnonlin(problem);
    
    
    % compute calibrated measurements
    D_hat = reshape(theta_est(1:9), 3, 3);
    o_hat = theta_est(10:12);

    % plot the differences in initial euler angles and the final euler
    % angles
    yaw_final = theta_est(13:2:end);
    roll_final = theta_est(14:2:end);

    % figure;
    % subplot(1,2,1)
    % plot(1:N, rad2deg(yaw), 1:N, rad2deg(yaw_final));
    % legend('initial', 'final')
    % title('yaw angle')

    % subplot(1,2,2)
    % plot(1:N, rad2deg(roll), 1:N, rad2deg(roll_final));
    % legend('initial', 'final')
    % title('roll angle')



end



function [z, J] = costFunctionMag(theta, y, N)
% Copyright 2024 by Chuan Huang
%   Cost function for magnetometer calibration
%   Inputs:
%           theta: (9 x 1 + 3 + 2N) vector  
%                   1 : 9    Distortion matrix D
%                   10 : 12  Bias b
%                   13 : end [yaw_1, roll_1, yaw_2, roll_2,...]
%               y: 3 N x 1 vector   measurements from accelerometer
%               N: 1 x 1 scalar     number of measurements
%   Outputs:
%               z: 3 N x 1 vector   residual
    
        D = reshape(theta(1:9), 3, 3);
        o = theta(10:12);

        ori = reshape(theta(13:end), 2, N);




        m_b = [sin(ori(1, :)); ...
               cos(ori(1, :)) .* cos(ori(2, :)); ...
               -cos(ori(1, :)) .* sin(ori(2, :));];
    
        
        
        z = y - D * m_b - o;
        z = z(:);


        J = zeros(3 * N, 12 + 2 * N);
        for i = 1 : N
            J(3 * i - 2 : 3 * i, 1:9) = -kron(m_b(:, i)', eye(3));
            J(3 * i - 2 : 3 * i, 10:12) = -eye(3);
            J(3 * i - 2 : 3 * i, 12 + 2 * i - 1) = -D * [cos(ori(1,i)); -sin(ori(1,i)) * cos(ori(2,i)); sin(ori(1,i)) * sin(ori(2,i))];
            J(3 * i - 2 : 3 * i, 12 + 2 * i) = -D * [0; -cos(ori(1,i)) * sin(ori(2,i)); -cos(ori(1,i)) * cos(ori(2,i))];
        end




    end