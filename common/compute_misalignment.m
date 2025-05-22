function [R, o_gyro] = compute_misalignment(mag, gyro, setting)
% Copyright (C) 2024 by Chuan Huang
% This function computes the misalignment between the magnetometer and the gyroscope biases
% Inputs:
%       mag: 3 x N matrix of magnetometer readings
%       gyro: 3 x N matrix of gyroscope readings
%       setting: struct with the following fields    
% Output:
%       R: 3 x 3 matrix, misalignment matrix
%       o_gyro: 3 x 1 vector, gyroscope bias
% Reference:
%      [1] On Misalignment Between Magnetometer and Inertial Sensors, Yuanxin Wu, et al., 2016

    N = size(mag, 2);
    dt = 1/setting.freq;
    % Get interval length in seconds
    MM = floor(setting.alignment_interval  / dt);
    % Compute W_k and M_k based on Equation (12) and (13)
    W = zeros(3 * length(1:MM:N-MM), 9);
    M = zeros(3 * length(1:MM:N-MM), 3);
    for i = 1 : MM : N - MM
        W_k = zeros(3, 9);
        M_k = zeros(3, 3);
        for j = 0 : MM - 1
            W_k = W_k + dt / 2 * (kron(gyro(:, i + j)', vect2skew(mag(:, i + j))) + ...
                                  kron(gyro(:, i + j + 1)',  vect2skew(mag(:, i + j + 1))));
            M_k = M_k - dt / 2 * (vect2skew(mag(:, i + j)) + vect2skew(mag(:, i + j + 1)));
        end
        W(3*i-2:3*i, :)= W_k;
        M(3*i-2:3*i, :)= M_k;
    end


    A = [W, M];

    % Compute differences of magnetometer readings between time k and k + 1
    diff_mag = zeros(3 * length(1:MM:N-MM), 1);
    for i = 1 : MM : N - MM
        diff_mag(3*i-2:3*i, :) =  mag(: , i + MM) - mag(:, i);
    end


    

    % compute initial R and biases
    % initial_value = A \ diff_mag;
    % R0 = reshape(initial_value(1:9), 3, 3);
    % R0 = Gram_Schmidt(R0);
    % if(det(R0) < 0)
    %     R0 = -R0;
    % end
  

    % o_gyro0 = initial_value(10:12);


    % x0.R = R0;
    % x0.h = o_gyro0;

    x0.R = eye(3);
    x0.h = zeros(3, 1);

    % set up nonlinear least square optimization
    Manifold.R = rotationsfactory(3, 1);
    Manifold.h = euclideanfactory(3, 1);
    problem.M = productmanifold(Manifold);
    problem.cost = @(x) cost(x, A, diff_mag);
    problem = manoptAD(problem);
    options.maxiter = 30;
    options.verbosity = 0;
    [sol, ~, ~] = trustregions(problem, x0, options);

    R = sol.R;
    o_gyro = sol.h;


end


function f = cost(x, A, diff_mag)
    R = x.R;
    h = x.h;

    f =  sum((A * [R(:); h(:)] - diff_mag).^2, 'all');
end




function Q = Gram_Schmidt(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1 : n
        v = A(:, j);
        for i = 1 : j - 1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end