function [D_hat, o_hat] = ellipsoid_fitting(y)
% Copyright (C) 2024 by Chuan Huang
% Input: 
%       y: 3 x N matrix of magnetometer readings
% Output: 
%       D_hat: 3 x 3 matrix, o_hat: 3 x 1 vector
% Reference:
% [1] Calibration of a magnetometer in combination with inertial sensors, Manon Kok, et al., 2012.
% [2] Magnetometer Calibration Using Inertial Sensors, Manon Kok, et al., 2016.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The estimated D_hat and o_hat are not optimal in terms of bringing the magnetometer readings to the unit sphere.
% In fact, they are just an approximated solution to the optimization problem || ||D^{-1} (y - o)||_2^2 - 1 ||_2^2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    N = size(y, 2);
    M = zeros(N, 13);

    for i = 1 : size(y, 2)
        M(i, :) = [kron(y(:,i)', y(:, i)') y(:, i)' 1];
    end


    A = sdpvar(3, 3);
    b = sdpvar(3, 1);
    c = sdpvar(1, 1);
    J = M * [A(:); b; c];

    % This part requires YALMIP and MOSEK
    % Optimization problem in [2] (Eq. 23)
    Objective = J' * J;
    Constraints = [A>=0, A(1,1)==1];
    ops = sdpsettings('solver','mosek', 'verbose', 0);

    sol = optimize(Constraints,Objective,ops);
    
    if sol.problem == 0
        A_hat = value(A);
        b_hat = value(b);
        c_hat = value(c);
    else
        display('Hmm, something went wrong!');
        sol.info
        yalmiperror(sol.problem)
    end


    beta = 1/(1/4 * b_hat' * (A_hat \ b_hat) - c_hat);
    
    % Here Eq. 23a in [1] is used because that in [2] is incorrect.
    D_inv = chol(beta*A_hat);
    
    D_hat = inv(D_inv);
    o_hat = -1/2 * A_hat \ b_hat;


end