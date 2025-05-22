function JL = leftSO3Jacobian(theta)
% This function calculates the left Jacobian matrix for the SO(3) group.
% It takes an input parameter theta and returns the left Jacobian matrix JL.
% Input: theta - 3x1 vector
% Output: JL - 3x3 matrix
% Reference: Sola, J. A micro Lie theory for state estimation in robotics. 2018.

    a = norm(theta);
    if(a < 1e-8)
        JL = eye(3);
    else
        JL = eye(3) + (1 - cos(a))/a^2 * vect2skew(theta) + (a - sin(a))/a^3 * (vect2skew(theta))^2;
    end
end