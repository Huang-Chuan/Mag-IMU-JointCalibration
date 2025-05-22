function invJR = invRightSO3Jacobian(theta)
% This function calculates the inverse right Jacobian matrix for the SO(3) group.
% It takes an input parameter theta and returns the inverse right Jacobian matrix JL.
% Input: theta - 3x1 vector
% Output: JL - 3x3 matrix
% Reference: Sola, J. A micro Lie theory for state estimation in robotics. 2018.

a = norm(theta);
if(a < 1e-9)
    invJR = eye(3);
else
    invJR = eye(3) + 1/2 * vect2skew(theta) + (1/a^2 - (1+cos(a))/(2*a*sin(a))) * (vect2skew(theta))^2;
end
end