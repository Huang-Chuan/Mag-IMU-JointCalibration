function w = log_SO3(R)
% This function computes the logarithm of a rotation matrix R in SO(3)
% Reference: Exploring SO3 logarithmic map: degeneracies and derivatives
    theta = acos((trace(R) - 1) / 2);
    if abs(3 - trace(R)) < 1e-8
        w = 0.5*(1+1/6*theta^2+7/360*theta^4)*[R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
    else
        w = theta / (2 * sin(theta)) *  [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
    end
end
