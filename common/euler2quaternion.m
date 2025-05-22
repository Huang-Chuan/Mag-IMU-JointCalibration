function q = euler2quaternion(euler)
% Copyright (C) 2024 by Chuan Huang
% Convert Euler angles to quaternion (Rotation sequence: ZYX)
% Input: euler, a 3x1 vector of Euler angles in radians [yaw pitch roll]
% Output: q, a 4x1 quaternion
    phi = euler(1);
    theta = euler(2);
    psi = euler(3);
    q = [cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
         sin(psi/2)*cos(theta/2)*cos(phi/2) - cos(psi/2)*sin(theta/2)*sin(phi/2);
         cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
         cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2)];
end