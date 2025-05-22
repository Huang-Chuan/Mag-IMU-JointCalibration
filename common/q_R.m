function q_R_matrix=q_R(q)
% Copyright (C) 2022 by Frida Viset
% Comments: Chuan Huang
% This function is used to calculate the right matrix of a quaternion
% q1 \times q2 = q_R(q2) q1
% Reference: eq. 17 in Quaternion Kinematics for the error-state KF, https://arxiv.org/abs/1711.02508, [v1] Fri, 3 Nov 2017 11:53:43 UTC (1,416 KB)
% Input: q, a 4x1 quaternion
% Output: q_R_matrix, a 4x4 matrix

qx=q(2);
qy=q(3);
qz=q(4);
OMEGA=[0 -qx -qy -qz;
    qx 0 qz -qy;
    qy -qz 0 qx;
    qz qy -qx 0];
q_R_matrix=(q(1)*eye(4)+OMEGA);
end