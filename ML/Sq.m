function S=Sq(q)
% Copyright 2024 by Chuan Huang
% This function is the Jacobian, \partial{q \times \delta q} / \partial {\delta \theta}
% Input: q, a 4x1 quaternion
% Output: S, a 4x3 matrix
% Reference: Eq. 281 in Quaternion Kinematics for the error-state KF, Joan Sola, 2017
   q0=q(1);   q1=q(2);   q2=q(3);   q3=q(4);
   S=1/2*[-q1 -q2 -q3;
       q0 -q3  q2;
       q3  q0 -q1;
      -q2  q1  q0];
end
