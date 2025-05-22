function rotm = q2r(q)
% q2r - Convert quaternion to rotation matrix
% Copyright 2024 by Chuan Huang
% Input:
%       q: 4 x 1 quaternion
% Output:
%       rotm: 3 x 3 rotation matrix


    q_w = q(1);
    q_x = q(2);
    q_y = q(3);
    q_z = q(4);

         

    rotm = [2*(q_w^2+q_x^2) - 1,2*(q_x*q_y-q_w*q_z),2*(q_x*q_z+q_w*q_y);
           2*(q_x*q_y+q_w*q_z), 2*(q_w^2+q_y^2) - 1,2*(q_y*q_z-q_w*q_x);
           2*(q_x*q_z-q_w*q_y), 2*(q_y*q_z+q_w*q_x), 2*(q_w^2+q_z^2) - 1];


end