function q_conjugated=q_conj(q)
% Copyright (C) 2024 by Chuan Huang
% This function is used to calculate the conjugated quaternion
% Input: q, a 4x1 quaternion
% Output: q_conjugated, a 4x4 matrix
    q_conjugated = [q(1); -q(2); -q(3); -q(4)];
end