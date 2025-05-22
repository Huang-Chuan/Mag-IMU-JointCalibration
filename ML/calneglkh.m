function [total_neg_lkh] = calneglkh(IMU, MAG, theta, x0, P0, Ra, Rw, Rm, setting)
    
    dt = 1 / setting.freq;

    % low pass filter
    downsample_rate = setting.downsample_rate;



    D = reshape(theta(1:9), 3, 3);
    h = theta(10:12);
    dip_angle = theta(13);
    acc_bias = theta(14:16);
    gyro_bias = theta(17:19);


    % remove bias from measurements
    IMU(:, 1:3) = IMU(:, 1:3) - acc_bias.';
    IMU(:, 4:6) = IMU(:, 4:6) - gyro_bias.';
    MAG = MAG - h.';

    x = x0;
    P = P0;

    g0 = setting.g0;
    m0 = [0; cos(dip_angle); -sin(dip_angle)];
    
    total_neg_lkh = 0;


    R = blkdiag(Ra, Rm);


    for i = 1 : size(IMU, 1)
        acc = IMU(i, 1 : 3).';
        mag = MAG(i, :).';
        
        Rn2b = q2r(x(1:4)).';

        
        if (mod(i, downsample_rate) == 0)
            %H = [getHg(x, g0); getHm(x, m0, D)];
            H = getH(x, g0, m0, D);
            delta_z = [acc - Rn2b * g0; mag - (D * Rn2b * m0)];
            [x, P, neg_lkh] = measurement_update(x, P, H, delta_z, R);
        else
            neg_lkh = 0;
        end
        
        % this is to make sure the likelihood is accumulated when the filter is stable
        if i > setting.warmup
            total_neg_lkh = total_neg_lkh + neg_lkh;
        end 
        
        
       [x, P] = time_update(x, P, IMU(i, 4 : 6).', Rw, dt);

    end
end




function Hg = getHg(x, g0)
    q = x(1:4);
    [Q0, Q1, Q2, Q3] = dQqdq(q);
    Hq = [Q0.' * g0, Q1.' * g0, Q2.' * g0, Q3.' * g0];
    HDeltax = Sq(q);
    Hg = Hq * HDeltax;
end




function Hm = getHm(x, m0, D)
    q = x(1:4);
    [Q0, Q1, Q2, Q3] = dQqdq(q);
    Hq = [D * Q0.' * m0, D * Q1.' * m0, D * Q2.' * m0, D * Q3.' * m0];
    HDeltax = Sq(q);
    Hm = Hq * HDeltax;
end


function H = getH(x, g0, m0, D)
    q = x(1:4);
    [Q0, Q1, Q2, Q3] = dQqdq(q);
    HDeltax = Sq(q);
    
    Hg = [Q0.' * g0, Q1.' * g0, Q2.' * g0, Q3.' * g0] * HDeltax;
    Hm = [D * Q0.' * m0, D * Q1.' * m0, D * Q2.' * m0, D * Q3.' * m0] * HDeltax;

    H = [Hg; Hm];
end




function [x, P] = time_update(x, P, gyro, R_omega, dt)
    % update quaternion
    x(1:4) = q_R(rotvec2quat(gyro * dt)) *  x(1:4);
    
    rotm = SO3_exp(gyro * dt);

    F = rotm';
    %Fi = eye(3) * dt;
    P = F * P * F' + dt^2 * R_omega;

end


function [x, P, neg_log] = measurement_update(x, P, H, delta_z, R)
    S = H * P * H' + R;
    %SI = inv(S);
    K = (P * H') / S;   

    neg_log = 1/2 * delta_z' * (S \ delta_z) + log(det(S));

    delta_x = K * delta_z;

    delta_q = delta_x(1:3);

    x(1:4) = q_R([1; 1/2*delta_q]) *  x(1:4);
    x(1:4) = x(1:4) / norm(x(1:4)); 

    P = (eye(size(P)) - K * H) * P;
end