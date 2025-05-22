function opt_val = gauss_newton_opt(init_value, data, settings)

    sigma_acc = settings.accNoiseDensity * sqrt(settings.freq /settings.preintegration_num);
    sigma_mag = settings.magNoiseDensity * sqrt(settings.freq /settings.preintegration_num);
    sigma_gyro = settings.gyroNoiseDensity * sqrt(settings.freq);
    preintegration_num = settings.preintegration_num;
             g = settings.g;

    % parse initial values
    D_hat = init_value.D(:);
    R_hat_Gauss = init_value.init_ori;
    dip_angle_hat = init_value.dip_angle;
    acc_bias_hat = init_value.acc_bias;
    mag_bias_hat = init_value.mag_bias;
    gyro_bias_est_err = zeros(3, 1);
    nominal_gyro_bias = init_value.gyro_bias;

    % parse data



    y_acc = data.y_acc;
    y_mag = data.y_mag;
    y_gyro = data.y_gyro;


    % IMU pre-integration
    [deltaR_ij, noise_cov, dRij_db] = gyroPreIntegration(y_gyro(:,1:end-1), nominal_gyro_bias,  sigma_gyro*eye(3), preintegration_num, settings.dt);
    
    B_array = zeros(3, 3, size(deltaR_ij, 3));
    
    for i = 1 : size(deltaR_ij, 3)
        B_array(:,:,i) = chol(noise_cov(:,:,i));    
    end

    % downsample accelorometer and magnetometer data
    y_acc_preInt = y_acc(:, 1:preintegration_num:end);
    y_mag_preInt = y_mag(:, 1:preintegration_num:end);


    for i = 1 : size(deltaR_ij, 3)
        R_hat_Gauss(:,:, i+1) = R_hat_Gauss(:,:, i) * deltaR_ij(:,:,i);
    end

    num_ori = size(R_hat_Gauss, 3);
    %res = zeros(9*num_ori - 3, 1);
    %grad_res = zeros(9*num_ori - 3, 3*num_ori + 9 + 9 + 1);


    dx = 1e6;

    while (norm(dx) > 1e-6)
        D_hat_mat = reshape(D_hat, 3, 3);
        m_hat = [0; cos(dip_angle_hat); -sin(dip_angle_hat)];
   
        % check if imu pre-integration is needed
        % if(norm(gyro_bias_est_err) > 0.4)
        %     [deltaR_ij, noise_cov, dRij_db] = gyroPreIntegration(y_gyro(:,1:end-1), nominal_gyro_bias,  sigma_gyro*eye(3), preintegration_num, settings.dt);
        % end

        % for j = 1 : num_ori
        %     grad_res(9*j-8:9*j-3, 3*j-2:3*j) = -[1/sigma_acc* vect2skew(-R_hat_Gauss(:,:,j)'*-g)*-R_hat_Gauss(:,:,j)'; ...
        %                                     1/sigma_mag * D_hat_mat*vect2skew(-R_hat_Gauss(:,:,j)'*m_hat)*-R_hat_Gauss(:,:,j)'];
        %     grad_res(9*j-8:9*j-6, 3*num_ori+9+1:3*num_ori+9+3) = -1/sigma_acc * eye(3);
        %     grad_res(9*j-5:9*j-3, 3*num_ori+1:3*num_ori+9) = -1/sigma_mag * kron((R_hat_Gauss(:,:,j)'*m_hat)',eye(3));
        %     grad_res(9*j-5:9*j-3, 3*num_ori+9+4:3*num_ori+9+6) = -1/sigma_mag * eye(3);
        %     grad_res(9*j-5:9*j-3, end) = -1/sigma_mag * D_hat_mat * R_hat_Gauss(:,:,j)' * [0; -sin(dip_angle_hat); -cos(dip_angle_hat)];

        %     res(9*j-8:9*j-3) = [1/sigma_acc*(y_acc_preInt(:,j) - R_hat_Gauss(:,:,j)' * -g - acc_bias_hat);...
        %                     1/sigma_mag*(y_mag_preInt(:,j) - D_hat_mat * R_hat_Gauss(:,:,j)' * m_hat - mag_bias_hat)];
        %     if j < num_ori
        %         B = squeeze(B_array(:,:,j));
        %         temp = log_SO3((deltaR_ij(:,:,j) * SO3_exp(dRij_db(:,:,j) * gyro_bias_est_err))' * R_hat_Gauss(:,:,j)'*R_hat_Gauss(:,:,j+1));

        %         res(9*j-2:9*j) = B\temp;
        %         grad_res(9*j-2:9*j, 3*j-2:3*j) = -B\(invRightSO3Jacobian(temp)* R_hat_Gauss(:,:,j+1)');
        %         grad_res(9*j-2:9*j, 3*j+1:3*j+3) = B\(invRightSO3Jacobian(temp)* R_hat_Gauss(:,:,j+1)');
        %         grad_res(9*j-2:9*j, 3*num_ori+9+7:3*num_ori+9+9) = B\(-invRightSO3Jacobian(temp)* expm(vect2skew(temp))' * rightSO3Jacobian(dRij_db(:,:,j) * gyro_bias_est_err) * dRij_db(:,:,j));
        %     end
        % end

        [grad_res, res] = construct_sparse_grad_and_res(num_ori, R_hat_Gauss, g, sigma_acc, sigma_mag, D_hat_mat, m_hat, acc_bias_hat, mag_bias_hat, dip_angle_hat, B_array, deltaR_ij, dRij_db, gyro_bias_est_err, y_acc_preInt, y_mag_preInt);

        [CC, RR] = qr(grad_res, -res);
        dx = RR\CC;
        
        
        % update the variables
        for j = 1 : num_ori
            R_hat_Gauss(:,:,j) = SO3_exp((dx(3*j-2:3*j))) * R_hat_Gauss(:,:, j);
        end
        D_hat = D_hat + dx(3*num_ori+1:3*num_ori+9);
        acc_bias_hat = acc_bias_hat + dx(3*num_ori+9+1:3*num_ori+9+3);
        mag_bias_hat = mag_bias_hat + dx(3*num_ori+9+4:3*num_ori+9+6);
        gyro_bias_est_err = gyro_bias_est_err + dx(3*num_ori+9+7:3*num_ori+9+9);
        gyro_bias_hat = nominal_gyro_bias + gyro_bias_est_err;
        dip_angle_hat = dip_angle_hat + dx(end);


    end


    if(settings.compute_covariance)
        % construct measurement covariance matrix
        %[~,RR] = qr(grad_res);
        % remove all zero rows 
        RR = RR(1:size(grad_res,2),:);
        % since we are interested in computing the covariance of the calibration parameters, we adopt the approach in ISAM 
        RR = RR(end-18:end,end-18:end);
        opt_val.covar = full( (RR' * RR) \ eye(size(RR))) ;
        % reorder the covariance matrix to match the order of the parameters [acc_bias; gyro_bias; mag_bias; D]
        opt_val.covar = opt_val.covar([10:12, 16:18, 1:9,13:15], [10:12, 16:18, 1:9, 13:15]);
    end



    opt_val.acc_bias = acc_bias_hat;
    opt_val.mag_bias = mag_bias_hat;
    opt_val.D = reshape(D_hat, 3, 3);
    %opt_val.R_hat = R_hat_Gauss;
    opt_val.dip_angle = dip_angle_hat;
    opt_val.gyro_bias = gyro_bias_hat;
    
    %display(R_hat_Gauss(:,:,1));

end