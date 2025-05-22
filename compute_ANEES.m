function [ANEES] = compute_ANEES(calibres, calib_paras)
    ANEES = 0;

    acc_bias_err = [calibres.acc_bias] - [calib_paras.acc_bias];
    gyro_bias_err = [calibres.gyro_bias] - [calib_paras.gyro_bias];
    D_err = reshape([calibres.D] - [calib_paras.D], 9, []);
    mag_bias_err = [calibres.mag_bias] - [calib_paras.mag_bias];
    

   
    for i = 1 : length(calib_paras)
        err_vec = [ acc_bias_err(:, i); gyro_bias_err(:, i); D_err(:, i); mag_bias_err(:, i)];
        ANEES = ANEES + err_vec.' * (calibres(i).covar \ err_vec);
    end


    ANEES = ANEES / (length(calibres) * length(err_vec));
end