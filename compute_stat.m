function [stat] = compute_stat(calibres, calib_paras)
% compute mean square error
%N = length(calibres);
acc_bias_err = sqrt(mean(sum(([calibres.acc_bias] - [calib_paras.acc_bias]).^2)));
disp(['acc_bias_err: ', num2str(acc_bias_err)]);

gyro_bias_err = sqrt(mean(sum(([calibres.gyro_bias] - [calib_paras.gyro_bias]).^2)));
disp(['gyro_bias_err: ', num2str(gyro_bias_err)]);

mag_bias_err = sqrt(mean(sum(([calibres.mag_bias] - [calib_paras.mag_bias]).^2)));
disp(['mag_bias_err: ', num2str(mag_bias_err)]);

D_err =  sqrt(mean(sum(reshape(([calibres.D] - [calib_paras.D]).^2, 9, []))));
disp(['D_err: ', num2str(D_err)]);


stat.acc_bias_err = acc_bias_err;
stat.gyro_bias_err = gyro_bias_err;
stat.mag_bias_err = mag_bias_err;
stat.D_err = D_err;

end