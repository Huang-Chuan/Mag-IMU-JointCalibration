function [calib_res] = calibrate_real(dataset, settings)
    % Copyright (C) 2024 by Chuan Huang
    % Input:

    imu = [dataset.data.y_acc; dataset.data.y_gyro ]';
    mag = dataset.data.y_mag';
    [result, covs, innovations, Ss, ANIS] = jointCalibration_wu(imu, mag, settings);


    % 
    % figure;
    % t = tiledlayout(3,1);
    % nexttile;
    % plot(arrayfun(@(x)x.mag(1),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.mag(1)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.mag(1)),Ss), 'r--')
    % legend('innovation','3 \sigma')
    % 
    % nexttile;
    % plot(arrayfun(@(x)x.mag(2),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.mag(2)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.mag(2)),Ss), 'r--')
    % nexttile;
    % plot(arrayfun(@(x)x.mag(3),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.mag(3)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.mag(3)),Ss), 'r--')
    % xlabel('time [s]')
    % 
    % 
    % title(t,'Magnetic Field Innovation')
    % 
    % 
    % 
    % figure;
    % innovations = innovations(~arrayfun(@(x) isempty(x.acc), innovations));
    % Ss = Ss(~arrayfun(@(x) isempty(x.acc), Ss));
    % t = tiledlayout(3,1);
    % nexttile;
    % plot(arrayfun(@(x)x.acc(1),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.acc(1)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.acc(1)),Ss), 'r--')
    % nexttile
    % plot(arrayfun(@(x)x.acc(2),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.acc(2)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.acc(2)),Ss), 'r--')
    % legend('innovation','3 \sigma')
    % nexttile;
    % plot(arrayfun(@(x)x.acc(3),innovations));
    % hold on
    % plot(arrayfun(@(x) 3*sqrt(x.acc(3)),Ss), 'r--')
    % plot(arrayfun(@(x) -3*sqrt(x.acc(3)),Ss), 'r--')
    % xlabel('time [s]')
    % 
    % title(t,'Accelerometer Innovation')




    calib_res.acc_bias = result.acc_bias;
    calib_res.gyro_bias = result.gyro_bias;
    calib_res.mag_bias = result.h;
    calib_res.D = result.S * norm(result.m0);
    
    % rescale the uncertainty
    if(settings.compute_covariance)
        calib_res.covar = result.uncertainty;
        calib_res.covar(10:18, 10:18) = calib_res.covar(10:18, 10:18) * norm(result.m0)^2;
    end
end
