function [calibres] = calibrate_real(datasets, settings)

    settings.preintegration_num = settings.downsample_rate;
    settings.g = -settings.g;

    nominal_gyro_bias = settings.nominal_gyro_bias;


    init_value.init_ori = q2r(settings.ML.init_pose);
    init_value.dip_angle = deg2rad(settings.dip_angle); % the dip angle of the magnetic field in rad
    init_value.acc_bias = zeros(3, 1);                 % the initial value of the accelerometer bias
    init_value.mag_bias = zeros(3, 1);                 % the initial value of the magnetometer bias
    init_value.gyro_bias = nominal_gyro_bias;          % the initial value of the gyro bias (rad/s) IMPORTANT: SET IT TO THE AVERAGE VALUE OF THE GYROSCOPE MEASUREMENTS WHEN 
                                                       % THE SENSOR IS STATIONARY

    data = datasets.data;
    %[A,~,expmfs] = magcal(data.y_mag');
    %init_value.D = expmfs*inv(A');
    
    [DI, o_hat] = ml_estimate(data.y_mag, settings);
    [R_D, ~] = compute_misalignment(DI\(data.y_mag - o_hat), data.y_gyro - nominal_gyro_bias, settings);

    init_value.D = DI * R_D;
    init_value.mag_bias = o_hat;



    opt_val = gauss_newton_opt(init_value, data, settings);
    calibres = opt_val;

end