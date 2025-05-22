function [datasets] = lowpass_data(datasets, settings)
    
    downsample_rate = settings.downsample_rate;
    
    
    for i = 1 : length(datasets)
        IMU = [datasets(i).data.y_acc' datasets(i).data.y_gyro'];
        MAG = datasets(i).data.y_mag';
        
        % low pass filter
        if (downsample_rate > 1)
            num_padding = settings.num_padding;
            acc = [flip(IMU(1:num_padding, 1:3)); IMU(:, 1:3); flip(IMU(end-num_padding+1:end, 1:3))];
            filtered_acc = lowpass(acc, (settings.freq/downsample_rate)/2, settings.freq);
            IMU(:, 1:3) = filtered_acc(num_padding+1:end-num_padding, :);
            mag = [flip(MAG(1:num_padding, :)); MAG; flip(MAG(end-num_padding+1:end, :))];
            filtered_mag = lowpass(mag, (settings.freq/downsample_rate)/2, settings.freq);
            MAG = filtered_mag(num_padding+1:end-num_padding, :);
        end

        
        datasets(i).data.y_acc = IMU(:, 1:3).';
        datasets(i).data.y_gyro = IMU(:, 4:6).';
        datasets(i).data.y_mag = MAG.';
    end






end