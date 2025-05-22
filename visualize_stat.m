function [] = visualize_stat(ekf_times, statEKF, jointmap_times, statMAP, ML_times, statML, sweeplist, type)
    if strcmp(type, 'freq')
        xlabel_str = 'Sampling Frequency (Hz)';
    else
        xlabel_str = 'Sampling Rate Ratio';
    end

    num_sweeps = length(sweeplist);
    methods = {'Proposed', 'Kok et al. (20 cores)', 'Wu et al.'};  % new order
    data = zeros(num_sweeps, length(methods));

    % Reordered colors to match: Proposed, Kok et al., Wu et al.
    bar_colors = [0.4660 0.6740 0.1880;  % Joint MAP (Green)
                  0.8500 0.3250 0.0980;  % ML (Orange)
                  0 0.4470 0.7410];      % EKF (Blue)

    % Line plot (keep order as is if you want legend to match original)
    figure('Position', [100, 100, 800, 380]);
    loglog(sweeplist, ekf_times, 'o-', 'LineWidth', 2, 'DisplayName', 'Wu et al.', 'Color', bar_colors(3,:));
    hold on;
    loglog(sweeplist, jointmap_times, 's-', 'LineWidth', 2, 'DisplayName', 'Proposed', 'Color', bar_colors(1,:));
    loglog(sweeplist, ML_times, 'x-', 'LineWidth', 2, 'DisplayName', 'Kok et al. (20 cores)', 'Color', bar_colors(2,:));
    xlabel(xlabel_str);
    ylabel('Computation Time (s)');
    legend('Location', 'best');
    grid on;
    ax = gca;
    ax.FontSize = 14;
    set(gca, 'LooseInset', max(get(gca,'TightInset'), 0.02)); 
    % Magnetometer Bias Error
    figure('Position', [100, 100, 400, 400]);
    data(:,1) = [statMAP.mag_bias_err];   % Proposed
    data(:,2) = [statML.mag_bias_err];    % Kok et al.
    data(:,3) = [statEKF.mag_bias_err];   % Wu et al.
    b = bar(data, 'grouped');

    for k = 1:length(b)
        b(k).FaceColor = bar_colors(k,:);
    end

    xlabel(xlabel_str);
    ylabel('Magnetometer Bias Error ($\mu$T)', 'Interpreter', 'latex');
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.XTick = 1:length(sweeplist);
    ax.XTickLabel = string(sweeplist);
    ylim([0.5*min(data(:)), 2*max(data(:))]); 
    set(gca, 'YScale', 'log');

    % D^m Error
    figure('Position', [100, 100, 400, 400]);
    data(:,1) = [statMAP.D_err];     % Proposed
    data(:,2) = [statML.D_err];      % Kok et al.
    data(:,3) = [statEKF.D_err];     % Wu et al.
    b = bar(data, 'grouped');

    for k = 1:length(b)
        b(k).FaceColor = bar_colors(k,:);
    end

    xlabel(xlabel_str);
    ylabel('$D^m$ Error', 'Interpreter', 'latex');
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.XTick = 1:length(sweeplist);
    ax.XTickLabel = string(sweeplist);
    ylim([0.5*min(data(:)), 2*max(data(:))]); 
    set(gca, 'YScale', 'log');

    % Accelerometer Bias Error
    figure('Position', [100, 100, 400, 400]);
    data(:,1) = [statMAP.acc_bias_err];  % Proposed
    data(:,2) = [statML.acc_bias_err];   % Kok et al.
    data(:,3) = [statEKF.acc_bias_err];  % Wu et al.
    b = bar(data, 'grouped');

    for k = 1:length(b)
        b(k).FaceColor = bar_colors(k,:);
    end

    xlabel(xlabel_str);
    ylabel('Accelerometer Bias Error ($\mathrm{m/s}^2$)', 'Interpreter', 'latex');
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.XTick = 1:length(sweeplist);
    ax.XTickLabel = string(sweeplist);
    ylim([0.1*min(data(:)), 1.2*max(data(:))]); 
    set(gca, 'YScale', 'log');

    % Gyroscope Bias Error
    figure('Position', [100, 100, 400, 400]);
    data(:,1) = rad2deg([statMAP.gyro_bias_err]);  % Proposed
    data(:,2) = rad2deg([statML.gyro_bias_err]);   % Kok et al.
    data(:,3) = rad2deg([statEKF.gyro_bias_err]);  % Wu et al.
    b = bar(data, 'grouped');

    for k = 1:length(b)
        b(k).FaceColor = bar_colors(k,:);
    end

    xlabel(xlabel_str);
    ylabel('Gyroscope Bias Error ($^\circ$/s)', 'Interpreter', 'latex');
    grid on;
    ax = gca;
    ax.FontSize = 14;
    ax.XTick = 1:length(sweeplist);
    ax.XTickLabel = string(sweeplist);
    ylim([0.8*min(data(:)), 2*max(data(:))]); 
    set(gca, 'YScale', 'log');
end
