%% Td analysis functions
% Bernardo AO

function [vel_grid, ratio_grid] = theta_delta_ratio(pideals, lfps, opt)
    % Inputs:
    % pideals: cell of 1d arrays with the pideals object for each session
    % lfps: cell of 1d arrays with the LFP object for each session
    % % opt | extra parametes ( Refer to ITPC_analysis)

    % Outputs:
    % vel_grid: 2d array with the mean velocity for each grid point
    % ratio_grid: 2d array with the mean theta delta ratio
    %   for each grid point
    
    n_sess = numel(lfps);
    vel_vec = [];
    ratio_vec = [];
    x_vec = [];
    y_vec = [];

    for n = 1:n_sess
        lfp = lfps{n};
        pideal = pideals{n};

        log_ratio_pt = TD_speed(pideal, lfp, opt.n_cycles, opt.win);
        vel_vec = [vel_vec; pideal.v];
        ratio_vec = [ratio_vec; log_ratio_pt];

        x_vec = [x_vec; pideal.x];
        y_vec = [y_vec; pideal.y];
    end

    vel_grid = plot_heatmap(x_vec, y_vec, vel_vec, "velocity", opt);
    ratio_grid = plot_heatmap(x_vec, y_vec, ratio_vec, "ratio", opt);
    %plot_speed_td(vel_vec,ratio_vec, opt)
end

function conv_result = conv_LFP(lfp, fs, freq_band, n_cycles)
    % Inputs:
    % lfp: 1d array with the LFP signal
    % fs: sampling time_vec
    % freq_band: desired freq_band [start, end]
    % n_cycles: number of cycles for the morlet wavelet
    % 
    % Outputs:
    % conv_result: convolution of the LFP with a morlet wavelett
    
    % Morlet wavelet
    center_freq = round(mean(freq_band));              
    
    w_len = n_cycles / center_freq; 
    t = -w_len/2 : 1/fs : w_len/2;
    s = n_cycles / (2*pi*center_freq); % std of Gaussian
    morlet_wavelet = exp(2*1i*pi*center_freq*t) .* exp(-t.^2/(2*s^2));
    morlet_wavelet = morlet_wavelet / sqrt(sum(abs(morlet_wavelet).^2));
    
    % Time domain
    %conv_result = conv(lfp, morlet_wavelet, 'same');

    % FFT
    delay = floor(length(morlet_wavelet) / 2);
    conv_result = circshift(fftfilt(morlet_wavelet, lfp), -delay);
end

function log_ratio_pt = TD_speed(pideal, lfp, n_cycles, win)
    % Inputs:
    % pideal: pideal object 
    % lfp: LFP object
    % n_cycles: number of cycles for the morlet wavelet
    % win: window for each side to smooth the result

    % Outputs:
    % log_ratio_pt: array with the value of the theta delta
    %   ratio for each pideal value

    theta_band = [6, 12];
    delta_band = [1, 6];

    lfpd = Data(lfp);
    lfpt = EEG.toSecond(lfp, 'Range');
    pt = pideal.t;
    fs = EEG.Fs(lfp);

    tp = abs(conv_LFP(lfpd, fs, theta_band, n_cycles)).^2;
    dp = abs(conv_LFP(lfpd, fs, delta_band, n_cycles)).^2;

    log_ratio = log(tp ./ dp);
    

    log_ratio_pt = zeros(size(pt));
    t = 1;
    T = length(lfpt);
    for i = 1:length(pt)
        while t < T && abs(lfpt(t + 1) - pt(i)) < abs(lfpt(t) - pt(i))
            t = t + 1;
        end
        if t-win > 0
            start_ind = t-win;
        else
            start_ind = 1;
        end
        if t+win < T
            end_ind = t+win;
        else
            end_ind = T;
        end
        log_ratio_pt(i) = mean(log_ratio(start_ind:end_ind));
    end
end

function plot_speed_td(v,r,opt)
    % Plots a scatter and gets the correlation of a velocity vector and a
    % ratio vector (v,r)
    downsample = opt.downsample;
    name = opt.fig_name;
    color = opt.color;

    % Calculate correlation
    corr_coef = corr(v, r);
    
    p = polyfit(v, r, 1);  
    r_fit = polyval(p, v);    

    % Plot
    figure;
    d_idx = randi(length(v),1,downsample);
    h = scatter(v(d_idx),r(d_idx),"filled");
    h.MarkerFaceAlpha = 0.2;

    hold on;
    plot(v, r_fit, 'r-', 'LineWidth', 2, 'Color',color);  
    
    text(10,  8, "corr = " + num2str(corr_coef,2))
    xlabel("velocity")
    ylabel("log theta/delta")
    ylim([-10, 10])
    xlim([0,15])

    t_name = name + " TD ratio velocity Correlation ";
    title(t_name)
    saveas(gcf, fullfile(opt.save_path, t_name + opt.ext));

end

function z_grid = plot_heatmap(x_vec, y_vec, z_vec, name, opt)
    % Plots a heatmap given x,y, and z values

    x_edges = -23:23;
    y_edges = -35:35;
    x_bins = length(x_edges);
    y_bins = length(y_edges);

    % Digitize x and y vectors into bin indices
    [~, ~, x_idx] = histcounts(x_vec, x_edges);
    [~, ~, y_idx] = histcounts(y_vec, y_edges);

    % Remove data outside the bins
    valid = x_idx > 0 & y_idx > 0;
    x_idx = x_idx(valid);
    y_idx = y_idx(valid);
    z_vec = z_vec(valid);

    % Convert 2D indices to linear indices for accumarray
    linear_idx = sub2ind([y_bins, x_bins], y_idx, x_idx);

    % Compute sum and count of z_vec in each bin
    vel_sum = accumarray(linear_idx, z_vec, [y_bins*x_bins, 1], @sum, NaN);
    vel_count = accumarray(linear_idx, 1, [y_bins*x_bins, 1], @sum, NaN);

    % Compute average 
    vel_avg = vel_sum ./ vel_count;

    % Reshape into 2D grid
    z_grid = reshape(vel_avg, [y_bins, x_bins]);

    % Create x and y centers for plotting
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

    % Plot
    figure;
    imagesc(x_centers, y_centers, z_grid);
    axis xy; % flip y-axis so low values are at the bottom
    t_name = opt.fig_name + " Average " + name + " heatmap";
    title(t_name);
    axis off
    colormap(flipud(hot))
    a = colorbar;
    a.Label.String = name;

    if name == "theta power" || name == "delta power"
        m = mean(z_grid,"all","omitmissing");
        s = std(z_grid, [],"all","omitmissing");
        clim([0,m+3*s]);
    end
    saveas(gcf, fullfile(opt.save_path, t_name + opt.ext));
end

% Old code
%{

 lfp_raw = zeros(n_trials, len_t);
            lfp_raw(ct, ti) = lfp(align_indx + t);
        r = log(theta_power(trial,:) ./ delta_power(trial,:));
        figure();
        yyaxis left
        plot(time_window_vec, r,color="k")
        yyaxis right
        plot(time_window_vec, theta_power(trial,:),color="blue")
        hold on
        plot(time_window_vec, delta_power(trial,:),color="red")
        hold off

first_t = {};
    s = 1;
    for trial = 1:n_trials
        basal_r = min(ratio(trial,time_window_vec < 0));
        t_small_td = time_window_vec(ratio(trial,:) < basal_r);
        if ~isempty(t_small_td)
            first_t{s} = t_small_td(1);
            s = s + 1;
        end
    end

function [ratio_mean, min_time, max_time] = theta_delta_ratio(lfps, ...
    align_indices, time_window, opt)
    % Inputs:
    % Same as ITPC_analysis

    save_path = "W:\Lorena\Analysis_scripts\Bernardo_code\plots";
    ext = ".png";
    plot_ts = "max";
    theta_band = [6, 12];
    delta_band = [1, 6];

    time_window_ind = int32(time_window(1)*opt.fs: opt.fs*opt.step : time_window(2)*opt.fs);
    time_window_vec = time_window(1):opt.step:time_window(2);

    len_t = length(time_window_ind);
    n_sessions = numel(lfps);
    n_trials = calc_n_trials(align_indices);

    theta_power = zeros(n_trials, len_t); 
    delta_power = zeros(n_trials, len_t); 
    

    %% delta theta power
    trial = 1;
    for n = 1:n_sessions
        lfp = lfps{n};
        align_indx = int32(align_indices{n});
        n_t = length(align_indx);
        ct = trial:trial + n_t - 1;

        
        tp = abs(conv_LFP(lfp, opt.fs, theta_band, 10)).^2;
        dp = abs(conv_LFP(lfp, opt.fs, delta_band, 10)).^2;

        ti = 0;
        for t = time_window_ind
            ti = ti + 1;
            theta_power(ct, ti) = tp(align_indx + t);
            delta_power(ct, ti) = dp(align_indx + t);
        end

        trial = trial + n_t;
    end
    t_mean = mean(theta_power, 1,"omitmissing");
    d_mean = mean(delta_power, 1,"omitmissing");
    ratio_mean = log(t_mean ./ d_mean);

    % find the time of the min ratio
    ratio = log(theta_power ./ delta_power);
    [~,min_ind] = min(ratio,[],2);
    min_time = time_window_vec(min_ind);
    
    [~,max_ind] = max(delta_power,[],2);
    max_time = time_window_vec(max_ind);

    %% Plot
    figure('Name','Theta delta ratio' + opt.band_name);
    subplot(2,1,1)
    p1 = plot(time_window_vec, t_mean, color="blue");
    hold on
    p2 = plot(time_window_vec, d_mean, color="red");
    xline(0,LineStyle=":",Color="#343a40")
    legend([p1 p2], {"theta","delta"})
    xlim(time_window)
    ylim([0,0.5e-4])
    hold off    
    box off
    ylabel("Power")

    subplot(2,1,2)
    plot(time_window_vec, ratio_mean, color=opt.color)
    hold on
    xline(0,LineStyle=":",Color="#343a40")
    if plot_ts == "max"
        ts = max_time;
    else
        ts = min_time;
    end
    for t = ts 
        xline(t,LineStyle=":",Color="green")
    end
    hold off
    xlabel("time [s]")
    xlim(time_window)
    %ylim([-2 2])
    ylabel("log ratio")
    box off
    t_name = opt.band_name + ' delta ratio ' + opt.align_name;
    
    sgtitle(t_name)
    saveas(gcf, fullfile(save_path, t_name + ext));
end

% plot session
figure;
yyaxis left
plot(log_ratio_pt, Color="blue")
hold on
yyaxis right
plot(pideal.v,color="red")
hold off
legend("log_ratio","v")

figure;
h = scatter(pideal.v,log_ratio_pt, Color=opt.color);
h.MarkerFaceAlpha = 0.4;
ylim([0,1e-4])

error("done")

%}