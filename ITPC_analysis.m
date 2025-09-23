%%%% ITPC analysis
%%% Bernardo AO


function [band_itpc, band_itpc_t, n_trials] = ITPC_analysis(lfps, ...
    align_indices, time_window, opt)
    % This script runs an ITPC on the lfp alligned to align_indx.
    % The analysis it performs are: 
    % Calculate phase: Plots the phase of the lfp as a function of time
    % Calculate ITPC: Plots the phase of each trial for the align indx
    % ITPC with respect to time: Plots the ITPC for a time window centered 
    % around the align indx

    % Inputs:
    % lfps: cell of 1d arrays with the LFP signal for each session
    % align_indx: cell of array with the indices to align the lfp 
    % time_window: array with the desired time window around the align 
    % times to study [s]

    % opt | extra parametes:
    % fs: sampling time_vec
    % step: size of the step for ITPC for the time window
    % freq_band: desired freq_band [start, end]
    % n_cycles: number of cycles for the morlet wavelet
    % band_name: name of the band for the plots
    % align_name: name of the alignment
    % color: color for the plots
    % show_fig: logical

    % Outputs:
    % band_itpc: ITPC at the align time
    % band_itpc_t: ITPC around the align time given the time window
    % n_trials: totla number of trials detected


    save_path = "W:\Lorena\Analysis_scripts\Bernardo_code\plots";
    ext = ".png";
    n_rand = 1000;

    time_window_ind = int32(time_window(1)*opt.fs: opt.fs*opt.step : time_window(2)*opt.fs);
    time_window_vec = time_window(1):opt.step:time_window(2);
    len_t = length(time_window_ind);

    n_sessions = numel(lfps);
    n_trials = calc_n_trials(align_indices);
    thetas_all = zeros(n_trials, 1); 
    thetas_all_t = zeros(n_trials, len_t); 
    thetas_rand = zeros(n_trials, n_rand); 
    
    trial = 1;
    for n = 1:n_sessions
        lfp = lfps{n};
        align_indx = int32(align_indices{n});
        n_t = length(align_indx);
        ct = trial:trial + n_t - 1;

        % Calculate phase
        
        %time_vec = (0:length(lfp)-1) / opt.fs;
        band_phase = get_phase_LFP(lfp, opt.fs, opt.freq_band, opt.n_cycles);
        thetas = band_phase(align_indx);
        thetas_all(ct) = thetas;
        
        % Calculate significance
        parfor nr = 1:n_rand
            align_indx_r = mod(align_indx + randi(length(band_phase)), ...
                length(band_phase));
            thetas_rand(ct,nr) = band_phase(align_indx_r);
        end

        % ITPC with respect to time        
        ti = 0;
        for t = time_window_ind
            ti = ti + 1;
            thetas_all_t(ct, ti) = band_phase(align_indx + t);
        end
        
        trial = trial + n_t;
    end

    % Calculate ITPC
    band_itpc = ITPC(thetas_all);

    itpc_rand = zeros(n_rand,1);
    for nr = 1:n_rand
        itpc_rand(nr) = ITPC(thetas_rand(:,nr));
    end
    mu_r = mean(itpc_rand); 
    std_r = std(itpc_rand);

    band_itpc_t = zeros(len_t, 1);
    for t = 1:len_t
        band_itpc_t(t) = ITPC(thetas_all_t(:, t));
    end

    %% Plottig

    % phase plot
    %{
    figure;
        ax = gca();
        yyaxis left
        plot(time_vec, lfp, Color="red")
        ylabel('LFP');
        yyaxis right
        plot(time_vec, band_phase, LineStyle="--",Color="k");
        xlabel('Time (s)');
        ylabel(opt.band_name + " phase (radians)");
        ax.YAxis(1).Color = "red";
        ax.YAxis(2).Color = "k";
        title(opt.band_name + " phase vs Time");
    %}

    % angle plot
    
    if opt.show_fig
        figure('Name','Polar' + opt.band_name);
    else
        figure('Visible', 'off');
    end
    t_name = opt.band_name + " phase at " + opt.align_name;
    polarhistogram(thetas_all, 10,Normalization="probability",FaceColor=opt.color)
    ax = gca;          % Get current polar axes
    r_limits = ax.RLim(2);
    text(pi/4,  1.1*r_limits, "ITPC = " + num2str(band_itpc) + newline + "n = " + n_trials)
    hold on;  

    for k = 1:length(thetas_all)
        x = [thetas_all(k) thetas_all(k)];
        y = [0 0.1];
        
        polarplot(x, y, 'Color', "#808080", 'LineWidth', 1);
        rticklabels("")
    end
    hold off;
    
    title(t_name)
    saveas(gcf, fullfile(save_path, t_name + ext));

    % ITPC with time plot
    if opt.show_fig
        figure('Name','time' + opt.band_name);
    else
        figure('Visible', 'off');
    end

    plot(time_window_vec, band_itpc_t, color=opt.color)
    hold on
    fill([time_window_vec(1), time_window_vec(end), ...
        time_window_vec(end), time_window_vec(1)], ...
        [0, 0, mu_r + std_r, mu_r + std_r], ...
        [0.2,0.22,0.25], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

    xline(0,LineStyle=":",Color="#343a40")

    xlabel("time [s]")
    ylabel("ITPC" + ", n = " + n_trials)
    ylim([0,0.5])
    box off
    t_name = opt.band_name + " ITPC around " + opt.align_name;
    title(t_name)
    saveas(gcf, fullfile(save_path, t_name + ext));

    %% Normalize for output
    band_itpc = (band_itpc - mu_r) / std_r;
    band_itpc_t = (band_itpc_t - mu_r) / std_r;
end



% Helper functions

function n_trials = calc_n_trials(cell)
    % Inputs:
    % cell: array with each session, and trials
    % 
    % Outputs:
    % n_trials: total number of trials for each session
    n_trials = 0;
    for session = 1:numel(cell)
        n_trials = n_trials + length(cell{session});
    end
end

function phase_LFP = get_phase_LFP(lfp, fs, freq_band, n_cycles)
    % Inputs:
    % lfp: 1d array with the LFP signal
    % fs: sampling time_vec
    % freq_band: desired freq_band [start, end]
    % n_cycles: number of cycles for the morlet wavelet
    % 
    % Outputs:
    % phase_LFP: phase of the LFP for each time bin
    
    % Morlet wavelet
    center_freq = round(mean(freq_band));              

    w_len = 6 / center_freq; 
    t = -w_len/2 : 1/fs : w_len/2;
    s = n_cycles / (2*pi*center_freq); % std of Gaussian
    morlet_wavelet = exp(2*1i*pi*center_freq*t) .* exp(-t.^2/(2*s^2));
    
    % Convolve LFP with Morlet wavelet and extract phase
    conv_result = conv(lfp, morlet_wavelet, 'same');
    phase_LFP = angle(conv_result);    

end

function itpc = ITPC(phase)
    % Inputs:
    % phase: an array shape (trials)
    %
    % Outputs:
    % itpc: inter trial phase clustering metric

    itpc = abs(mean(exp(phase*1i)));

end

