%%%%%%%%%%% CTG HRV Analysis 
% Converts FHR (BPM) to RR intervals sample-by-sample: RR = 60000 / FHR
% Segments: non-overlapping 1-minute windows. Will edit later for N = # of
% beats for segments later. Dr. Pini said 1-m is fine for now. 

%% File Selection
[file_name, file_path] = uigetfile('*.mat', 'Select the CTG .mat file');

if isequal(file_name, 0)
    disp('User selected Cancel');
    return;
else
    full_path = fullfile(file_path, file_name);
    load(full_path);
    fprintf('Loaded: %s\n', file_name)                                                                                                  
                                                                        
    %%%%%%%%%%% Parameters: Use same as Original Change Fs. 
    %%%%%%%%%%% Dopp. Ultrasound: Non-invasive Fetal Reading
    fs          = 4;              % CTG sampling frequency (Hz)
    perc_change = 0.25;           % Max allowed fractional beat-to-beat change
    n_beats     = 5;              % Window size for outlier detection
    std_x       = 2;              % Std multiplier for outlier 
    seg_len_s   = 60;             % Segment length in seconds (1 minute)
    seg_len_smp = seg_len_s * fs; % Segment length in samples (= 240)

    [~, file_base, ~] = fileparts(file_name);

    %% Split signal into non-overlapping 1-minute segments
    fhr_signal = signal_interp_dec(:);  % ensure column vector
    n_total    = length(fhr_signal);
    n_segs     = floor(n_total / seg_len_smp);

    fprintf('Total samples: %d | Segment length: %d samples | Segments: %d\n', ...
        n_total, seg_len_smp, n_segs);

    % Pre-allocate output arrays
    fRR_N_initial = nan(1, n_segs);
    fRR_N_final   = nan(1, n_segs); 
    fHR_final     = nan(1, n_segs);
    fSD           = nan(1, n_segs);
    fRMSSD        = nan(1, n_segs);
    fe            = nan(1, n_segs);

    %%%%%%%% Main loop over segments
    for pz = 1:n_segs

        %%%%%%%%%% Extract segment
        seg_start = (pz - 1) * seg_len_smp + 1;                         
        seg_end   =  pz      * seg_len_smp;
        fhr_seg   = fhr_signal(seg_start:seg_end);  % FHR in BPM, 240 samples

        %%%%%%%%%% Limits from Studies 
        fhr_seg(fhr_seg < 50 | fhr_seg > 200) = NaN;

        %%%%%%%%%% Skip segment if too few (valid) samples
        if sum(~isnan(fhr_seg)) < 5
            fprintf('Segment %d: too few valid FHR samples, skipping.\n', pz);
            continue
        end

        %%%%%%%%% Convert FHR to RR intervals (ms), sample-by-sample
        fRR = 60000 ./ fhr_seg;  % RR in ms, one value per 4 Hz sample

        fRR_N_initial(pz) = sum(~isnan(fRR));

        %%%%%%%%% Outlier cleaning on RR series (no interp)

        %%%%%%%%% Rolling-window mean ± std_x*SD bounds
        fRR_mean      = nan(size(fRR));
        fRR_std_upper = nan(size(fRR));
        fRR_std_lower = nan(size(fRR));

        for s = 1:length(fRR)
            if s <= n_beats
                win = fRR(1:n_beats);
            else
                win = fRR(s-n_beats:s);
            end
            win_valid = win(~isnan(win));
            if length(win_valid) < 2, continue; end
            fRR_mean(s)      = mean(win_valid);
            fRR_std_upper(s) = mean(win_valid) + std_x * std(win_valid);
            fRR_std_lower(s) = mean(win_valid) - std_x * std(win_valid);
        end

        outlier_idx = fRR < fRR_std_lower | fRR > fRR_std_upper;
        fRR_clean   = fRR;
        fRR_clean(outlier_idx) = NaN;

        %%%%%%%%%%%%% Beat-to-beat Percentage Change filter
        change_sig     = (fRR_clean(2:end) - fRR_clean(1:end-1)) ./ fRR_clean(1:end-1);
        change_sig_idx = change_sig < -perc_change | change_sig > perc_change;
        change_sig_idx = logical([1; change_sig_idx]);
        fRR_clean(change_sig_idx) = NaN;

        %%%%%%%%%%%%% Global 3-SD filter
        g_mean = mean(fRR_clean, 'omitmissing');
        g_std  = std(fRR_clean,  'omitmissing');
        fRR_clean(fRR_clean < g_mean - 3*g_std | ...
                  fRR_clean > g_mean + 3*g_std) = NaN;

        fRR_N_final(pz) = sum(~isnan(fRR_clean));

        if fRR_N_final(pz) < 5
            fprintf('Segment %d: too few clean RR values (%d), metrics set to NaN.\n', pz, fRR_N_final(pz));
            continue
        end

        %%%%%%%%%%% HRV metrics 
        fHR_final(pz) = 60000 / mean(fRR_clean, 'omitmissing');  % mean HR [BPM]
        fSD(pz)       = std(fRR_clean, 'omitmissing');            % SDNN [ms]

        diff_fRR = diff(fRR_clean);
        diff_fRR(isnan(diff_fRR)) = [];
        if ~isempty(diff_fRR)
            fRMSSD(pz) = sqrt(mean(diff_fRR .^ 2));              % RMSSD [ms] Note: Ask about low values here & in SampEn.
        end

        fe(pz) = sampen(fRR_clean, 2, 0.2);                      % Sample Entropy- Stick to default tolerance 2/ All 2's. Articles say values around 0.9-1 are good. But IDK. Check big dawg.

        %%%%%%%%%% Plotting
        t_seg = (0:seg_len_smp-1)' / fs / 60;   % Time Axis in Minutes

        figure('WindowState', 'maximized');

        %%%%%%%%% Panel 1: Raw FHR signal
        subplot(2,2,1)
        plot(t_seg, fhr_seg, 'b', 'LineWidth', 1.5);
        xlabel('Time [min]')
        ylabel('FHR [BPM]')
        title('Raw FHR segment')
        ylim([50 200])

        %%%%%%%% Panel 2: Cleaned RR series over time
        subplot(2,2,3)
        plot(t_seg, fRR_clean, 'k', 'LineWidth', 1.5);
        xlabel('Time [min]')
        ylabel('RR interval [ms]')
        title('Cleaned RR series')

        %%%%%%%% Panel 3: Histogram of cleaned RR intervals
        subplot(2,2,2)
        histogram(fRR_clean, 'BinLimits', [250 750], 'BinWidth', 5);
        xlabel('RR interval [ms]')
        ylabel('Count')
        title(sprintf('Segment %d / %d', pz, n_segs))

        %%%%%%%% Metrics text box in histogram panel
        metrics_str = sprintf('N detected: %d\nN cleaned: %d\nfHR: %.1f BPM\nfSDNN: %.1f ms\nfRMSSD: %.1f ms\nfSampEn: %.3f', ...
            fRR_N_initial(pz), fRR_N_final(pz), fHR_final(pz), ...
            fSD(pz), fRMSSD(pz), fe(pz));
        ax = gca;
        text(ax, 260, ax.YLim(2)*0.95, metrics_str, ...
            'FontSize', 10, 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'black', 'Color', 'white');

        %%%%%%%%% Panel 4: Poincaré plot
        subplot(2,2,4)
        scatter(fRR_clean(1:end-1), fRR_clean(2:end), 'blue');
        xlim([250 750]); ylim([250 750]);
        xlabel('RR_n [ms]')
        ylabel('RR_{n+1} [ms]')
        title('Poincaré plot')

        sgtitle(sprintf('CTG segment %d / %d   (%.1f – %.1f min)', ...
            pz, n_segs, (pz-1)*seg_len_s/60, pz*seg_len_s/60))

        pause
        close all
        clear fhr_seg fRR fRR_clean fRR_mean fRR_std_lower fRR_std_upper

    end

    %%%%%%%%%% Save results to Excel: Fetal 
    segment_N   = (1:n_segs)';
    f_final_mat = table(segment_N, fRR_N_initial', fRR_N_final', ...
                        fHR_final', fSD', fRMSSD', fe', ...
        'VariableNames', {'Segment_idx', 'N_detected_RR', 'N_cleaned_RR', ...
                          'fHR_BPM', 'fSDNN_ms', 'fRMSSD_ms', 'fSampEn'});

    out_name = [file_base '_CTG_HRV.xlsx'];
    writetable(f_final_mat, out_name, 'WriteRowNames', false);
    fprintf('Results saved to: %s\n', out_name);

end