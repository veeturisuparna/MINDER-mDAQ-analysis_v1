% =========================================================================
% main_EDA_Analysis.m
% =========================================================================
%
% Main script to run the complete Electrodermal Activity (EDA) analysis pipeline.
%
% Workflow:
% 1. Sets configuration parameters (participant, file paths, analysis settings).
% 2. Loads raw data from a DC-EDA device and a BIOPAC system.
% 3. Calls `processEDAData` to clean, resample, filter, and decompose signals, plotting each step.
% 4. Calls analysis functions to compare the two signals (SCL and SCR).
% 5. Calls plotting functions to visualize the final results.
%
% To run:
% - Set the `config.participant` variable.
% - Ensure data files are in the specified `config.dataPath`.
% - Run this script from the MATLAB command window or editor.
%
function main_EDA_Analysis()
    % close all; clear; clc; % Uncomment to clear workspace and close figures
    %% --- 1. Configuration ---
    % General Settings
    config.participant = 'P2'; % <-- Change to P2, P3, etc. as needed
    config.dataPath = '.'; % Assumes data files are in the same directory as the script
    % Signal Processing Settings
    config.resampleRateHz = 100;     % Common sampling rate for both signals
    config.filterCutoffHz.low = 1.5; % Low-pass filter for general smoothing
    config.filterCutoffHz.high = 10; % Low-pass filter before normalization
    config.sclCutoffHz = 0.05;        % Cutoff frequency to separate SCL (tonic) from SCR (phasic)
    
    % --- Artifact Removal Settings (MODIFIED FOR Z-SCORE NORMALIZATION) ---
    config.artifact.sd_window_sec = 60;
    % Threshold is now for Z-SCORED data. A value of 0.01 means "flag a window 
    % if its std dev is less than 1% of the overall signal's std dev".
    config.artifact.sd_threshold = 0.01; 
    config.artifact.slope_window_sec = 5;
    config.artifact.slope_variation = 0.20;
    config.artifact.absolute_slope_tolerance = 1e-4;

    % Analysis Settings
    config.analysis.dtw_window_min = 0.1;      % Sliding window duration for DTW (minutes)
    config.analysis.dtw_overlap_min = 0;     % Overlap for DTW windows (minutes)
    config.analysis.dtw_downsample_factor = 10; % Downsample factor to speed up DTW
    config.analysis.coherence_freq_max = 0.05; % Max frequency for mean coherence calculation
    config.analysis.scr_cross_window = 0.5;   % 1 second window for moving-average baseline
    config.analysis.scr_cross_std_factor = 0.1; % Threshold is this factor times the signal's std dev
    config.analysis.scr_min_zc_dist_sec = 0.5; % Minimum distance between ZCR events (seconds)
    
    % Event Definitions (in minutes)
    config.events = {
        'DB',   [5, 10];
        'MM1',  [14, 16];
        'MM2',  [17, 19.5];
        'MM3',  [20.5, 25];
    };
    %% --- 2. Load Data ---
    fprintf('--- Loading data for participant: %s ---\n', config.participant);
    
    % Load custom DC-EDA data
    dcFile = fullfile(config.dataPath, sprintf('%s.csv', config.participant));
    if ~exist(dcFile, 'file')
        error('DC-EDA file not found: %s', dcFile);
    end
    EDA_In = readtable(dcFile);
    % Load BIOPAC data from .mat file
    biopacFile = fullfile(config.dataPath, sprintf('%s (BIOPAC_only_All_Sensors_Ch1Sync).mat', config.participant));
    if ~exist(biopacFile, 'file')
        error('BIOPAC file not found: %s', biopacFile);
    end
    
    loaded_data = load(biopacFile);
    field_names = fieldnames(loaded_data);
    if isempty(field_names)
        error('BIOPAC .mat file appears to be empty or in an unexpected format: %s', biopacFile);
    end
    BIOPAC_raw = loaded_data.(field_names{1});
    fprintf('  - Loaded variable ''%s'' from BIOPAC .mat file.\n', field_names{1});
    
    %% --- 3. Process Signals ---
    fprintf('--- Processing EDA signals ---\n');
    [dc_eda, biopac_eda] = processEDAData(EDA_In, BIOPAC_raw, config);
    %% --- 4. Perform Analysis ---
    fprintf('--- Analyzing SCL and SCR components ---\n');
    scl_results = analyzeSCL(dc_eda, biopac_eda, config.analysis);
    scr_results = analyzeSCR(dc_eda, biopac_eda, config); % Pass full config
    event_results = analyzeEvents(dc_eda, biopac_eda, config.events);
    %% --- 5. Plot Results ---
    fprintf('--- Generating final plots ---\n');
    plotSignalOverview(dc_eda, biopac_eda, scl_results, scr_results);
    plotEventAnalysis(event_results, config.events);
    plotDecomposition(dc_eda, biopac_eda);
    
    fprintf('--- Analysis complete for %s ---\n', config.participant);
end

% =========================================================================
% processEDAData.m (UPDATED FUNCTION)
% =========================================================================
function [dc_eda, biopac_eda] = processEDAData(EDA_In, BIOPAC_raw, config)
    % --- Prepare Initial Signals ---
    % Process DC-EDA data based on the columns found in the input table
    if any(strcmp(EDA_In.Properties.VariableNames, 'VOLTAGE_MV'))
        fprintf('  - Processing DC-EDA data from VOLTAGE_MV and ISI_US columns.\n');
        isi_us = EDA_In.ISI_US;
        
        voltage_mv = EDA_In.VOLTAGE_MV;
        DC_EDA_intermediate = fillmissing(voltage_mv / 100, 'linear');
        signal_DCEDA_raw = 3300./((3.3-(2*DC_EDA_intermediate/1000))*825);
        signal_DCEDA_raw = fillmissing(signal_DCEDA_raw, 'linear'); % Ensure no NaNs from conversion
        
        time_DCEDA_raw = cumsum(isi_us) / (1e6 * 60); % Convert microseconds -> minutes
    elseif any(strcmp(EDA_In.Properties.VariableNames, 'eda'))
        fprintf('  - Processing DC-EDA data from eda and isi/ISI_US columns.\n');
        if any(strcmp(EDA_In.Properties.VariableNames, 'isi'))
            isi_us = EDA_In.isi;
        elseif any(strcmp(EDA_In.Properties.VariableNames, 'ISI_US'))
            isi_us = EDA_In.ISI_US;
        else
            error('Could not find ISI column (expected "isi" or "ISI_US") in DC-EDA file.');
        end
        signal_DCEDA_raw = fillmissing(EDA_In.eda, 'linear');
        time_DCEDA_raw = cumsum(isi_us) / (1e6 * 60); % Convert microseconds -> minutes
    elseif any(strcmp(EDA_In.Properties.VariableNames, 'EDA_uS'))
        fprintf('  - Processing DC-EDA data from EDA_uS and ISI_US columns (original format).\n');
        time_DCEDA_raw = cumsum(EDA_In.ISI_US) / (1e6 * 60); % microseconds -> minutes
        signal_DCEDA_raw = fillmissing(EDA_In.EDA_uS, 'linear');
    else
        error('Could not find expected EDA columns in the DC-EDA file. Looked for "VOLTAGE_MV", "eda", or "EDA_uS".');
    end
    
    % Extract BIOPAC signal and create time vector in minutes
    fs_biopac_original = 1000; % BIOPAC AcqKnowledge saves at 1kHz
    signal_BIOPAC_raw = fillmissing(BIOPAC_raw(:,5), 'linear'); % Using column 5 as per your script
    time_BIOPAC_raw = (0:length(signal_BIOPAC_raw)-1)' / (fs_biopac_original * 60);
    
    % --- Resample to a Common Time Base ---
    t_max = max(time_DCEDA_raw(end), time_BIOPAC_raw(end));
    fs_resample_min = config.resampleRateHz * 60; % Target sample rate in samples/min
    time_uniform = (0 : 1/fs_resample_min : t_max)';
    target_len = length(time_uniform);
    
    signal_DCEDA_resampled = resample(signal_DCEDA_raw, time_DCEDA_raw, fs_resample_min);
    signal_BIOPAC_resampled = resample(signal_BIOPAC_raw, time_BIOPAC_raw, fs_resample_min);
    
    % --- Pad or Truncate signals to ensure they have the same length ---
    fprintf('  - Standardizing signal lengths to %d points.\n', target_len);
    len_dc = length(signal_DCEDA_resampled);
    if len_dc > target_len
        signal_DCEDA_resampled = signal_DCEDA_resampled(1:target_len);
    elseif len_dc < target_len
        padding = repmat(signal_DCEDA_resampled(end), target_len - len_dc, 1);
        signal_DCEDA_resampled = [signal_DCEDA_resampled; padding];
    end
    len_bio = length(signal_BIOPAC_resampled);
    if len_bio > target_len
        signal_BIOPAC_resampled = signal_BIOPAC_resampled(1:target_len);
    elseif len_bio < target_len
        padding = repmat(signal_BIOPAC_resampled(end), target_len - len_bio, 1);
        signal_BIOPAC_resampled = [signal_BIOPAC_resampled; padding];
    end
    
    % --- Store time vector and raw resampled signals for plots ---
    dc_eda.time = time_uniform;
    biopac_eda.time = time_uniform;
    dc_eda.raw = signal_DCEDA_resampled;
    biopac_eda.raw = signal_BIOPAC_resampled;
    
    % --- NEW: NORMALIZE signals via z-score before artifact detection ---
    fprintf('  - Normalizing signals via z-score before artifact detection...\n');
    dc_raw_z = zscore(dc_eda.raw);
    biopac_raw_z = zscore(biopac_eda.raw);
    
    % --- Remove Artifacts (Applied to Z-Scored Data) ---
    fprintf('  - Applying artifact removal...\n');
    fs_hz = config.resampleRateHz;
    
    % Define a nested helper function for artifact detection and correction
    function cleaned_signal = detectAndCorrect(signal, artifact_config)
        is_artifact = movstd(signal, round(artifact_config.sd_window_sec * fs_hz)) < artifact_config.sd_threshold;
        slopes = gradient(signal) * fs_hz; % Slope in units per second
        is_artifact = is_artifact | (abs(slopes) > artifact_config.absolute_slope_tolerance);
        signal_with_nans = signal;
        signal_with_nans(is_artifact) = NaN;
        cleaned_signal = fillmissing(signal_with_nans, 'linear', 'EndValues', 'nearest');
        if any(isnan(cleaned_signal))
            cleaned_signal = fillmissing(cleaned_signal, 'nearest');
        end
    end
    
    % Apply artifact correction to the Z-SCORED signals
    dc_cleaned_z = detectAndCorrect(dc_raw_z, config.artifact);
    biopac_cleaned_z = detectAndCorrect(biopac_raw_z, config.artifact);
    
    % --- NEW: DE-NORMALIZE the cleaned signals back to their original scale ---
    fprintf('  - De-normalizing signals back to original scale...\n');
    % To de-normalize: new_value = (z_scored_value * original_std) + original_mean
    dc_eda.cleaned = (dc_cleaned_z * std(dc_eda.raw, 'omitnan')) + mean(dc_eda.raw, 'omitnan');
    biopac_eda.cleaned = (biopac_cleaned_z * std(biopac_eda.raw, 'omitnan')) + mean(biopac_eda.raw, 'omitnan');

    % --- Filtering ---
    [b_low, a_low] = butter(8, config.filterCutoffHz.low / (fs_hz/2), 'low');
    [b_high, a_high] = butter(4, config.filterCutoffHz.high / (fs_hz/2), 'high');
    dc_eda.filtered_1_5Hz = filtfilt(b_low, a_low, dc_eda.cleaned);
    dc_eda.filtered_10Hz = filtfilt(b_high, a_high, dc_eda.cleaned);
    
    biopac_eda.filtered_1_5Hz = filtfilt(b_low, a_low, biopac_eda.cleaned);
    biopac_eda.filtered_10Hz = filtfilt(b_high, a_high, biopac_eda.cleaned);
    
    % --- SCL/SCR Decomposition using FFT Masking ---
    [dc_eda.scl, dc_eda.scr] = decomposeEDA(dc_eda.filtered_1_5Hz, fs_hz, 0.05);
    [biopac_eda.scl, biopac_eda.scr] = decomposeEDA(biopac_eda.filtered_1_5Hz, fs_hz, 0.05);
    
    % --- Normalization for Event Analysis ---
    % Min-Max normalization on the 1.5Hz filtered signal
    dc_eda.normalized = (dc_eda.filtered_1_5Hz - min(dc_eda.filtered_1_5Hz)) / (max(dc_eda.filtered_1_5Hz) - min(dc_eda.filtered_1_5Hz));
    biopac_eda.normalized = (biopac_eda.filtered_1_5Hz - min(biopac_eda.filtered_1_5Hz)) / (max(biopac_eda.filtered_1_5Hz) - min(biopac_eda.filtered_1_5Hz));
end

% =========================================================================
% decomposeEDA.m (Sub-function)
% =========================================================================
function [scl, scr] = decomposeEDA(signal, fs, scl_cutoff_hz)
    signal_fft = fft(signal);
    n = length(signal);
    
    cutoff_index = round(scl_cutoff_hz / (fs / n));
    freq_mask = zeros(size(signal_fft));
    freq_mask(1:cutoff_index) = 1;
    freq_mask(end-cutoff_index+1:end) = 1; % Include negative frequencies
    
    scl_fft = signal_fft .* freq_mask;
    scl = real(ifft(scl_fft));
    scr = signal - scl;
end

% =========================================================================
% analyzeSCL.m
% =========================================================================
function results = analyzeSCL(dc_eda, biopac_eda, analysis_config)
    scl_dc_z = zscore(dc_eda.scl);
    scl_biopac_z = zscore(biopac_eda.scl);
    
    [xc_scl, lags_scl] = xcorr(scl_dc_z, scl_biopac_z, 'coeff');
    [results.max_corr, idx_max] = max(xc_scl);
    best_lag_samples = lags_scl(idx_max);
    
    fs_hz = length(biopac_eda.time) / (biopac_eda.time(end) * 60);
    fs_min = fs_hz * 60;
    results.lag_min = best_lag_samples / fs_min;
    
    scl_dc_shifted = circshift(dc_eda.scl, -best_lag_samples);
    
    [coh_scl, f_scl] = mscohere(dc_eda.scl, biopac_eda.scl, hamming(256), [], [], fs_hz);
    results.avg_coherence = mean(coh_scl(f_scl <= analysis_config.coherence_freq_max));
    
    window_len = round(analysis_config.dtw_window_min * fs_min);
    step = window_len - round(analysis_config.dtw_overlap_min * fs_min);
    ds_factor = analysis_config.dtw_downsample_factor;
    
    num_windows = floor((length(biopac_eda.time) - window_len) / step);
    if num_windows < 1
        warning('Signal is too short for the specified DTW window size.');
        results.avg_dtw = NaN;
    else
        normalized_dtw_vals = zeros(1, num_windows);
        for i = 1:num_windows
            idx_start = (i - 1) * step + 1;
            idx_end = idx_start + window_len - 1;
            
            seg1 = scl_dc_shifted(idx_start:idx_end);
            seg2 = biopac_eda.scl(idx_start:idx_end);
            
            seg1_norm = zscore(seg1);
            seg2_norm = zscore(seg2);
            
            if any(isnan(seg1_norm)), seg1_norm(:) = 0; end
            if any(isnan(seg2_norm)), seg2_norm(:) = 0; end
            
            seg1_ds = downsample(seg1_norm, ds_factor);
            seg2_ds = downsample(seg2_norm, ds_factor);
    
            dtw_dist = dtw(seg1_ds, seg2_ds);
            normalized_dtw_vals(i) = dtw_dist / length(seg1_ds);
        end
        results.avg_dtw = mean(normalized_dtw_vals, 'omitnan');
    end
    
    results.display_str = sprintf('SCL Analysis | Cross-corr: %.2f | Lag: %.2f min | Mean Coh (f<%.2fHz): %.2f | Norm. DTW: %.2f', ...
                          results.max_corr, results.lag_min, analysis_config.coherence_freq_max, results.avg_coherence, results.avg_dtw);
    disp(results.display_str);
end

% =========================================================================
% analyzeSCR.m
% =========================================================================
function results = analyzeSCR(dc_eda, biopac_eda, config)
    fs = config.resampleRateHz;
    win_samples = fs * config.analysis.scr_cross_window;
    baseline_dc     = movmean(dc_eda.scr, win_samples);
    baseline_biopac = movmean(biopac_eda.scr, win_samples);
    
    scr_dc_centered     = dc_eda.scr - baseline_dc;
    scr_biopac_centered = biopac_eda.scr - baseline_biopac;
    
    thresh_dc     = std(scr_dc_centered) * 0.07;
    thresh_biopac = std(scr_biopac_centered) * 0.009;
    
    initial_idx_dc     = find(diff(sign(scr_dc_centered)) ~= 0 & ...
                              abs(scr_dc_centered(1:end-1)) > thresh_dc);
    initial_idx_biopac = find(diff(sign(scr_biopac_centered)) ~= 0 & ...
                              abs(scr_biopac_centered(1:end-1)) > thresh_biopac);
                      
    min_dist_samples = round(config.analysis.scr_min_zc_dist_sec * fs);
    if ~isempty(initial_idx_dc)
        final_idx_dc = initial_idx_dc(1);
        for i = 2:length(initial_idx_dc)
            if (initial_idx_dc(i) - final_idx_dc(end)) > min_dist_samples
                final_idx_dc = [final_idx_dc; initial_idx_dc(i)];
            end
        end
        idx_dc = final_idx_dc;
    else
        idx_dc = [];
    end
    if ~isempty(initial_idx_biopac)
        final_idx_biopac = initial_idx_biopac(1);
        for i = 2:length(initial_idx_biopac)
            if (initial_idx_biopac(i) - final_idx_biopac(end)) > min_dist_samples
                final_idx_biopac = [final_idx_biopac; initial_idx_biopac(i)];
            end
        end
        idx_biopac = final_idx_biopac;
    else
        idx_biopac = [];
    end
    
    results.num_dc_cross     = numel(idx_dc);
    results.num_biopac_cross = numel(idx_biopac);
    
    if results.num_biopac_cross > 0
        results.deltaZC = abs(results.num_dc_cross - results.num_biopac_cross) / results.num_biopac_cross;
    else
        results.deltaZC = NaN; 
    end
    
    results.display_str = sprintf('SCR Analysis | BIOPAC ZC: %d | DC-EDA ZC: %d | ΔZC: %.3f', ...
            results.num_biopac_cross, results.num_dc_cross, results.deltaZC);
    disp(results.display_str);
    
    results.dc_cross_idx     = idx_dc;
    results.biopac_cross_idx = idx_biopac;
    results.dc_cross_times     = dc_eda.time(idx_dc);
    results.biopac_cross_times = biopac_eda.time(idx_biopac);
end

% =========================================================================
% analyzeEvents.m
% =========================================================================
function results = analyzeEvents(dc_eda, biopac_eda, events)
    numEvents = size(events, 1);
    
    results.means_bio = zeros(numEvents, 1);
    results.stds_bio = zeros(numEvents, 1);
    results.means_dc = zeros(numEvents, 1);
    results.stds_dc = zeros(numEvents, 1);
    
    for i = 1:numEvents
        eventTime = events{i, 2};
        idx = biopac_eda.time >= eventTime(1) & biopac_eda.time <= eventTime(2);
        
        results.means_bio(i) = mean(biopac_eda.normalized(idx));
        results.stds_bio(i) = std(biopac_eda.normalized(idx));
        results.means_dc(i) = mean(dc_eda.normalized(idx));
        results.stds_dc(i) = std(dc_eda.normalized(idx));
    end
    
    baseline_bio = results.means_bio(1);
    baseline_dc = results.means_dc(1);
    results.perc_change_bio = ((results.means_bio(2:end) - baseline_bio) ./ baseline_bio) * 100;
    results.perc_change_dc = ((results.means_dc(2:end) - baseline_dc) ./ baseline_dc) * 100;
end

% =========================================================================
% plotSignalOverview.m
% =========================================================================
function plotSignalOverview(dc_eda, biopac_eda, scl_results, scr_results)
    figure('Name', 'EDA Signal Analysis Overview', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    
    subplot(2, 2, 1);
    hold on;
    plot(dc_eda.time, dc_eda.raw, 'r-', 'DisplayName', 'DC-EDA Raw');
    plot(dc_eda.time, dc_eda.cleaned, 'r--', 'LineWidth', 1.5, 'DisplayName', 'DC-EDA Cleaned');
    plot(biopac_eda.time, biopac_eda.raw, 'b-', 'DisplayName', 'BIOPAC Raw');
    plot(biopac_eda.time, biopac_eda.cleaned, 'b--', 'LineWidth', 1.5, 'DisplayName', 'BIOPAC Cleaned');
    hold off;
    grid on;
    title('Raw vs. Cleaned Signals (Artifacts Corrected)');
    xlabel('Time (minutes)');
    ylabel('EDA (\muS)');
    legend('show', 'Location', 'best');
    
    subplot(2, 2, 2);
    hold on;
    yyaxis left;
    plot(biopac_eda.time, biopac_eda.filtered_10Hz, 'b', 'DisplayName', 'BIOPAC');
    ylabel('BIOPAC EDA (\muS)');
    yyaxis right;
    plot(dc_eda.time, dc_eda.filtered_10Hz, 'r', 'DisplayName', 'DC-EDA');
    ylabel('DC-EDA EDA (\muS)');
    hold off;
    grid on;
    title('10 Hz Low-pass Filtered Signals');
    xlabel('Time (minutes)');
    legend('show', 'Location', 'best');
    
    subplot(2, 2, 3);
    hold on;
    yyaxis left;
    plot(biopac_eda.time, biopac_eda.scl, 'b', 'DisplayName', 'BIOPAC SCL');
    ylabel('BIOPAC SCL (\muS)');
    yyaxis right;
    plot(dc_eda.time, dc_eda.scl, 'r', 'DisplayName', 'DC-EDA SCL');
    ylabel('DC-EDA SCL (\muS)');
    hold off;
    grid on;
    title('Tonic Component (SCL) Comparison');
    xlabel('Time (minutes)');
    legend('show', 'Location', 'best');
    annotation('textbox', [0.1, 0.05, 0.4, 0.05], 'String', scl_results.display_str, ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 9);

    subplot(2, 2, 4);
    hold on;
    yyaxis left;
    plot(biopac_eda.time, biopac_eda.scr, 'b', 'DisplayName', 'BIOPAC SCR');
    plot(scr_results.biopac_cross_times, biopac_eda.scr(scr_results.biopac_cross_idx), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', 'BIOPAC ZCR');
    ylabel('BIOPAC SCR (\muS)');
    
    yyaxis right;
    plot(dc_eda.time, dc_eda.scr, 'r', 'DisplayName', 'DC-EDA SCR');
    plot(scr_results.dc_cross_times, dc_eda.scr(scr_results.dc_cross_idx), 'ks', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'DC-EDA ZCR');
    ylabel('DC-EDA SCR (\muS)');
    
    hold off;
    grid on;
    title('Phasic Component (SCR) Comparison with Zero-Crossings');
    xlabel('Time (minutes)');
    legend('show', 'Location', 'best');
    annotation('textbox', [0.55, 0.05, 0.4, 0.05], 'String', scr_results.display_str, ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 9);
end

% =========================================================================
% plotEventAnalysis.m
% =========================================================================
function plotEventAnalysis(results, events)
    figure('Name', 'Event-Based Analysis', 'NumberTitle', 'off', 'Position', [200, 200, 1000, 500]);
    event_labels = events(:,1);
    
    subplot(1, 2, 1);
    bar_data = [results.means_bio, results.means_dc];
    b = bar(bar_data, 'grouped');
    hold on;
    
    [ngroups, nbars] = size(bar_data);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        err = (i==1) * results.stds_bio + (i==2) * results.stds_dc;
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, bar_data(:,i), err, 'k', 'linestyle', 'none', 'LineWidth', 1);
    end
    
    hold off;
    set(gca, 'XTickLabel', event_labels);
    ylabel('Mean Normalized EDA (0-1)');
    xlabel('Event');
    legend('BIOPAC', 'DC-EDA', 'Location', 'northwest');
    title({'Mean \pm SD of Normalized EDA', 'across Events'});
    grid on;
    
    subplot(1, 2, 2);
    perc_change_data = [results.perc_change_bio, results.perc_change_dc];
    bar(perc_change_data, 'grouped');
    
    set(gca, 'XTickLabel', event_labels(2:end));
    ylabel('Percentage Change from Baseline (%)');
    xlabel('Event');
    legend('BIOPAC', 'DC-EDA', 'Location', 'northwest');
    title({'Percentage Change in Mean EDA', sprintf('from Baseline (%s)', event_labels{1})});
    grid on;
end

% =========================================================================
% plotDecomposition.m
% =========================================================================
function plotDecomposition(dc_eda, biopac_eda)
    figure('Name', 'EDA Signal Decomposition', 'NumberTitle', 'off', 'Position', [300, 300, 1000, 700]);
    
    subplot(2, 1, 1);
    hold on;
    plot(biopac_eda.time, biopac_eda.filtered_10Hz, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Overall EDA (Filtered)');
    plot(biopac_eda.time, biopac_eda.scl, 'b--', 'LineWidth', 1.5, 'DisplayName', 'SCL (Tonic)');
    plot(biopac_eda.time, biopac_eda.scr, 'g:', 'LineWidth', 1.5, 'DisplayName', 'SCR (Phasic)');
    hold off;
    grid on;
    title('BIOPAC Signal Decomposition');
    xlabel('Time (minutes)');
    ylabel('EDA (\muS)');
    legend('show', 'Location', 'best');
    
    subplot(2, 1, 2);
    hold on;
    plot(dc_eda.time, dc_eda.filtered_10Hz, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Overall EDA (Filtered)');
    plot(dc_eda.time, dc_eda.scl, 'r--', 'LineWidth', 1.5, 'DisplayName', 'SCL (Tonic)');
    plot(dc_eda.time, dc_eda.scr, 'm:', 'LineWidth', 1.5, 'DisplayName', 'SCR (Phasic)');
    hold off;
    grid on;
    title('DC-EDA Signal Decomposition');
    xlabel('Time (minutes)');
    ylabel('EDA (\muS)');
    legend('show', 'Location', 'best');
end
