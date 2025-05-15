%This script processes multiple simulation files and calculated the mean alignment



clear;

file_patterns_smooth = { ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_smooth_2f_run%d.mat',  [1 2 3 4], [28 28 38 19]; ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_smooth_1f_run%d.mat',  [1 2 3 4], [20 23 28 30] ...
    };
file_patterns_canyon = { ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_canyon_1f_run%d.mat',  [1 2 3 4 5 6], [19 32 41 30 26 27]; ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_canyon_2f_run%d.mat',  [3 4 5 6 7 8], [32 14 18 36 37 43]; ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_canyon_2f_run%d.mat',  [1 2],         [19 16]; ...
    %'/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_canyon_1f_run%d.mat',  [7 8],         [22 36] ...
    };
file_patterns_ridge = { ...
     '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_ridge_2f_run%d.mat',  [1 2 3 4], [60 29 33 35]; ...
     '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_ridge_1f_run%d.mat',  [1 2 3 4], [24 31 26 28] ...
    };

file_patterns = [file_patterns_smooth; file_patterns_canyon; file_patterns_ridge];

all_mat_files_to_load = {}; all_runs = []; all_topographies = {};
all_speeds = []; all_startframes = []; all_directions = {};


for p = 1:size(file_patterns,1)
    pattern_ref = file_patterns{p,1}; runs = file_patterns{p,2}; starts = file_patterns{p,3};
    if length(runs) ~= length(starts); error('Mismatch runs/starts: %s', pattern_ref); end
    parts = split(pattern_ref, '/'); filename_part = parts{end}; name_parts = split(filename_part, '_');
    feat_str = name_parts{2}; speed_str = name_parts{3}; speed_val = str2double(speed_str(1));
    for r = 1:length(runs)
        current_run = runs(r); current_start = starts(r);
        matFilename = sprintf(pattern_ref, current_run);
        all_mat_files_to_load{end+1} = matFilename; all_runs(end+1) = current_run;
        all_topographies{end+1} = feat_str; all_speeds(end+1) = speed_val;
        all_startframes(end+1) = current_start;
        
        if mod(current_run, 2) == 1
            all_directions{end+1} = 'prograde';
        else
            all_directions{end+1} = 'retrograde';
        end
    end
end


numFiles = length(all_mat_files_to_load);
numFramesToPlot = 500;
smoothing_window = 20;

results_mean = cell(numFiles, 1);


f = 2 * 2 * pi * (8 / 60);


for i = 1:numFiles
    matFilename = all_mat_files_to_load{i};
    run_num = all_runs(i); topo = all_topographies{i}; speed = all_speeds(i); direction = all_directions{i};
    fprintf('Processing file %d/%d: %s (Run %d, Topo: %s, Speed: %df, Dir: %s)\n', ...
            i, numFiles, matFilename, run_num, topo, speed, direction);
    try
        loaded_data = load(matFilename, 'x_masked_all', 'y_masked_all', ...
                           'masked_normInterpolated_dzy_all', 'masked_normU_all', ...
                           'masked_normInterpolated_dzx_all', 'masked_normV_all');
        x_masked_all = loaded_data.x_masked_all; y_masked_all = loaded_data.y_masked_all;
        dzy_all = loaded_data.masked_normInterpolated_dzy_all; normU_all = loaded_data.masked_normU_all;
        dzx_all = loaded_data.masked_normInterpolated_dzx_all; normV_all = loaded_data.masked_normV_all;
    catch ME
        warning('Could not load variables in %s. Skipping. Error: %s', matFilename, ME.message);
        results_mean{i} = NaN;
        continue;
    end
    num_actual_timesteps = length(dzy_all);
    if num_actual_timesteps == 0
        warning('Empty data in %s. Skipping.', matFilename);
        results_mean{i} = NaN;
        continue;
    end


    mean_metric_run = nan(1, num_actual_timesteps);
    for t_idx = 1:num_actual_timesteps
        try
            x_current_masked = x_masked_all{t_idx}; y_current_masked = y_masked_all{t_idx};
            dzy_current = dzy_all{t_idx}; normU_current = normU_all{t_idx};
            dzx_current = dzx_all{t_idx}; normV_current = normV_all{t_idx};
        catch ME_extract
             warning('Error extracting t=%d in %s. Skip step. Err: %s', t_idx, matFilename, ME_extract.message); continue;
        end
        if isempty(x_current_masked) || isempty(dzy_current) || isempty(normU_current) || isempty(dzx_current) || isempty(normV_current); continue; end
        try
            cosTheta_current_masked = -dzy_current .* normU_current + dzx_current .* normV_current;
            cosTheta_current_masked = max(min(cosTheta_current_masked, 1), -1);
        catch ME_cosTheta
            warning('Error calculating cosTheta t=%d in %s. Skip step. Err: %s', t_idx, matFilename, ME_cosTheta.message); continue;
        end

        try
            if strcmp(topo, 'canyon') && speed == 2
                if contains(matFilename, '28_canyon_2f'); spatial_mask = (x_current_masked >= -8 & x_current_masked <= 8) & (y_current_masked >= 6 & y_current_masked <= 14);
                elseif contains(matFilename, '29_canyon_2f'); spatial_mask = (x_current_masked >= -8 & x_current_masked <= 8) & (y_current_masked >= 6 & y_current_masked <= 14.5);
                else; spatial_mask = false(size(x_current_masked)); warning('Filename %s matched canyon/2f but not specific 28/29 pattern.', matFilename); end
            elseif strcmp(topo, 'canyon') && speed == 1
                 if contains(matFilename, '28_canyon_1f'); spatial_mask = (x_current_masked >= -8 & x_current_masked <= 8) & (y_current_masked >= 6 & y_current_masked <= 14);
                 elseif contains(matFilename, '29_canyon_1f'); spatial_mask = (x_current_masked >= -8 & x_current_masked <= 8) & (y_current_masked >= 6 & y_current_masked <= 14.5);
                 else; spatial_mask = false(size(x_current_masked)); warning('Filename %s matched canyon/1f but not specific 28/29 pattern.', matFilename); end
            elseif strcmp(topo, 'smooth'); spatial_mask = (x_current_masked >= -7 & x_current_masked <= 8) & (y_current_masked >= 7.5 & y_current_masked <= 14.5);
            elseif strcmp(topo, 'ridge'); spatial_mask = (x_current_masked >= -7 & x_current_masked <= 8) & (y_current_masked >= 6 & y_current_masked <= 15);
            else; error('Unknown topography/speed combo: %s', matFilename); end
            % --- End ROI Mask ---
            cosTheta_in_roi = cosTheta_current_masked(spatial_mask);
            alignment_metric = abs(cosTheta_in_roi);
        catch ME_mask_metric
            warning('Error masking/metric t=%d in %s. Skip step. Err: %s', t_idx, matFilename, ME_mask_metric.message); continue;
        end
        if (~isempty(alignment_metric))
            mean_metric_run(t_idx) = mean(alignment_metric, 'omitnan');
        end
    end
    results_mean{i} = mean_metric_run;
end

% non-dimensional time
max_relative_time = (numFramesToPlot / 5 - 1/5) * f;
common_time_grid = linspace(0, max_relative_time, numFramesToPlot);

% for std plotting
shading_alpha = 0.2;

% creating subgroups for plotting
% Speed, Direction, Label, Color, LineStyle
subgroups = { ...
    1, 'prograde', '-1rpm, Prograde', 'b', '-'; ... 
    1, 'retrograde', '+1rpm, Retrograde', 'r', '-'; ...
    2, 'prograde', '-2rpm, Prograde', 'b', '--'; ...
    2, 'retrograde', '+2rpm, Retrograde', 'r', '--' ...
};
num_subgroups = size(subgroups, 1);

unique_topos = {'smooth', 'canyon', 'ridge'}; % Process in this order

for t = 1:length(unique_topos)
    current_topo = unique_topos{t};
    fprintf('\nGenerating plot for topography: %s\n', current_topo);

    figure('Renderer', 'painters', 'Position', [100, 100, 1000, 600]); hold on; % Figure for individual runs % Create a new figure for this topography
    hold on;
    legend_entries = {};
    plot_handles = [];

    for sg = 1:num_subgroups
        target_speed = subgroups{sg, 1};
        target_direction = subgroups{sg, 2};
        subgroup_label = subgroups{sg, 3};
        subgroup_color = subgroups{sg, 4};
        subgroup_linestyle = subgroups{sg, 5};

        fprintf('  Processing subgroup: %s\n', subgroup_label);

        % Find indices for the subgroups (topo, speed, direction)
        indices_subgroup = find(strcmp(all_topographies, current_topo) & ...
                                (all_speeds == target_speed) & ...
                                strcmp(all_directions, target_direction));

        if isempty(indices_subgroup)
            fprintf('    WARNING: No runs found for subgroup: %s\n', subgroup_label);
            continue; % Skip to next subgroup if no runs match
        end

        interpolated_means_group = nan(numFramesToPlot, length(indices_subgroup)); % Store interpolated MEANS
        valid_runs_in_group = 0;

        for k = 1:length(indices_subgroup)
            idx = indices_subgroup(k);
            if isscalar(results_mean{idx}) && isnan(results_mean{idx}); continue; end % Use results_mean

            
            full_mean = results_mean{idx};
            startframe = all_startframes(idx);
            num_available_frames = length(full_mean);
            start_idx = startframe;
            end_idx = min(startframe + numFramesToPlot - 1, num_available_frames);
            if start_idx < 1 || start_idx > num_available_frames || start_idx > end_idx; continue; end
            mean_to_plot = full_mean(start_idx:end_idx);
            frames_plot = start_idx:end_idx;
            xTime_abs = (frames_plot / 5 - 1/5) * f;
            t_start_f = (start_idx / 5 - 1/5) * f;
            xTime_relative = xTime_abs - t_start_f;
            valid_indices_plot = find(~isnan(mean_to_plot));
            if isempty(valid_indices_plot); continue; end
            time_valid   = xTime_relative(valid_indices_plot);
            mean_valid = mean_to_plot(valid_indices_plot); 
            
            [time_valid_unique, ia, ~] = unique(time_valid);
            mean_valid_unique = mean_valid(ia); % Use mean data
            if length(time_valid_unique) < 2; continue; end
            interpolated_means_group(:, k) = interp1(time_valid_unique, mean_valid_unique, common_time_grid, 'linear', NaN);
            valid_runs_in_group = valid_runs_in_group + 1;
        end

        % Calculate and Plot Average of mean +/- Std Dev of alignment
        if valid_runs_in_group > 1
            avg_of_means = mean(interpolated_means_group, 2, 'omitnan');
            std_of_means = std(interpolated_means_group, 0, 2, 'omitnan');
            upper_bound = avg_of_means + std_of_means;
            lower_bound = avg_of_means - std_of_means;
            avg_of_means_smoothed = movmean(avg_of_means, smoothing_window, 'omitnan');
            valid_plot_indices = find(~isnan(avg_of_means_smoothed));
            if isempty(valid_plot_indices);
                 fprintf('    WARNING: Average calculation resulted in NaNs for %s\n', subgroup_label);
                 continue;
            end
            time_grid_valid = common_time_grid(valid_plot_indices);
            avg_valid = avg_of_means_smoothed(valid_plot_indices);
            upper_valid = upper_bound(valid_plot_indices);
            lower_valid = lower_bound(valid_plot_indices);
            upper_valid = min(upper_valid, 1); lower_valid = max(lower_valid, 0);

            % std dev shading
            fill([time_grid_valid, fliplr(time_grid_valid)], ...
                 [lower_valid', fliplr(upper_valid')], ...
                 subgroup_color, 'FaceAlpha', shading_alpha, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
            h = plot(time_grid_valid, avg_valid, 'Color', subgroup_color, ...
                     'LineWidth', 2, 'LineStyle', subgroup_linestyle);
            plot_handles(end+1) = h;
            legend_entries{end+1} = sprintf('%s', subgroup_label);
            fprintf('    Plotted Avg of Mean +/- Std Dev for %s\n', subgroup_label);

        elseif valid_runs_in_group == 1
             avg_of_means = interpolated_means_group(:, ~all(isnan(interpolated_means_group)));
             avg_of_means_smoothed = movmean(avg_of_means, smoothing_window, 'omitnan');
             valid_plot_indices = find(~isnan(avg_of_means_smoothed));
             if isempty(valid_plot_indices);
                 fprintf('    WARNING: Single run calculation resulted in NaNs for %s\n', subgroup_label);
                 continue;
             end
             time_grid_valid = common_time_grid(valid_plot_indices);
             avg_valid = avg_of_means_smoothed(valid_plot_indices);

             h = plot(time_grid_valid, avg_valid, 'Color', subgroup_color, ...
                      'LineWidth', 2, 'LineStyle', subgroup_linestyle);
             plot_handles(end+1) = h;
             legend_entries{end+1} = sprintf('%s (1 run)', subgroup_label);
             fprintf('    Plotted single run mean (no Std Dev) for %s\n', subgroup_label);
        else
            fprintf('    Skipping plot for %s (no valid runs)\n', subgroup_label);
        end

    end

    % final plotting
    if ~isempty(plot_handles)
        xlabel(['Elapsed $t \cdot f$ since start'], 'interpreter', 'latex');
        ylabel(['Mean Alignment $|\cos{\theta}|$'], 'interpreter', 'latex');
        %title(sprintf('Average of Mean ± Std Dev Alignment Metric - Topography: %s', current_topo)); % Updated Title
        legend(plot_handles, legend_entries, 'Location', 'best');
        grid on;
        ylim([0 1]); yticks(0:0.2:1);
        hold off;

    else
        close(gcf);
        fprintf('WARNING: No data plotted for topography %s.\n', current_topo);
    end

end