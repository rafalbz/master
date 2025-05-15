%
%This script processes multiple files and calculates the total kinetic energy (TKE) and RMS velocity for different experimental conditions




file_patterns_smooth = { ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_smooth_2f_run%d.mat',  [1 2 3 4], [28 28 38 19]; ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_smooth_1f_run%d.mat',  [1 2 3 4], [20 23 28 30] ...
    };
file_patterns_canyon = { ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_canyon_1f_run%d.mat',  [1 2 3 4 5 6], [19 32 41 30 26 27]; ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_canyon_2f_run%d.mat',  [1 2], [19 16]; ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_canyon_2f_run%d.mat',  [3 4 5 6 7 8], [32 14 18 36 37 43]; ...
    '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/29_canyon_1f_run%d.mat',  [7 8], [22 36] ...
    };
file_patterns_ridge = { ...
     '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_ridge_2f_run%d.mat',  [1 2 3 4], [60 29 33 35]; ...
     '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/results/28_ridge_1f_run%d.mat',  [1 2 3 4], [24 31 26 28] ...
    };


file_patterns = [file_patterns_smooth; file_patterns_canyon ;file_patterns_ridge];


all_files = {};
all_runs = [];
all_startframes = [];
for p = 1:size(file_patterns,1)
    pattern = file_patterns{p,1};
    runs = file_patterns{p,2};
    starts = file_patterns{p,3};
    if length(runs) ~= length(starts)
        error('Mismatch between number of runs and start frames for pattern: %s', pattern);
    end
    for r = 1:length(runs)
        all_files{end+1} = sprintf(pattern, runs(r));
        all_runs(end+1) = runs(r);
        all_startframes(end+1) = starts(r);
    end
end


numFiles = length(all_files);
numFramesToAnalyze = 500;
smoothing_window = 30;

figure('Renderer', 'painters', 'Position', [100, 100, 1000, 600]); hold on;
colors = lines(numFiles);
all_normalized_energies = nan(numFiles, numFramesToAnalyze);
all_delta_A = nan(numFiles, 1);

f = 2 * 2 * pi * (8 / 60);

fprintf('Processing %d files...\n', numFiles);

for i = 1:numFiles
    matFilename = all_files{i};
    run = all_runs(i);
    startframe = all_startframes(i);
    fprintf('Processing file %d/%d: %s (Run %d, Start %d)\n', i, numFiles, matFilename, run, startframe);

    try
        load(matFilename, 'U', 'V', 'x', 'y');
    catch ME
        warning('Could not load file %s. Skipping. Error: %s', matFilename, ME.message);
        all_delta_A(i) = NaN; % Ensure it's NaN if file load fails
        continue; % Skip to next file
    end


    if isempty(x) || isempty(y) || isempty(x{1}) || isempty(y{1})
        warning('Coordinate data x{1} or y{1} is missing/empty for file %s. Skipping file.', matFilename);
        all_normalized_energies(i, :) = NaN;
        all_delta_A(i) = NaN;
        continue;
    end
    
    x_grid = x{1};
    y_grid = y{1};

    % dx and dy
    current_dx = 1;
    current_dy = 1;
    if size(x_grid, 2) > 1
        current_dx = abs(x_grid(1, 2) - x_grid(1, 1));
    else
        warning('x_grid has only one column for file %s. dx assumed to be 1. Area calculation might be incorrect.', matFilename);
    end
    if size(y_grid, 1) > 1
        current_dy = abs(y_grid(2, 1) - y_grid(1, 1));
    else
        warning('y_grid has only one row for file %s. dy assumed to be 1. Area calculation might be incorrect.', matFilename);
    end
    current_delta_A_val = current_dx * current_dy;
    if current_delta_A_val == 0
        warning('Calculated delta_A is zero for file %s. Check grid data. Defaulting to 1 to avoid issues, but results will be incorrect.', matFilename);
        current_delta_A_val = 1; 
    end
    all_delta_A(i) = current_delta_A_val;

    total_energy = nan(1, numFramesToAnalyze);

    for idx = 1:numFramesToAnalyze
        frameIndex = startframe + idx - 1;

        % --- Robustness Checks ---
        if frameIndex < 1 || frameIndex > length(U) || isempty(U{frameIndex})
             warning('Invalid or empty U data for frame %d in file %s. Skipping frame.', frameIndex, matFilename);
             continue; % Skip this frame index
        end
         if frameIndex < 1 || frameIndex > length(V) || isempty(V{frameIndex})
             warning('Invalid or empty V data for frame %d in file %s. Skipping frame.', frameIndex, matFilename);
             continue; % Skip this frame index
        end
         if isempty(x) || isempty(y) || isempty(x{1}) || isempty(y{1})
            warning('Coordinate data x{1} or y{1} is missing/empty for file %s. Skipping file.', matFilename);
            total_energy(:) = NaN;
            break;
         end

        x_grid = x{1};
        y_grid = y{1};
        U_frame = U{frameIndex};
        V_frame = V{frameIndex};

        % Mask for each topographic region
        try
            if contains(matFilename, '28_canyon_2f')
                mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14);
            elseif contains(matFilename, '28_canyon_1f')
                mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14);
            elseif contains(matFilename, '29_canyon_2f')
                mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14.5);
            elseif contains(matFilename, '29_canyon_1f')
                mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14.5);
            elseif contains(matFilename, '29_smooth_2f')
                mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15.5);
            elseif contains(matFilename, '29_smooth_1f')
                mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15.5);
            elseif contains(matFilename, '28_ridge_2f')
                mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15);
            elseif contains(matFilename, '28_ridge_1f')
                mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15);
            else
                error('Unknown file or run configuration for file: %s', matFilename);
            end
        catch ME_mask
             warning('Error creating mask for frame %d in file %s. Skipping frame. Error: %s', frameIndex, matFilename, ME_mask.message);
             continue;
        end

        if ~isequal(size(mask), size(U_frame))
             warning('Mask size [%s] does not match U_frame size [%s] for frame %d in file %s. Skipping frame.', ...
                     num2str(size(mask)), num2str(size(U_frame)), frameIndex, matFilename);
             continue;
        end

        U_masked = U_frame .* mask;
        V_masked = V_frame .* mask;
        xenergy = 0.5 * (U_masked.^2 + V_masked.^2);
        total_energy(idx) = sum(xenergy(:), 'omitnan') * current_delta_A_val;
    end 


    max_energy = max(total_energy);
    if ~isempty(max_energy) && max_energy > 0 && ~isnan(max_energy)
        normalized_energy = total_energy ;
    else
        normalized_energy = nan(1, numFramesToAnalyze); % Assign NaN if max_energy is invalid or zero
        warning('Max energy is zero, NaN, or empty for file %s. Normalization resulted in NaNs.', matFilename);
    end
    all_normalized_energies(i, :) = normalized_energy;

    % Smooth the normalized energy (handle NaNs)
    total_energy_smooth = movmean(normalized_energy, 1, 'omitnan');

    % Plot individual run (only if data is valid)
    if any(~isnan(total_energy_smooth))
        xData = (1:numFramesToAnalyze);
        xPlot = (xData / 5 - 1 / 5) * f; % Non-dimensional time
        % xPlot = xData;
%{
        if mod(run,2)==1
            plot(xPlot, total_energy_smooth, ':', 'Color', colors(i,:) ,'LineWidth', 1, ...
                'DisplayName', sprintf('Prograde (Run %d, Start %d)', run, startframe));
        else
            plot(xPlot, total_energy_smooth, '--', 'Color', colors(i,:) ,'LineWidth', 1, ...
                'DisplayName', sprintf('Retrograde (Run %d, Start %d)', run, startframe));
        end
%}
    end
end


is_prograde = mod(all_runs, 2) == 1;
is_retrograde = mod(all_runs, 2) == 0;
idx_prograde_all = find(is_prograde);
idx_retrograde_all = find(is_retrograde);


% Use the safe calculation function (assuming it's defined later or copy it here if needed)
calculate_mean_smooth = @(indices, data, window) ...
     ( ~isempty(indices) && sum(~isnan(data(indices,:)), 'all') > 0 ) ...
     * movmean(mean(data(indices,:), 1, 'omitnan'), window, 'omitnan');

mean_prograde_all_smooth = calculate_mean_smooth(idx_prograde_all, all_normalized_energies, smoothing_window);
mean_retrograde_all_smooth = calculate_mean_smooth(idx_retrograde_all, all_normalized_energies, smoothing_window);


xData = (1:numFramesToAnalyze);
xPlot = (xData / 5 - 1 / 5) * f; % Non-dimensional time

%if any(mean_prograde_all_smooth)
%    plot(xPlot, mean_prograde_all_smooth, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Overall Mean Prograde');
%end
%if any(mean_retrograde_all_smooth)
%    plot(xPlot, mean_retrograde_all_smooth, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Overall Mean Retrograde');
%end

% Find indices for each category
is_1f = contains(all_files, '_1f_');
is_2f = contains(all_files, '_2f_');
is_canyon = contains(all_files, 'canyon');
is_ridge = contains(all_files, 'ridge');
is_smooth = contains(all_files, 'smooth');
% is_prograde and is_retrograde are already defined

idx_1f_pro_can = find(is_1f & is_prograde & is_canyon);
idx_1f_pro_smo = find(is_1f & is_prograde & is_smooth);
idx_1f_pro_rid = find(is_1f & is_prograde & is_ridge);
idx_1f_retro_can = find(is_1f & is_retrograde & is_canyon);
idx_1f_retro_smo = find(is_1f & is_retrograde & is_smooth);
idx_1f_retro_rid = find(is_1f & is_retrograde & is_ridge);
idx_2f_pro_can = find(is_2f & is_prograde & is_canyon);
idx_2f_pro_smo = find(is_2f & is_prograde & is_smooth);
idx_2f_pro_rid = find(is_2f & is_prograde & is_ridge);
idx_2f_retro_can = find(is_2f & is_retrograde & is_canyon);
idx_2f_retro_smo = find(is_2f & is_retrograde & is_smooth);
idx_2f_retro_rid = find(is_2f & is_retrograde & is_ridge);

idx_1f_retro = find(is_1f & is_retrograde);
idx_1f_pro = find(is_1f & is_prograde);
idx_2f_pro = find(is_2f & is_prograde);
idx_2f_retro = find(is_2f & is_retrograde);

% Calculate smoothed means (using the same safe function)
mean_1f_pro_can_smooth = calculate_mean_smooth(idx_1f_pro_can, all_normalized_energies, smoothing_window);
mean_1f_pro_smo_smooth = calculate_mean_smooth(idx_1f_pro_smo, all_normalized_energies, smoothing_window);
mean_1f_pro_rid_smooth = calculate_mean_smooth(idx_1f_pro_rid, all_normalized_energies, smoothing_window);
mean_1f_retro_can_smooth = calculate_mean_smooth(idx_1f_retro_can, all_normalized_energies, smoothing_window);
mean_1f_retro_smo_smooth = calculate_mean_smooth(idx_1f_retro_smo, all_normalized_energies, smoothing_window);
mean_1f_retro_rid_smooth = calculate_mean_smooth(idx_1f_retro_rid, all_normalized_energies, smoothing_window);
mean_2f_pro_can_smooth = calculate_mean_smooth(idx_2f_pro_can, all_normalized_energies, smoothing_window);
mean_2f_pro_smo_smooth = calculate_mean_smooth(idx_2f_pro_smo, all_normalized_energies, smoothing_window);
mean_2f_pro_rid_smooth = calculate_mean_smooth(idx_2f_pro_rid, all_normalized_energies, smoothing_window);
mean_2f_retro_can_smooth = calculate_mean_smooth(idx_2f_retro_can, all_normalized_energies, smoothing_window);
mean_2f_retro_smo_smooth = calculate_mean_smooth(idx_2f_retro_smo, all_normalized_energies, smoothing_window);
mean_2f_retro_rid_smooth = calculate_mean_smooth(idx_2f_retro_rid, all_normalized_energies, smoothing_window);

tke_conversion_factor = 1;

% Helper function to plot individual runs for a given category
plot_individual_runs = @(indices, data, x_axis, color, style, smooth_win, conv_factor) ...
    arrayfun(@(idx) plot(x_axis, movmean(data(idx,:), smooth_win, 'omitnan') * conv_factor, ...
                         'Color', color, 'LineStyle', style, ...
                         'LineWidth', 0.1, 'HandleVisibility', 'off'), ...
             indices(any(~isnan(movmean(data(indices,:), smooth_win, 'omitnan')), 2)));

% --- Plot individual runs for 1f smooth ---
% Plot Prograde runs in blue
%plot_individual_runs(idx_1f_pro_smo, all_normalized_energies, xPlot, 'b', '-', 1,tke_conversion_factor);
%plot_individual_runs(idx_2f_pro_smo, all_normalized_energies, xPlot, 'b', '-', 1,tke_conversion_factor);
% Plot Retrograde runs in red
%plot_individual_runs(idx_1f_retro_smo, all_normalized_energies, xPlot, 'r', '-', 1,tke_conversion_factor);
%plot_individual_runs(idx_2f_retro_smo, all_normalized_energies, xPlot, 'r', '-', 1,tke_conversion_factor);
%plot_individual_runs(idx_1f_pro_can, all_normalized_energies, xPlot, 'b', '-', 1, tke_conversion_factor);
%plot_individual_runs(idx_1f_retro_can, all_normalized_energies, xPlot, 'r', '-', 1, tke_conversion_factor);
plot_individual_runs(idx_2f_pro_can, all_normalized_energies, xPlot, 'b', '-', 1, tke_conversion_factor);
plot_individual_runs(idx_2f_retro_can, all_normalized_energies, xPlot, 'r', '-', 1, tke_conversion_factor);
%plot_individual_runs(idx_1f_pro_rid, all_normalized_energies, xPlot, 'b', '-', 1, tke_conversion_factor);
%plot_individual_runs(idx_1f_retro_rid, all_normalized_energies, xPlot, 'r', '-', 1, tke_conversion_factor);
%plot_individual_runs(idx_2f_pro_rid, all_normalized_energies, xPlot, 'b', '-', 1, tke_conversion_factor);
%plot_individual_runs(idx_2f_retro_rid, all_normalized_energies, xPlot, 'r', '-', 1, tke_conversion_factor);

% Plot smoothed means if they contain valid data
%plot(xPlot, mean_1f_pro_can_smooth * tke_conversion_factor, 'b-', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde canyon');
%plot(xPlot, mean_1f_pro_smo_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde smooth');
%plot(xPlot, mean_1f_pro_rid_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde ridge');
%plot(xPlot, mean_1f_retro_can_smooth * tke_conversion_factor, 'r-', 'LineWidth', 4, 'DisplayName', ' -1rpm, Retrograde canyon');
%plot(xPlot, mean_1f_retro_smo_smooth, 'r-', 'LineWidth', 4, 'DisplayName', '-1rpm, Retrograde smooth');
%plot(xPlot, mean_1f_retro_rid_smooth, 'r-', 'LineWidth', 4, 'DisplayName', '-1rpm, Retrograde ridge');
plot(xPlot, mean_2f_pro_can_smooth * tke_conversion_factor, 'b-', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde canyon');
%plot(xPlot, mean_2f_pro_smo_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde smooth');
%plot(xPlot, mean_2f_pro_rid_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde ridge');
plot(xPlot, mean_2f_retro_can_smooth * tke_conversion_factor, 'r-', 'LineWidth', 4, 'DisplayName', ' -2rpm, Retrograde canyon');
%plot(xPlot, mean_2f_retro_smo_smooth, 'r-', 'LineWidth', 4, 'DisplayName', '-2rpm, Retrograde smooth');
%plot(xPlot, mean_2f_retro_rid_smooth, 'r-', 'LineWidth', 4, 'DisplayName', ' -2rpm, Retrograde ridge');

%mean_1f_retro_smooth = calculate_mean_smooth(idx_1f_retro, all_normalized_energies, smoothing_window);
%mean_2f_pro_smooth = calculate_mean_smooth(idx_2f_pro, all_normalized_energies, smoothing_window);
%mean_2f_retro_smooth = calculate_mean_smooth(idx_2f_retro, all_normalized_energies, smoothing_window);
%plot(xPlot, mean_1f_pro_smooth, 'b:', 'LineWidth', 3, 'DisplayName', 'Mean 1f Prograde');
%plot(xPlot, mean_1f_retro_smooth, 'r:', 'LineWidth', 3, 'DisplayName', 'Mean 1f Retrograde');
%plot(xPlot, mean_2f_pro_smooth, 'b-', 'LineWidth', 3, 'DisplayName', 'Mean 2f Prograde');
%plot(xPlot, mean_2f_retro_smooth, 'r-', 'LineWidth', 3, 'DisplayName', 'Mean 2f Retrograde');


hold off; % Now turn hold off for the first figure
xlabel(['Elapsed $t \cdot f$ since start'], 'interpreter', 'latex');
ylabel('Integral of Total Kinetic Energy [cm^4s^{-2}]');
%title('Integral of Total Kinetic Energy');
legend('show', 'Location', 'best');
grid on;
set(gca, 'YMinorGrid', 'on');
set(gca, 'YScale', 'log');


% 
% RMS Velocity
% 

fprintf('\nCalculating RMS velocity...\n');

numFiles = length(all_files);
numFramesToAnalyze = size(all_normalized_energies, 2);
num_points_in_mask = zeros(numFiles, 1);

% all_average_energies will store mean_spatial(0.5*(U^2+V^2))*delta_A
all_average_energies_intermediate = nan(numFiles, numFramesToAnalyze); 
all_rms_velocity = nan(numFiles, numFramesToAnalyze);


% --- Loop to determine number of points in mask for each file ---
for i = 1:numFiles
    matFilename = all_files{i};

    if isnan(all_delta_A(i))
        warning('delta_A is NaN for file %s (index %d). Skipping RMS calculation.', matFilename, i);
        all_average_energies_intermediate(i, :) = nan(1, numFramesToAnalyze);
        all_rms_velocity(i, :) = nan(1, numFramesToAnalyze);
        continue;
    end

    try
        load(matFilename, 'x', 'y');
        if isempty(x) || isempty(y) || isempty(x{1}) || isempty(y{1})
             warning('Coordinate data x{1} or y{1} is missing/empty for file %s. Cannot calculate RMS velocity.', matFilename);
             num_points_in_mask(i) = 0;
             all_average_energies_intermediate(i, :) = nan(1, numFramesToAnalyze);
             all_rms_velocity(i, :) = nan(1, numFramesToAnalyze);
             continue;
        end
        x_grid = x{1};
        y_grid = y{1};

        if contains(matFilename, '28_canyon_2f')
            mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14);
        elseif contains(matFilename, '28_canyon_1f')
            mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14);
        elseif contains(matFilename, '29_canyon_2f')
            mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14.5);
        elseif contains(matFilename, '29_canyon_1f')
            mask = (x_grid >= -8 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 14.5);
        elseif contains(matFilename, '29_smooth_2f')
            mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15.5);
        elseif contains(matFilename, '29_smooth_1f')
            mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15.5);
        elseif contains(matFilename, '28_ridge_2f')
            mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15);
        elseif contains(matFilename, '28_ridge_1f')
            mask = (x_grid >= -7 & x_grid <= 8) & (y_grid >= 6 & y_grid <= 15);
        else
            warning('Unknown file type for mask calculation: %s', matFilename);
            mask = false(size(x_grid));
        end

        num_points_in_mask(i) = sum(mask(:));

    catch ME
        warning('Could not load file %s to determine mask size for RMS. Skipping. Error: %s', matFilename, ME.message);
        num_points_in_mask(i) = 0;
        all_average_energies_intermediate(i, :) = nan(1, numFramesToAnalyze);
        all_rms_velocity(i, :) = nan(1, numFramesToAnalyze);
        continue;
    end

    % RMS Velocity
    % stores sum(0.5 * (U^2+V^2)) * delta_A(i)
    % where delta_A(i) is in cm^2 and (U^2+V^2) is in cm^2/s^2
    % So all_normalized_energies is in cm^4/s^2
    if num_points_in_mask(i) > 0 && all_delta_A(i) > 0
        % all_average_energies_intermediate is (Integral of KE) / N_points
        % Units: cm^4/s^2
        all_average_energies_intermediate(i, :) = all_normalized_energies(i, :) / num_points_in_mask(i);

        % Calculate mean_spatial(U^2 + V^2)
        % mean_spatial(U^2+V^2) = ( (Integral of 0.5*(U^2+V^2)) / N_points * 2) / delta_A(i)
        % Units: ( (cm^4/s^2) * 2 ) / cm^2 = cm^2/s^2
        mean_sq_velocity_spatial_avg = (all_average_energies_intermediate(i, :) * 2) / all_delta_A(i);

        % Take the square root, ensuring non-negativity
        mean_sq_velocity_spatial_avg(mean_sq_velocity_spatial_avg < 0) = 0;
        % all_rms_velocity units: sqrt(cm^2/s^2) = cm/s
        all_rms_velocity(i, :) = sqrt(mean_sq_velocity_spatial_avg);

    else
        all_average_energies_intermediate(i, :) = nan(1, numFramesToAnalyze);
        all_rms_velocity(i, :) = nan(1, numFramesToAnalyze);
         if num_points_in_mask(i) == 0 && exist('ME','var') && ~contains(ME.identifier, 'FileNotFound') 
             warning('Number of points in mask is zero for file %s. RMS velocity set to NaN.', matFilename);
         elseif num_points_in_mask(i) == 0 && ~exist('ME','var')
             warning('Number of points in mask is zero for file %s. RMS velocity set to NaN.', matFilename);
         elseif all_delta_A(i) <= 0
             warning('delta_A is zero or negative for file %s (index %d). RMS velocity set to NaN.', matFilename, i);
         end
    end
    if exist('ME','var'), clear ME; end

end



figure('Renderer', 'painters', 'Position', [100, 100, 1000, 600]);
hold on;
smoothing_window = 10;

rms_vel_1f_pro_can_smooth = calculate_mean_smooth(idx_1f_pro_can, all_rms_velocity, smoothing_window);
rms_vel_1f_pro_smo_smooth = calculate_mean_smooth(idx_1f_pro_smo, all_rms_velocity, smoothing_window);
rms_vel_1f_pro_rid_smooth = calculate_mean_smooth(idx_1f_pro_rid, all_rms_velocity, smoothing_window);
rms_vel_1f_retro_can_smooth = calculate_mean_smooth(idx_1f_retro_can, all_rms_velocity, smoothing_window);
rms_vel_1f_retro_smo_smooth = calculate_mean_smooth(idx_1f_retro_smo, all_rms_velocity, smoothing_window);
rms_vel_1f_retro_rid_smooth = calculate_mean_smooth(idx_1f_retro_rid, all_rms_velocity, smoothing_window);
rms_vel_2f_pro_can_smooth = calculate_mean_smooth(idx_2f_pro_can, all_rms_velocity, smoothing_window);
rms_vel_2f_pro_smo_smooth = calculate_mean_smooth(idx_2f_pro_smo, all_rms_velocity, smoothing_window);
rms_vel_2f_pro_rid_smooth = calculate_mean_smooth(idx_2f_pro_rid, all_rms_velocity, smoothing_window);
rms_vel_2f_retro_can_smooth = calculate_mean_smooth(idx_2f_retro_can, all_rms_velocity, smoothing_window);
rms_vel_2f_retro_smo_smooth = calculate_mean_smooth(idx_2f_retro_smo, all_rms_velocity, smoothing_window);
rms_vel_2f_retro_rid_smooth = calculate_mean_smooth(idx_2f_retro_rid, all_rms_velocity, smoothing_window);

% --- Replicate Plotting Logic using rms_vel_* variables ---
plot(xPlot, rms_vel_1f_pro_can_smooth, 'b--', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde canyon');
%plot(xPlot, rms_vel_1f_pro_smo_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde smooth');
%plot(xPlot, rms_vel_1f_pro_rid_smooth, 'b-.', 'LineWidth', 4, 'DisplayName', '+1rpm, Prograde ridge');
plot(xPlot, rms_vel_1f_retro_can_smooth, 'r--', 'LineWidth', 4, 'DisplayName', '-1rpm, Retrograde canyon');
%plot(xPlot, rms_vel_1f_retro_smo_smooth, 'r-', 'LineWidth', 4, 'DisplayName', '-1rpm, Retrograde smooth');
%plot(xPlot, rms_vel_1f_retro_rid_smooth, 'r-.', 'LineWidth', 4, 'DisplayName', '-1rpm, Retrograde ridge');
plot(xPlot, rms_vel_2f_pro_can_smooth, 'b--', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde canyon');
%plot(xPlot, rms_vel_2f_pro_smo_smooth, 'b-', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde smooth');
%plot(xPlot, rms_vel_2f_pro_rid_smooth, 'b--', 'LineWidth', 4, 'DisplayName', '+2rpm, Prograde ridge');
plot(xPlot, rms_vel_2f_retro_can_smooth, 'r--', 'LineWidth', 4, 'DisplayName', '-2rpm, Retrograde canyon');
%plot(xPlot, rms_vel_2f_retro_smo_smooth, 'r-', 'LineWidth', 4, 'DisplayName', '-2rpm, Retrograde smooth');
%plot(xPlot, rms_vel_2f_retro_rid_smooth, 'r--', 'LineWidth', 4, 'DisplayName', '-2rpm, Retrograde ridge');



hold off;
xlabel('Elapsed $t \cdot f$ since start', 'interpreter', 'latex');
ylabel('RMS Velocity [cms^{-1}]');
%title('RMS Velocity');
legend('show', 'Location', 'best');
grid on;
set(gca, 'YMinorGrid', 'on');
% set(gca, 'YScale', 'log');
