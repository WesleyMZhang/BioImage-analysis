% Add Bio-Formats library to the MATLAB path
addpath('C:\Users\mzhang658\OneDrive - Georgia Institute of Technology\Desktop\Phage-biofilm codes_HS\bfmatlab');  % Adjust the path as needed

% Define the .nd2 file path
nd2_filepath = 'E:\01_11_ Ti_0\Ti_0.nd2';  % Adjust the path as needed

% Initialize the reader
try
    reader = bfGetReader(nd2_filepath);
catch ME
    fprintf('Error initializing bfGetReader: %s\n', ME.message);
    return;  % Exit the script if the file cannot be opened
end

% Get total series count
seriesCount = reader.getSeriesCount();
fprintf('Total number of series in the file: %d\n', seriesCount);

% Option 1: Input the series range to analyze interactively
series_range = input('Enter the series range to analyze (e.g., [1 30] for series 1 to 30): ');

% Option 2: Alternatively, define the series range directly here
% Uncomment the line below if you want to skip the interactive input
% series_range = [1 30];  % Adjust the range as needed

% Validate the input
if numel(series_range) ~= 2 || series_range(1) < 1 || series_range(2) > seriesCount
    error('Invalid series range. Ensure that it is within [1, %d].', seriesCount);
end

% Define starting planes to test
starting_planes = [1, 5, 6];  % Starting from 1st, 5th, and 6th planes

% Prompt the user to input sensitivity values for each channel
fprintf('\nEnter binarization sensitivity values for each channel (between 0 and 1):\n');
sensitivity_channel_1 = input('Sensitivity for Channel 1 (e.g., 0.4): ');
sensitivity_channel_2 = input('Sensitivity for Channel 2 (e.g., 0.4): ');

% Validate sensitivity values
if sensitivity_channel_1 <= 0 || sensitivity_channel_1 >= 1
    error('Invalid sensitivity value for Channel 1. It must be between 0 and 1.');
end
if sensitivity_channel_2 <= 0 || sensitivity_channel_2 >= 1
    error('Invalid sensitivity value for Channel 2. It must be between 0 and 1.');
end

% Initialize data structures to store results for each starting plane
numSeries = series_range(2) - series_range(1) + 1;
max_time_points = 0;  % To keep track of the maximum number of time points

% Initialize cell arrays to store results per series for each starting plane
counts_diff_per_series = cell(numSeries, length(starting_planes));
bacteria_count_channel_2_per_series = cell(numSeries, length(starting_planes));
bacteria_count_channel_1_per_series = cell(numSeries, length(starting_planes));  % New cell array for Channel 1 counts

% Process each series in the specified range
idx = 0;
for seriesIdx = series_range(1):series_range(2)
    idx = idx + 1;
    fprintf('\nProcessing Series %d of %d\n', seriesIdx, seriesCount);

    try
        % Set the current series (0-indexed)
        reader.setSeries(seriesIdx - 1);
    catch ME
        warning('Error setting series %d: %s. Skipping...', seriesIdx, ME.message);
        continue;  % Skip this series and continue with the next one
    end

    % Get metadata for this series
    omeMeta = reader.getMetadataStore();  % Access metadata store
    try
        z_size = omeMeta.getPixelsSizeZ(seriesIdx - 1).getValue();  % Number of z-stacks
        t_size = omeMeta.getPixelsSizeT(seriesIdx - 1).getValue();  % Number of time points
        channel_count = omeMeta.getChannelCount(seriesIdx - 1);     % Number of channels (e.g., 2 channels)
    catch ME
        warning('Error retrieving metadata for series %d: %s. Skipping...', seriesIdx, ME.message);
        continue;  % Skip this series and continue with the next one
    end

    fprintf('Series %d - Z size: %d, T size: %d, Channels: %d\n', seriesIdx, z_size, t_size, channel_count);

    % Update max_time_points
    if t_size > max_time_points
        max_time_points = t_size;
    end

    if z_size >= min(starting_planes)
        % Initialize arrays to store counts per time point, per z-plane
        counts_channel_1_z = zeros(t_size, z_size);
        counts_channel_2_z = zeros(t_size, z_size);

        % Initialize counters for successful and failed planes
        successful_planes = 0;
        failed_planes = 0;

        for timeIdx = 1:t_size
            % Process each z-plane
            for zIdx = 1:z_size
                try
                    % Compute plane indices for both channels
                    % Planes are 0-indexed in Bio-Formats
                    planeIdx_channel_1 = reader.getIndex(zIdx - 1, 0, timeIdx - 1) + 1;  % Channel 1
                    planeIdx_channel_2 = reader.getIndex(zIdx - 1, 1, timeIdx - 1) + 1;  % Channel 2

                    % Load the image for both channels
                    I_channel_1 = bfGetPlane(reader, planeIdx_channel_1);
                    I_channel_2 = bfGetPlane(reader, planeIdx_channel_2);

                    % Normalize images if necessary
                    if isa(I_channel_1, 'uint16') || isa(I_channel_1, 'uint32')
                        I_channel_1 = mat2gray(I_channel_1);  % Convert to double in [0,1]
                    end
                    if isa(I_channel_2, 'uint16') || isa(I_channel_2, 'uint32')
                        I_channel_2 = mat2gray(I_channel_2);
                    end

                    % Binarize and count bacteria for both channels using specified sensitivities
                    I_channel_1_bn = binarize_adapt(I_channel_1, sensitivity_channel_1);
                    I_channel_2_bn = binarize_adapt(I_channel_2, sensitivity_channel_2);

                    % Count the non-zero pixels (bacteria count) for both channels
                    counts_channel_1_z(timeIdx, zIdx) = nnz(I_channel_1_bn);
                    counts_channel_2_z(timeIdx, zIdx) = nnz(I_channel_2_bn);

                    % Increment successful planes counter
                    successful_planes = successful_planes + 1;
                catch ME
                    % Increment failed planes counter
                    failed_planes = failed_planes + 1;
                    warning('Error processing z-plane %d for series %d: %s', zIdx, seriesIdx, ME.message);
                    % Continue to next z-plane
                end
            end
        end

        % For each starting plane, compute cumulative counts and store results
        for sp_idx = 1:length(starting_planes)
            start_plane = starting_planes(sp_idx);
            if z_size >= start_plane
                % Sum counts from starting plane to the last plane
                sum_counts_channel_1 = sum(counts_channel_1_z(:, start_plane:end), 2);
                sum_counts_channel_2 = sum(counts_channel_2_z(:, start_plane:end), 2);

                % Compute the difference for each time point
                counts_diff = sum_counts_channel_1 - sum_counts_channel_2;
                bacteria_count_channel_2 = sum_counts_channel_2;
                bacteria_count_channel_1 = sum_counts_channel_1;  % Store Channel 1 counts

                % Store the counts difference and Channel counts for this series and starting plane
                counts_diff_per_series{idx, sp_idx} = counts_diff;
                bacteria_count_channel_2_per_series{idx, sp_idx} = bacteria_count_channel_2;
                bacteria_count_channel_1_per_series{idx, sp_idx} = bacteria_count_channel_1;
            else
                warning('Series %d has fewer than %d z-planes, skipping starting plane %d...', seriesIdx, start_plane, start_plane);
                counts_diff_per_series{idx, sp_idx} = [];
                bacteria_count_channel_2_per_series{idx, sp_idx} = [];
                bacteria_count_channel_1_per_series{idx, sp_idx} = [];
            end
        end

        % Print summary of successful and failed planes
        fprintf('Series %d: %d successful planes, %d failed planes.\n', seriesIdx, successful_planes, failed_planes);
    else
        warning('Series %d has fewer than %d z-planes, skipping...', seriesIdx, min(starting_planes));
    end
end

% Close the reader after processing
reader.close();

%% Plotting the results for different starting planes

% Determine the number of groups (every three series)
group_size = 3;
num_groups = ceil(numSeries / group_size);

% Generate a colormap with enough distinct colors using 'distinguishable_colors'
if exist('distinguishable_colors', 'file')
    colors = distinguishable_colors(num_groups);
else
    % Alternative: use hsv colormap if distinguishable_colors is not available
    colors = hsv(num_groups);
end

% Custom legend labels
% Assuming you have up to 10 groups: 1 Control + 9 PFU sets
legend_labels = {'Control Set', 'PFU_1', 'PFU_2', 'PFU_3', 'PFU_4', 'PFU_5', 'PFU_6', 'PFU_7', 'PFU_8', 'PFU_9'};

% Time indices
time_indices = 1:max_time_points;

% Plotting
figure;

for sp_idx = 1:length(starting_planes)
    start_plane = starting_planes(sp_idx);

    % Initialize matrices to collect counts difference and Channel counts
    group_counts_diff = cell(num_groups, 1);
    group_counts_channel_2 = cell(num_groups, 1);
    group_counts_channel_1 = cell(num_groups, 1);

    % Loop over each group
    for group_idx = 1:num_groups
        % Determine the indices of the series in this group
        start_idx = (group_idx - 1) * group_size + 1;
        end_idx = min(group_idx * group_size, numSeries);

        % Initialize
        group_counts_diff{group_idx} = [];
        group_counts_channel_2{group_idx} = [];
        group_counts_channel_1{group_idx} = [];

        % Collect data from series in the group
        for series_i = start_idx:end_idx
            counts_diff = counts_diff_per_series{series_i, sp_idx};
            counts_channel_2 = bacteria_count_channel_2_per_series{series_i, sp_idx};
            counts_channel_1 = bacteria_count_channel_1_per_series{series_i, sp_idx};

            % Check if data exists for this series and starting plane
            if isempty(counts_diff)
                continue;
            end

            % Pad counts with NaNs if necessary
            if length(counts_diff) < max_time_points
                pad_len = max_time_points - length(counts_diff);
                counts_diff = [counts_diff; NaN(pad_len, 1)];
                counts_channel_2 = [counts_channel_2; NaN(pad_len, 1)];
                counts_channel_1 = [counts_channel_1; NaN(pad_len, 1)];
            end

            % Collect data
            group_counts_diff{group_idx} = [group_counts_diff{group_idx}, counts_diff];
            group_counts_channel_2{group_idx} = [group_counts_channel_2{group_idx}, counts_channel_2];
            group_counts_channel_1{group_idx} = [group_counts_channel_1{group_idx}, counts_channel_1];
        end
    end

    % Subplot arrangement
    % 1st column: live bacteria counting (Counts Diff)
    % 2nd column: Propidium Iodide plotting (Channel 2)
    % 3rd column: SYTO9 plotting (Channel 1)
    subplot_idx_diff = (sp_idx - 1) * 3 + 1;
    subplot_idx_channel_2 = subplot_idx_diff + 1;
    subplot_idx_channel_1 = subplot_idx_diff + 2;

    % Plot Mean Counts Difference (live bacteria counting)
    subplot(length(starting_planes), 3, subplot_idx_diff);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_diff{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series_in_group = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series_in_group);
        legend_label = legend_labels{min(group_idx, length(legend_labels))}; 
        
        % Set color: Control Set (group_idx=1) is black, others use colors array
        if group_idx == 1
            line_color = 'k';
        else
            line_color = colors(group_idx, :);
        end

        errorbar(time_indices, mean_data, sem_data, '-o', 'Color', line_color, 'DisplayName', legend_label);
    end
    xlabel('Time Point');
    ylabel('Mean Bacteria Count Difference (Ch1 - Ch2)');
    title('live bacteria counting');
    legend('show', 'Location', 'bestoutside');
    grid on;
    xlim([1, max_time_points]);
    xticks(1:max_time_points);
    ylim([0 4e7]);  % Fixed y-axis range

    % Plot Mean Channel 2 Counts (Propidium Iodide plotting)
    subplot(length(starting_planes), 3, subplot_idx_channel_2);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_channel_2{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series_in_group = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series_in_group);
        legend_label = legend_labels{min(group_idx, length(legend_labels))};
        
        % Set color: Control Set (group_idx=1) is black, others use colors array
        if group_idx == 1
            line_color = 'k';
        else
            line_color = colors(group_idx, :);
        end

        errorbar(time_indices, mean_data, sem_data, '-s', 'Color', line_color, 'DisplayName', legend_label);
    end
    xlabel('Time Point');
    ylabel('Mean Bacteria Count (Channel 2)');
    title('Propidium Iodide plotting');
    legend('show', 'Location', 'bestoutside');
    grid on;
    xlim([1, max_time_points]);
    xticks(1:max_time_points);
    ylim([0 4e7]);  % Fixed y-axis range

    % Plot Mean Channel 1 Counts (SYTO9 plotting)
    subplot(length(starting_planes), 3, subplot_idx_channel_1);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_channel_1{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series_in_group = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series_in_group);
        legend_label = legend_labels{min(group_idx, length(legend_labels))};
        
        % Set color: Control Set (group_idx=1) is black, others use colors array
        if group_idx == 1
            line_color = 'k';
        else
            line_color = colors(group_idx, :);
        end

        errorbar(time_indices, mean_data, sem_data, '-d', 'Color', line_color, 'DisplayName', legend_label);
    end
    xlabel('Time Point');
    ylabel('Mean Bacteria Count (Channel 1)');
    title('SYTO9 plotting');
    legend('show', 'Location', 'bestoutside');
    grid on;
    xlim([1, max_time_points]);
    xticks(1:max_time_points);
    ylim([0 4e7]);  % Fixed y-axis range
end

disp('Bacteria counting and plotting complete.');

% Save data to .mat file for easy re-plotting
save('plot_data.mat', 'counts_diff_per_series', 'bacteria_count_channel_1_per_series', 'bacteria_count_channel_2_per_series', 'time_indices', 'starting_planes', 'series_range', 'sensitivity_channel_1', 'sensitivity_channel_2');

%% Functions

function bw = binarize_adapt(I, sensitivity)
    % Normalize the image if necessary
    if isa(I, 'uint16') || isa(I, 'uint32')
        I = mat2gray(I);  % Convert to double in [0,1]
    end

    % Validate sensitivity value
    if sensitivity <= 0 || sensitivity >= 1
        error('Sensitivity must be between 0 and 1.');
    end

    % Apply adaptive thresholding with the specified sensitivity
    bw = imbinarize(I, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', sensitivity);
end
