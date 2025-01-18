% Load the saved data
load('C:\Users\mzhang658\OneDrive - Georgia Institute of Technology\Desktop\Phage-biofilm codes_HS\Ti_0_CFU_6.mat');

% Variables loaded include:
% counts_diff_per_series, bacteria_count_channel_1_per_series, bacteria_count_channel_2_per_series
% time_indices, starting_planes, series_range, sensitivity_channel_1, sensitivity_channel_2

% Determine number of series and groups
numSeries = series_range(2) - series_range(1) + 1;
group_size = 3;
num_groups = ceil(numSeries / group_size);

% Define custom legend labels
legend_labels = {'Control Set', 'PFU_1', 'PFU_2', 'PFU_3', 'PFU_4', 'PFU_5', 'PFU_6', 'PFU_7', 'PFU_8', 'PFU_9'};

% Generate colors
if exist('distinguishable_colors', 'file')
    colors = distinguishable_colors(num_groups);
else
    colors = hsv(num_groups);
end

% Create a new figure for plotting
figure;

% Loop over each starting plane
for sp_idx = 1:length(starting_planes)
    start_plane = starting_planes(sp_idx);
    
    % Initialize matrices to group data
    group_counts_diff = cell(num_groups, 1);
    group_counts_channel_2 = cell(num_groups, 1);
    group_counts_channel_1 = cell(num_groups, 1);
    
    % Collect data for each group
    for group_idx = 1:num_groups
        start_idx = (group_idx - 1) * group_size + 1;
        end_idx = min(group_idx * group_size, numSeries);

        group_counts_diff{group_idx} = [];
        group_counts_channel_2{group_idx} = [];
        group_counts_channel_1{group_idx} = [];
        
        for series_i = start_idx:end_idx
            counts_diff = counts_diff_per_series{series_i, sp_idx};
            counts_channel_2 = bacteria_count_channel_2_per_series{series_i, sp_idx};
            counts_channel_1 = bacteria_count_channel_1_per_series{series_i, sp_idx};

            if isempty(counts_diff)
                continue;
            end
            
            % Ensure all data vectors match the max time length
            max_time_points = length(time_indices);
            if length(counts_diff) < max_time_points
                pad_len = max_time_points - length(counts_diff);
                counts_diff = [counts_diff; NaN(pad_len, 1)];
                counts_channel_2 = [counts_channel_2; NaN(pad_len, 1)];
                counts_channel_1 = [counts_channel_1; NaN(pad_len, 1)];
            end

            group_counts_diff{group_idx} = [group_counts_diff{group_idx}, counts_diff];
            group_counts_channel_2{group_idx} = [group_counts_channel_2{group_idx}, counts_channel_2];
            group_counts_channel_1{group_idx} = [group_counts_channel_1{group_idx}, counts_channel_1];
        end
    end

    % Subplot indices
    subplot_idx_diff = (sp_idx - 1)*3 + 1;
    subplot_idx_channel_2 = subplot_idx_diff + 1;
    subplot_idx_channel_1 = subplot_idx_diff + 2;

    % Plot Live Bacteria Counting (Counts Diff)
    subplot(length(starting_planes), 3, subplot_idx_diff);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_diff{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series);
        legend_label = legend_labels{min(group_idx, length(legend_labels))};
        
        % Control Set is black, others use the colors array
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
    xlim([1, max(time_indices)]);
    xticks(1:max(time_indices));
    ylim([-1e7 4e7]); 

    % Plot Propidium Iodide (Channel 2)
    subplot(length(starting_planes), 3, subplot_idx_channel_2);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_channel_2{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series);
        legend_label = legend_labels{min(group_idx, length(legend_labels))};
        
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
    xlim([1, max(time_indices)]);
    xticks(1:max(time_indices));
    ylim([-1e7 4e7]);

    % Plot SYTO9 (Channel 1)
    subplot(length(starting_planes), 3, subplot_idx_channel_1);
    hold on;
    for group_idx = 1:num_groups
        group_data = group_counts_channel_1{group_idx};
        if isempty(group_data), continue; end
        mean_data = mean(group_data, 2, 'omitnan');
        std_data = std(group_data, 0, 2, 'omitnan');
        n_series = size(group_data, 2);
        sem_data = std_data ./ sqrt(n_series);
        legend_label = legend_labels{min(group_idx, length(legend_labels))};
        
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
    xlim([1, max(time_indices)]);
    xticks(1:max(time_indices));
    ylim([-1e7 4e7]);

end
