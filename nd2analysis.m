%% ND2 File Analysis Script

% Clear workspace and close figures
clear; close all; clc;

%% Setup Bio-Formats MATLAB Toolbox

% Path to the bfmatlab folder (Update this path)
bfmatlabPath = 'C:\Users\mzhang658\OneDrive - Georgia Institute of Technology\Desktop\Phage-biofilm codes_HS\bfmatlab';  % Replace with the actual path

% Add bfmatlab to MATLAB path
addpath(bfmatlabPath);

% Verify Bio-Formats installation
bfCheckJavaPath();

%% Define the Path to Your .nd2 File

% Path to your .nd2 file (Update this path)
nd2_file = 'E:\0908_Ti_0\Ti_0.nd2';  % Replace with the actual path

% Check if the .nd2 file exists
if exist(nd2_file, 'file') ~= 2
    error('ND2 file not found at the specified path:\n%s', nd2_file);
end

%% Process the ND2 File

% Create a Bio-Formats reader
reader = bfGetReader(nd2_file);

% Get the number of positions (series)
numSeries = reader.getSeriesCount();
fprintf('Number of positions (series): %d\n', numSeries);

% Get the number of time points
sizeT = reader.getSizeT();
fprintf('Number of time points: %d\n', sizeT);

% Get number of channels
sizeC = reader.getSizeC();

% Total duration (in hours)
total_time = 15;  % 15 hours

% Time between time points (assuming evenly spaced)
time_interval = total_time / (sizeT - 1);

% Number of groups (every 3 positions form a group)
numGroups = floor(numSeries / 3);

% Initialize arrays to store mean and std areas for each time point, group, and channel
mean_areas_ch1 = zeros(sizeT, numGroups);
std_areas_ch1 = zeros(sizeT, numGroups);
mean_areas_ch2 = zeros(sizeT, numGroups);
std_areas_ch2 = zeros(sizeT, numGroups);
times = zeros(sizeT, 1);

% Process data for each time point
for t = 1:sizeT
    % Assign time for this time point
    times(t) = (t - 1) * time_interval;  % Time in hours

    % Process data for each group
    for g = 1:numGroups
        % Positions in this group
        idx_start = (g - 1) * 3 + 1;
        idx_end = idx_start + 2;
        positions = idx_start:idx_end;
        
        % Initialize lists to collect areas for each channel in this group
        areas_ch1 = [];
        areas_ch2 = [];

        % Process data for each position in the group
        for s = positions
            % Ensure position index does not exceed numSeries
            if s > numSeries
                break;
            end

            % Set the reader to the current series (position)
            reader.setSeries(s - 1);  % Series index starts from 0
            fprintf('Processing position %d/%d at time point %d/%d in group %d/%d\n', s, numSeries, t, sizeT, g, numGroups);

            % Get dimension sizes for this series
            sizeZ = reader.getSizeZ();

            % Process images for each channel
            for c = 1:sizeC
                % Initialize list to collect areas for this channel, position, and time point
                areas = [];
                for z = 5:sizeZ  % Start from 5th slice in z-stack
                    % Get the index of the image plane
                    index = reader.getIndex(z - 1, c - 1, t - 1) + 1;
                    % Read the plane
                    plane = bfGetPlane(reader, index);
                    % Convert to double
                    I = double(plane);
                    % Binarize the image
                    bw = binarize_adapt(I);
                    % Remove small objects
                    bw = bwareaopen(bw, 10);
                    % Label connected components
                    CC = bwconncomp(bw);
                    % Measure areas of connected components
                    stats = regionprops(CC, 'Area');
                    % Collect areas
                    areas = [areas, [stats.Area]];
                end
                % Store areas for each channel
                if c == 1
                    areas_ch1 = [areas_ch1, areas];
                elseif c == 2
                    areas_ch2 = [areas_ch2, areas];
                end
            end
        end

        % Calculate mean and standard deviation for each channel in this group at this time point
        if ~isempty(areas_ch1)
            mean_areas_ch1(t, g) = mean(areas_ch1);
            std_areas_ch1(t, g) = std(areas_ch1);
        else
            mean_areas_ch1(t, g) = 0;
            std_areas_ch1(t, g) = 0;
        end
        if ~isempty(areas_ch2)
            mean_areas_ch2(t, g) = mean(areas_ch2);
            std_areas_ch2(t, g) = std(areas_ch2);
        else
            mean_areas_ch2(t, g) = 0;
            std_areas_ch2(t, g) = 0;
        end
    end
end

% Close the reader
reader.close();

%% Plot the Results for Channel 1

% Plot mean area over time for each group
figure;
hold on;
for g = 1:numGroups
    errorbar(times, mean_areas_ch1(:, g), std_areas_ch1(:, g), '-o', 'LineWidth', 1);
end
xlabel('Time (hours)');
ylabel('Mean Area (\mum^2)');
title('Mean Fluorescent Area Over Time - Channel 1 (Grouped)');
legend(arrayfun(@(x) sprintf('Group %d', x), 1:numGroups, 'UniformOutput', false), 'Location', 'BestOutside');
grid on;
hold off;

%% Plot the Results for Channel 2

% Plot mean area over time for each group
figure;
hold on;
for g = 1:numGroups
    errorbar(times, mean_areas_ch2(:, g), std_areas_ch2(:, g), '-o', 'LineWidth', 1);
end
xlabel('Time (hours)');
ylabel('Mean Area (\mum^2)');
title('Mean Fluorescent Area Over Time - Channel 2 (Grouped)');
legend(arrayfun(@(x) sprintf('Group %d', x), 1:numGroups, 'UniformOutput', false), 'Location', 'BestOutside');
grid on;
hold off;

disp('Processing and plotting complete.');

%% Function Definitions

function bw = binarize_adapt(I)
    % Adjust the binarization parameters if needed
    bw = imbinarize(I, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);
    % Remove small objects (noise)
    bw = bwareaopen(bw, 10);
    % Fill holes in the binary image
    bw = imfill(bw, 'holes');
end
