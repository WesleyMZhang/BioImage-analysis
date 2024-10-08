clear; close all;

% Define the parent folder
parentfolder = 'C:\Users\zhang\Dropbox (GaTech)\curtis_lab\Curtis_Lab Data from Students\Mingxian Zhang\Bacteria_phage_data\Ti = 0\CFU_7\Control';

% Get the list of PFU subfolders
pfu_folders = dir(parentfolder);
pfu_folders = pfu_folders([pfu_folders.isdir] & ~ismember({pfu_folders.name}, {'.', '..'}));

% Initialize arrays to store time and count data
times = [];
counts = [];

% Process each PFU subfolder
for p = 1:length(pfu_folders)
    pfu_folder = fullfile(parentfolder, pfu_folders(p).name);
    fprintf('Processing PFU folder: %s\n', pfu_folder);  % Debugging output

    % Get the list of time subfolders
    time_folders = dir(pfu_folder);
    time_folders = time_folders([time_folders.isdir] & ~ismember({time_folders.name}, {'.', '..'}));
    
    % Process each time subfolder
    for t = 1:length(time_folders)
        time_folder = fullfile(pfu_folder, time_folders(t).name);
        fprintf('Processing time folder: %s\n', time_folder);  % Debugging output
        try
            [time, count] = process_subfolder(time_folder);
            fprintf('Time: %d, Count: %d\n', time, count);  % Debugging output
            times = [times; time];
            counts = [counts; count];
        catch ME
            warning('Error processing folder %s: %s', time_folder, ME.message);
        end
    end
end

% Check if data is collected
if isempty(times) || isempty(counts)
    error('No data collected. Please check the folder structure and image files.');
end

% Calculate average counts per time point
unique_times = unique(times);
avg_counts = arrayfun(@(t) mean(counts(times == t)), unique_times);

% Plot the results
figure;
plot(unique_times, avg_counts, '-o');
xlabel('Time');
ylabel('Average Live Bacteria Count');
title('Average Live Bacteria Count Over Time');
grid on;

disp('Processing and plotting complete.');


% Define the binarization function (example using adaptive thresholding)
function bw = binarize_adapt(I)
    bw = imbinarize(I, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);
end

% Function to process and count live bacteria in a given subfolder
function [time, count] = process_subfolder(subfolder)
    MyFolderInfo = dir(fullfile(subfolder, '*.tif'));
    
    if isempty(MyFolderInfo)
        error('No .tif files found in the folder: %s', subfolder);
    end
    
    str = MyFolderInfo(end).name;
    str = convertCharsToStrings(str);

    % Identify relevant parameters
    num_channels = str2double(extractBetween(str, "c", ".tif"));
    num_planes = str2double(extractBetween(str, "z", "c"));

    % Initialize count
    live_bacteria_count = 0;

    % Process images
    for n = 1:num_planes
        try
            % Access the SYTO9 channel images (odd indices)
            if mod(n, 2) == 1
                image_file_g = fullfile(subfolder, MyFolderInfo(n).name);
                I_g = imread(image_file_g);

                % Binarize the image
                I_g_bn = binarize_adapt(I_g);

                % Count live bacteria
                live_bacteria_count = live_bacteria_count + sum(I_g_bn(:));
            end
        catch ME
            warning('Error loading image at plane %d: %s', n, ME.message);
        end
    end

    % Extract time information from the folder name
    [~, folder_name, ~] = fileparts(subfolder);
    time = str2double(extractAfter(folder_name, 't'));

    % Return the time and count
    count = live_bacteria_count;
end
