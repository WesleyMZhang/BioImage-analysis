% User inputs
clear; close all;
parentfolder = 'C:\Users\zhang\GaTech Dropbox\CoS\Physics\curtis_lab\Curtis_Lab Data from Students\Mingxian Zhang\Bacteria_phage_data'; % Replace with the actual path to your data

% Parameters
Ti_values = [0, 4, 8, 12]; % Different Ti values
Bac_values = [6, 7, 8]; % Different CFU values
PFU_values = 0:9; % Different PFU conditions (0 for no phage, 1-9 for 10^n PFU/mL)
timepoints = 1:15; % Different time points (1 to 15)

% Initialize arrays to store results
results = [];

% Process all Ti values
for Ti = Ti_values
    % Construct Ti directory path
    Ti_dir = fullfile(parentfolder, sprintf('Ti = %d', Ti));
    if ~isfolder(Ti_dir)
        disp(['Ti directory does not exist: ', Ti_dir]);
        continue;
    end
    
    for Bac = Bac_values
        % Construct CFU directory path
        CFU_dir = fullfile(Ti_dir, sprintf('CFU_%d', Bac));
        if ~isfolder(CFU_dir)
            disp(['CFU directory does not exist: ', CFU_dir]);
            continue;
        end
        
        for PFU = PFU_values
            % Construct PFU directory path
            if (PFU == 0)
                PFU_dir = fullfile(CFU_dir, 'Control');
            else
                PFU_dir = fullfile(CFU_dir, sprintf('PFU_%d', PFU));
            end
            
            if ~isfolder(PFU_dir)
                disp(['PFU directory does not exist: ', PFU_dir]);
                continue;
            end
            
            % Get point position subfolders
            point_subfolders = dir(PFU_dir);
            point_subfolders = point_subfolders([point_subfolders.isdir]);
            point_subfolders = point_subfolders(~ismember({point_subfolders.name}, {'.', '..'})); % Remove '.' and '..'
            
            for time = timepoints
                % Collect binarization percentages for this time point
                percentages = [];
                for point_subfolder = point_subfolders'
                    subfolder = fullfile(point_subfolder.folder, point_subfolder.name, sprintf('t%02d', time));
                    if isfolder(subfolder)
                        disp(['Processing subfolder: ', subfolder]);
                        [time_val, percentage] = process_subfolder(subfolder);
                        if ~isempty(percentage)
                            percentages = [percentages; percentage];
                        end
                    else
                        disp(['Subfolder does not exist: ', subfolder]);
                    end
                end
                
                % Store results (calculate mean and std)
                if ~isempty(percentages)
                    mean_percentage = mean(percentages);
                    std_percentage = std(percentages);
                    results = [results; Ti, Bac, PFU, time, mean_percentage, std_percentage];
                end
            end
        end
    end
end

% Check if results array is populated
if isempty(results)
    error('No data collected. Please check the folder structure and input parameters.');
end

% Plot results in subplots for each Ti-CFU combination
figure;
t = tiledlayout(numel(Ti_values), numel(Bac_values)); % Create a tiled layout
title(t, 'Binarization Percentage vs. Time for Different Ti, CFU, and PFU Conditions');

for i = 1:numel(Ti_values)
    Ti = Ti_values(i);
    for j = 1:numel(Bac_values)
        Bac = Bac_values(j);
        nexttile;
        hold on;
        colors = lines(numel(PFU_values));
        max_time = 0; % To track the maximum time point for setting x-axis limits
        for PFU = PFU_values
            subset = results(results(:, 1) == Ti & results(:, 2) == Bac & results(:, 3) == PFU, :);
            if ~isempty(subset)
                time_vals = unique(subset(:, 4));
                mean_vals = arrayfun(@(t) mean(subset(subset(:, 4) == t, 5)), time_vals);
                std_vals = arrayfun(@(t) std(subset(subset(:, 4) == t, 5)), time_vals);
                sem_vals = std_vals ./ sqrt(arrayfun(@(t) sum(subset(:, 4) == t), time_vals));
                errorbar(time_vals, mean_vals, sem_vals, 'Color', colors(PFU + 1, :), 'DisplayName', sprintf('PFU = %d', PFU));
                max_time = max(max_time, max(time_vals)); % Update maximum time point
            else
                warning('No data for Ti = %d, CFU = %d, PFU = %d', Ti, Bac, PFU);
            end
        end

        xlabel('Time');
        ylabel('Binarization Percentage (%)');
        title(sprintf('Ti = %d, CFU = %d', Ti, Bac));
        legend('show');
        grid on;
        xlim([0 max_time]); % Set x-axis limits from 0 to the maximum time point
        hold off;
    end
end

disp('Plotting complete.');

% Functions

function percentage = calculate_binarization_percentage(binary_image)
    percentage = sum(binary_image(:)) / numel(binary_image) * 100;
end

function [time, percentage] = process_subfolder(subfolder)
    MyFolderInfo = dir(subfolder);
    MyFolderInfo = MyFolderInfo(~ismember({MyFolderInfo.name}, {'.', '..'}));
    
    if isempty(MyFolderInfo)
        warning(['No images found in subfolder: ', subfolder]);
        time = [];
        percentage = [];
        return;
    end
    
    str = MyFolderInfo(end).name;
    str = convertCharsToStrings(str);

    num_planes = str2double(extractBetween(str, "z", "c"));

    percentages = zeros(num_planes, 1);
    for n = 1:min(num_planes, 10) % Process only the first 10 images
        try
            % Access the SYTO9 channel images (odd indices)
            image_file_g = fullfile(subfolder, MyFolderInfo(2*n-1).name);
            I_g = imread(image_file_g);
            
            % Binarize the image
            I_g_bn = binarize_adapt(I_g);
            
            % Calculate the binarization percentage
            percentages(n) = calculate_binarization_percentage(I_g_bn);
        catch ME
            warning('Error processing image at plane %d: %s', n, ME.message);
        end
    end

    % Use the maximum percentage in the first 10 images for this time point
    [max_percentage, idx] = max(percentages);
    if idx <= num_planes
        time = str2double(extractAfter(subfolder, 't')); % Assuming subfolder names like 't01', 't02', etc.
        percentage = max_percentage;
    else
        time = [];
        percentage = [];
    end
end

function binary_image = binarize_adapt(image)
    % Binarize a 16-bit image using adaptive thresholding
    % Convert the image into 16-bit unsigned integer data
    image = uint16(image);
    % Apply adaptive thresholding with a sensitivity of 0.4
    T = adaptthresh(image, 0.4, 'Statistic', 'mean');
    % Binarize the image using the adaptive threshold filter
    BW = imbinarize(image, T);
    % Perform morphological operations tuned to recognize bacteria (longitudinal and cross-sections)
    se = strel('disk', 2);
    BW = imclose(BW, se);
    BW = bwareaopen(BW, 12);
    binary_image = BW;
end
