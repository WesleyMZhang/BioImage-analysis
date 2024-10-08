% User inputs
clear; close all;
parentfolder = 'C:\Users\zhang\GaTech Dropbox\CoS\Physics\curtis_lab\Curtis_Lab Data from Students\Mingxian Zhang\Bacteria_phage_data\Ti = 12\CFU_7\PFU_4\xy44\t15'; % Add the path to the folder with raw data


% Populate the following variables to access the corresponding data
ti = 2; % Enter 1 for Ti0, 2 for Ti4, 3 for Ti8, 4 for Ti12
b = 3; % Enter 1 for 10^6 CFU/mL, 2 for 10^7 CFU/mL, 3 for 10^8 CFU/mL 
p = 0;  % Enter 0 for no phage control, n = 1 or 2 or .. or 9 for 10^n PFU/mL 
r = 3;  % Enter r = 1 or 2 or 3 for the corresponding replicate
t = 1;  % Enter t = 1 or 2 or .. or 15 for the corresponding timepoint

Ti = [0, 4, 8, 12];
Bac = [6, 7, 8];

% Access the relevant folder containing image data
stackdir = parentfolder;

MyFolderInfo = dir(stackdir);
MyFolderInfo = MyFolderInfo(~ismember({MyFolderInfo.name}, {'.', '..'})); % Remove the first two rows
str = MyFolderInfo(end).name;
str = convertCharsToStrings(str);

% Identify relevant parameters
num_channels = str2double(extractBetween(str, "c", ".tif"));
num_planes = str2double(extractBetween(str, "z", "c"));

% Preallocate a 3D array to store the image stack
info = imfinfo(fullfile(stackdir, MyFolderInfo(1).name));
imgHeight = info.Height;
imgWidth = info.Width;
imageStack = zeros(imgHeight, imgWidth, num_planes, 'uint16');
binarizedStack = zeros(imgHeight, imgWidth, num_planes, 'logical');

% Load images into the 3D array
for n = 1:num_planes
    try
        % Access the SYTO9 channel images (odd indices)
        image_file_g = fullfile(stackdir, MyFolderInfo(2*n-1).name);
        I_g = imread(image_file_g);
        
        % Store images in the stack
        imageStack(:, :, n) = I_g;  % Only store odd index, channel 1 image
        
        % Binarize the image
        I_g_bn = binarize_adapt(I_g);
        binarizedStack(:, :, n) = I_g_bn;
    catch ME
        warning('Error loading image at plane %d: %s', n, ME.message);
    end
end

% Create a figure and axes for displaying the images
hFig = figure('Name', 'Z-Stack Viewer - Odd Channels Only', 'NumberTitle', 'off');
hAx1 = subplot(1, 2, 1, 'Parent', hFig);
hAx2 = subplot(1, 2, 2, 'Parent', hFig);

% Display the first image with adjusted contrast
hImg1 = imagesc(imadjust(imageStack(:, :, 1)), 'Parent', hAx1);
hImg2 = imagesc(mat2gray(binarizedStack(:, :, 1), [0, 1]), 'Parent', hAx2);
axis(hAx1, 'off');
axis(hAx2, 'off');
title(hAx1, 'Original Image');
title(hAx2, 'Binarized Image');

% Create a slider for navigating through the stack
hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num_planes, ...
    'Value', 1, 'SliderStep', [1/(num_planes - 1), 1/(num_planes - 1)], ...
    'Position', [100, 50, 400, 20], 'Callback', @(src, event) updateImage(src, hImg1, hImg2, imageStack, binarizedStack, hAx1, hAx2));

disp('Use the slider to scroll through the odd-order channel images.');

% Function to update the displayed image based on slider value
function updateImage(hSlider, hImg1, hImg2, imageStack, binarizedStack, hAx1, hAx2)
    idx = round(get(hSlider, 'Value'));  % Get the current slider value
    I_fig = imageStack(:, :, idx);
    I_fig_bn = binarizedStack(:, :, idx);
    set(hImg1, 'CData', imadjust(I_fig));  % Update the original image data with adjusted contrast
    set(hImg2, 'CData', mat2gray(I_fig_bn, [0, 1]));  % Update the binarized image data
    title(hAx1, sprintf('Original Image - Plane %d', idx));
    title(hAx2, sprintf('Binarized Image - Plane %d', idx));
    drawnow;
end


% Define the binarization function
function binary_image = binarize_adapt(image)
    % Binarize a 16-bit image using adaptive thresholding
    % Convert the image into 16-bit unsigned integer data
    image = uint16(image);
    % Apply adaptive thresholding with a sensitivity of 0.4
    T = adaptthresh(image, 0.4, 'Statistic', 'mean');
    % Binarize the image using the adaptive threshold filter
    BW = imbinarize(image, T);
    % Perform morphological operations tuned to recognize bacteria (longitudinal
    % and cross-sections)
    se = strel('disk', 2);
    BW = imclose(BW, se);
    BW = bwareaopen(BW, 12);
    binary_image = BW;
end