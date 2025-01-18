% Add Bio-Formats library to the MATLAB path
addpath('C:\Users\mzhang658\OneDrive - Georgia Institute of Technology\Desktop\Phage-biofilm codes_HS\bfmatlab');  % Adjust the path to where bfmatlab is located

% Define the .nd2 file path
nd2_filepath = 'E:\01_11_ Ti_0\Ti_0.nd2';  % Adjust the file path as needed

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

% Prompt the user to select a series
seriesIdx = input(sprintf('Enter the series number to analyze (1 to %d): ', seriesCount));

% Validate the input
if seriesIdx < 1 || seriesIdx > seriesCount
    error('Invalid series index. Ensure that it is within [1, %d].', seriesCount);
end

% Set the current series (0-indexed)
reader.setSeries(seriesIdx - 1);

% Get metadata for this series
omeMeta = reader.getMetadataStore();  % Access metadata store
try
    z_size = omeMeta.getPixelsSizeZ(seriesIdx - 1).getValue();  % Number of z-stacks
    t_size = omeMeta.getPixelsSizeT(seriesIdx - 1).getValue();  % Number of time points
    channel_count = omeMeta.getChannelCount(seriesIdx - 1);     % Number of channels
catch ME
    warning('Error retrieving metadata for series %d: %s.', seriesIdx, ME.message);
    reader.close();
    return;
end

fprintf('Series %d - Z size: %d, T size: %d, Channels: %d\n', seriesIdx, z_size, t_size, channel_count);

% Prompt the user to select a time point
if t_size > 1
    timeIdx = input(sprintf('Enter the time point to analyze (1 to %d): ', t_size));
else
    timeIdx = 1;
end

% Validate the time index
if timeIdx < 1 || timeIdx > t_size
    error('Invalid time index. Ensure that it is within [1, %d].', t_size);
end

% Sensitivity values for each channel
% Adjust the sensitivity values as needed
sensitivity_values = zeros(1, channel_count);
sensitivity_values(1) = 0.05;  % Channel 1 (e.g., Green)
sensitivity_values(2) = 0.01;  % Channel 2 (e.g., Red)

% Initialize arrays to store the original and binarized images
z_stack_original = cell(z_size, channel_count);
z_stack_binarized = cell(z_size, channel_count);

% Loop through the z-stack for all channels
for chIdx = 1:channel_count
    for zIdx = 1:z_size
        try
            % Compute the plane index
            planeIdx = reader.getIndex(zIdx - 1, chIdx - 1, timeIdx - 1) + 1;

            % Load the image
            I = bfGetPlane(reader, planeIdx);

            % Normalize the image if necessary
            if isa(I, 'uint16') || isa(I, 'uint32')
                I_norm = mat2gray(I);  % Convert to double in [0,1]
            else
                I_norm = I;
            end

            % Apply binarization with the specified sensitivity
            I_bin = binarize_adapt(I_norm, sensitivity_values(chIdx));

            % Store the original and binarized images
            z_stack_original{zIdx, chIdx} = I_norm;
            z_stack_binarized{zIdx, chIdx} = I_bin;
        catch ME
            warning('Error processing z-plane %d, channel %d: %s', zIdx, chIdx, ME.message);
            z_stack_original{zIdx, chIdx} = [];
            z_stack_binarized{zIdx, chIdx} = [];
        end
    end
end

% Close the reader after processing
reader.close();

%% Create an interactive GUI with sliders
% Convert cell arrays to 4D matrices for easier indexing
% Dimensions: [height, width, z, channels]
% Get image dimensions from a sample image
sample_idx = find(~cellfun(@isempty, z_stack_original(:, 1)), 1);
if isempty(sample_idx)
    error('No images found to display.');
end
sample_image = z_stack_original{sample_idx, 1};
[img_height, img_width] = size(sample_image);

% Initialize matrices
orig_stack = zeros(img_height, img_width, z_size, channel_count);
bin_stack = zeros(img_height, img_width, z_size, channel_count);

for chIdx = 1:channel_count
    for zIdx = 1:z_size
        if ~isempty(z_stack_original{zIdx, chIdx})
            orig_stack(:, :, zIdx, chIdx) = z_stack_original{zIdx, chIdx};
            bin_stack(:, :, zIdx, chIdx) = z_stack_binarized{zIdx, chIdx};
        else
            orig_stack(:, :, zIdx, chIdx) = zeros(img_height, img_width);
            bin_stack(:, :, zIdx, chIdx) = zeros(img_height, img_width);
        end
    end
end

% Initialize the figure and display
fig = figure('Name', 'Binarized Z-stack Viewer', 'NumberTitle', 'off');

% Initial z-plane and channel indices
current_z = 1;
current_channel = 1;

% Store variables in the figure's app data for access in callbacks
setappdata(fig, 'orig_stack', orig_stack);
setappdata(fig, 'bin_stack', bin_stack);
setappdata(fig, 'z_size', z_size);
setappdata(fig, 'channel_count', channel_count);
setappdata(fig, 'current_z', current_z);
setappdata(fig, 'current_channel', current_channel);

% Display the initial images
subplot(1, 2, 1);
hOrigImage = imshow(orig_stack(:, :, current_z, current_channel), []);
title(sprintf('Original - Z-plane %d of %d, Channel %d of %d', current_z, z_size, current_channel, channel_count));

subplot(1, 2, 2);
hBinImage = imshow(bin_stack(:, :, current_z, current_channel), []);
title(sprintf('Binarized - Z-plane %d of %d, Channel %d of %d', current_z, z_size, current_channel, channel_count));

% Create sliders for z-plane and channel selection
% Slider for z-plane
hZSlider = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', z_size, 'Value', current_z, ...
    'SliderStep', [1/(z_size - 1), 1/(z_size - 1)], ...
    'Position', [150 20 300 20], ...
    'Callback', @zSliderCallback);

% Slider for channel
hChannelSlider = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', channel_count, 'Value', current_channel, ...
    'SliderStep', [1/(channel_count - 1), 1/(channel_count - 1)], ...
    'Position', [150 60 300 20], ...
    'Callback', @channelSliderCallback);

% Labels for sliders
uicontrol('Style', 'text', 'Position', [50 20 80 20], 'String', 'Z-plane');
uicontrol('Style', 'text', 'Position', [50 60 80 20], 'String', 'Channel');

% Store handles in the figure's app data
setappdata(fig, 'hOrigImage', hOrigImage);
setappdata(fig, 'hBinImage', hBinImage);
setappdata(fig, 'hZSlider', hZSlider);
setappdata(fig, 'hChannelSlider', hChannelSlider);

% Callback function for z-plane slider
function zSliderCallback(hObj, ~)
    fig = ancestor(hObj, 'figure');
    current_z = round(get(hObj, 'Value'));
    setappdata(fig, 'current_z', current_z);
    updateImage(fig);
end

% Callback function for channel slider
function channelSliderCallback(hObj, ~)
    fig = ancestor(hObj, 'figure');
    current_channel = round(get(hObj, 'Value'));
    setappdata(fig, 'current_channel', current_channel);
    updateImage(fig);
end

% Function to update the displayed images
function updateImage(fig)
    % Retrieve variables from app data
    orig_stack = getappdata(fig, 'orig_stack');
    bin_stack = getappdata(fig, 'bin_stack');
    current_z = getappdata(fig, 'current_z');
    current_channel = getappdata(fig, 'current_channel');
    z_size = getappdata(fig, 'z_size');
    channel_count = getappdata(fig, 'channel_count');
    hOrigImage = getappdata(fig, 'hOrigImage');
    hBinImage = getappdata(fig, 'hBinImage');
    
    % Update the original image
    set(hOrigImage, 'CData', orig_stack(:, :, current_z, current_channel));
    subplot(1, 2, 1);
    title(sprintf('Original - Z-plane %d of %d, Channel %d of %d', current_z, z_size, current_channel, channel_count));
    
    % Update the binarized image
    set(hBinImage, 'CData', bin_stack(:, :, current_z, current_channel));
    subplot(1, 2, 2);
    title(sprintf('Binarized - Z-plane %d of %d, Channel %d of %d', current_z, z_size, current_channel, channel_count));
end

%% Binarization function
function bw = binarize_adapt(I, sensitivity)
    % Normalize the image if necessary
    if isa(I, 'uint16') || isa(I, 'uint32')
        I = mat2gray(I);  % Convert to double in [0,1]
    end

    % Apply adaptive thresholding with the specified sensitivity
    bw = imbinarize(I, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', sensitivity);
end
