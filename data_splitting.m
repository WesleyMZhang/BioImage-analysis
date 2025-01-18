% Add Bio-Formats library to the MATLAB path
addpath('C:\Users\mzhang658\OneDrive - Georgia Institute of Technology\Desktop\Phage-biofilm codes_HS\bfmatlab'); % Replace with the actual path to bfmatlab

% Specify the ND2 file to read
nd2File = 'E:\0929_Ti_12\Ti_12.nd2'; % Replace with your ND2 file name

% Initialize the Bio-Formats reader
reader = bfGetReader(nd2File);

% Get the total number of series in the ND2 file
seriesCount = reader.getSeriesCount();

% Define the series ranges
ranges = {1:30, 31:60, 61:90};
outputFiles = {'series1_30.ome.tiff', 'series31_60.ome.tiff', 'series61_90.ome.tiff'};

% Loop over each range to extract and save the corresponding series
for idx = 1:length(ranges)
    seriesRange = ranges{idx};
    outputFile = outputFiles{idx};
    
    % Initialize an empty cell array to store images
    images = {};
    
    for s = seriesRange
        reader.setSeries(s - 1); % Series indexing starts from 0 in Bio-Formats
        numImages = reader.getImageCount();
        
        for i = 1:numImages
            % Read the image plane
            plane = bfGetPlane(reader, i);
            images{end+1} = plane;
        end
    end
    
    % Convert the cell array to a 3D array if all images are the same size
    if ~isempty(images)
        imageSize = size(images{1});
        numPlanes = length(images);
        imageStack = zeros([imageSize, numPlanes], class(images{1}));
        
        for k = 1:numPlanes
            imageStack(:,:,k) = images{k};
        end
        
        % Save the images to an OME-TIFF file
        bfsave(imageStack, outputFile);
        disp(['Saved ', outputFile]);
    else
        disp(['No images found for series range ', num2str(seriesRange)]);
    end
end

% Close the reader
reader.close();
