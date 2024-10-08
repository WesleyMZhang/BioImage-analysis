function binary_image = binarize_adapt(image)
%Binarize a 16-bit image using adaptive thresholding
%Convert the image into 16-bit unsigned integer data
image = uint16(image);
%Apply adaptive thresholding with a sensitivity of 0.4
T = adaptthresh(image,0.4,'Statistic','mean');
%Binarize the image using the adaptive threshold filter
BW = imbinarize(image,T);
%Perform morphological operations tuned to recognize bacteria (longitudinal
%ans cross-sections)
se = strel('disk', 2);
BW = imclose(BW, se);
BW = bwareaopen(BW, 12);
binary_image = BW;

