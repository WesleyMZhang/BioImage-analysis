% Goal is to get one clean stack to calculate the bacteria occupied voxels
%User inputs: 
%Section 1: Folder with raw data (parentfolder), Stack to process (ti,b,p,r,t) and 
%Section 2: Plane to visualize (plane2vis)
clear; close all;
parentfolder = 'C:\Users\zhang\Dropbox (GaTech)\curtis_lab\Curtis_Lab Data from Students\Hemaa Selvakumar\Data\Ti_4\8\4_1\T00001'; % Add the path to the folder with raw data

% Populate the following variables to access the corresponding data
ti=2; %Enter 1 for Ti0, 2 for Ti4, 3 for Ti8, 4 for Ti12
b=3; %Enter 1 for 10^6 CFU/mL, 2 for 10^7 CFU/mL, 3 for 10^8 CFU/mL 
p=0;  %Enter 0 for no phage control, n = 1 or 2 or .. or 9 for 10^n PFU/mL 
r=3;  %Enter r = 1 or 2 or 3 for the corresponsding replicate
t=1;  % Enter t = 1 or 2 or .. or 15 for the corresponding timepoint

Ti=[0,4,8,12];
Bac = [6,7,8];

%Access the relevant folder containing image data
stackdir = parentfolder;

MyFolderInfo = dir(stackdir);
MyFolderInfo = remove_first2rows(MyFolderInfo);
str = MyFolderInfo(end).name;
str = convertCharsToStrings(str);

% Identify relevant parameters
num_channels = str2double(extractBetween(str,"C0","Z"));
num_planes = str2double(extractBetween(str,"Z0",".tif"));

%Create a matrix to store the binarized image data
pixel_area = zeros(num_planes, 3);
stack_g = zeros(1024,1024,num_planes);

%Loop for binarizing the image stack with multiple channel data
for n = 1:num_planes
    %Access the SYTO9 channel images
    image_file_g = strcat(stackdir,'\',MyFolderInfo(n).name);
    I_g = imread(image_file_g);
    stack_g(:,:,n) = I_g;
    %Binarize the SYTO9 channel images
    I_bn_g = binarize_adapt(I_g);
    %Access the Propidium Iodide channel images
    image_file_r = strcat(stackdir,'\',MyFolderInfo(num_planes+n).name);
    I_r = imread(image_file_r);
    %Binarize the Propidium Iodide channel images
    I_bn_r = binarize_adapt(I_r);
    %Compute the number of non-zero pixels in each plane after binarization
    pixel_area(n,1)= sum(I_bn_g,[1 2]);
    pixel_area(n,2)= sum(I_bn_r,[1 2]);
    %Find the overlaping pixels in SYTO9 and Propidium iodide channels
    I_bn_o = I_bn_g.*I_bn_r;
    %Compute the number of overlaping pixels in each plane
    pixel_area(n,3)= sum(I_bn_o,[1 2]); 
end

%Estimate the plane 1-3 microns below the sample surface 
g_area = squeeze(pixel_area(:,1));
if b==1 
    plane_low=find_plane_low(g_area,1);
else if b==2
        plane_low=find_plane_low(g_area,1);
else if b==3
        if ti==2
            plane_low=find_plane_low(g_area,1);
        else
            plane_low=find_plane_low(g_area,2);
        end
    end
    end
end
%Set the total stack to be 20 microns - chosen based on overall data as
%highest common number of planes that can be obtained for this dataset given 
%the drift of automated micrscope stages/objectives 
stacksize = 19;
plane_high = plane_low-stacksize;
if plane_high<=0
    plane_high = 1;
end
%Obtain the area data for this clean image stack of uniform size
g = squeeze(pixel_area(plane_high:plane_low,1));
y = squeeze(pixel_area(plane_high:plane_low,3));
%Subtract the overlapping area from the SYTO9 channel area to obtain the
%live area. Sum over the area and multiply by 1 Z pizel step size to get
%the binarized pixel volume
vol = sum(g)-sum(y);
conf_vol=vol;
%One bacteria occupies approximately 80 voxels in this setup
conf_bacnum = conf_vol/80;

%% Intermediate figure generation for demonstration purposes
%Gather relevant data
%Plane2vis is user input','b=10^'p
plane2vis = 12; % Plane # from the lowest plane of the final stack
plane_fig = plane_low-plane2vis;
image_plane_fig_g = strcat(stackdir,'\',MyFolderInfo(plane_fig).name);
image_plane_fig_r = strcat(stackdir,'\',MyFolderInfo(num_planes+plane_fig).name);
%Generate the images
figure;
tile = tiledlayout(4,2);
tile.TileSpacing = 'compact';
if p==0
    p2print = int2str(0);
else
    p2print = strcat('10^',int2str(p));
end
title(tile,{'Intermediate outputs from image processing',...
    strcat('Ti = ', int2str(Ti(ti)),'hr, ',' b = 10^',int2str(Bac(b)),'CFU/mL, ',...
    ' p = ',p2print,'PFU/mL, ',' Replicate = ', int2str(r),', ',...
    ' Timepoint = ',sprintf('%02d ',t))});
%Tile 1 in the figure
%Area occupied by bacteria in the SYTO9 channel vs plane # of image stack
nexttile
plot(1:num_planes, squeeze(pixel_area(1:num_planes,1)));
hold on
xline(plane_high,'--r', {'Final stack','Highest plane'}); 
xline(plane_low,'--r', {'Final stack','Lowest plane'}); 
hold off
xlabel('Plane number')
ylabel('SYTO9 pixel area');
title('Raw data stack');
hold off

%Tile 2 in the figure
%Optimized stack of uniform size - Area occupied by bacteria in the SYTO9
%channel vs plane # of image stack
nexttile
plot(plane_high:plane_low, squeeze(pixel_area(plane_high:plane_low,1)));
hold on
xline(plane_fig,'--b', {'','Visualized','below'});
hold off
xlabel('Plane number')
ylabel('SYTO9 pixel area');
title({'Zoomed-in Final Stack','of Uniform Size'});
hold off

%Tile 3 in the figure
%Intensity-scaled SYTO9 raw image at plane that is to be visualized
nexttile
I_fig = imread(image_plane_fig_g);
%C = [min(I_fig, [],'all') max(I_fig,[],'all')];
imagesc(imadjust(I_fig));
hold on
axis image
xlabel('Pixel number')
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
ylabel('Pixel number');
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
title({'Intensity-scaled SYTO9 raw image',sprintf('Plane %d',plane_fig)});
hold off

%Tile 4 in the figure
%Binarized SYTO9 image at plane that is to be visualized
nexttile
I_fig_bn = binarize_adapt(I_fig);
imshow(mat2gray(I_fig_bn,[0,1]));
%h=I_fig_bn;im=h.CData;im=imresize(im,10);
hold on
axis on
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
xlabel('Pixel number')
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
ylabel('Pixel number');
title({'Binarized SYTO9 image',sprintf('Plane %d',plane_fig)});
hold off

%Tile 5 in the figure
%Intensity-scaled PI raw image at plane that is to be visualized
nexttile
I_fig_r = imread(image_plane_fig_r);
imagesc(imadjust(I_fig_r));
hold on
axis image
xlabel('Pixel number')
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
ylabel('Pixel number');
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
title({'Intensity-scaled PI raw image',sprintf('Plane %d',plane_fig)});
hold off

%Tile 6 in the figure
%Binarized PI image at plane that is to be visualized
nexttile
I_fig_bn_r = binarize_adapt(I_fig_r);
imshow(mat2gray(I_fig_bn_r,[0,1]));%.*65535);
hold on
axis on
xlabel('Pixel number')
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
ylabel('Pixel number');
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
title({'Binarized PI image',sprintf('Plane %d',plane_fig)});
hold off

%Tile 7 in the figure
%Estimated pixels occupied by live bacteria on the plane that is to be visualized
nexttile
I_fig_bn_o = I_fig_bn.*I_fig_bn_r;
I_fig_bn_l = I_fig_bn-I_fig_bn_o;
imshow(mat2gray(I_fig_bn_l,[0,1]));%.*65535);
hold on
axis on
xlabel('Pixel number')
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
ylabel('Pixel number');
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
title({'Area occupied by live bacteria',sprintf('Plane %d',plane_fig)});
hold off

%Tile 8 in the figure
%Estimated pixels occupied by dead bacteria on the plane that is to be visualized
nexttile
I_fig_bn_r = binarize_adapt(I_fig_r);
imshow(mat2gray(I_fig_bn_r,[0,1]));%.*65535);
hold on
axis on
xlabel('Pixel number')
xticks([600,1200,1800]);
xticklabels({'600','1200','1800'});
ylabel('Pixel number');
yticks([600,1200,1800]);
yticklabels({'600','1200','1800'});
title({'Area occupied by dead image',sprintf('Plane %d',plane_fig)});
hold off
%% 
%Gather relevant data
%Plane2vis is user input','b=10^'p
% plane2vis = 3; % Plane # from the lowest plane of the final stack
% plane_fig = plane_low-plane2vis;
% image_plane_fig_g = strcat(stackdir,'\',MyFolderInfo(plane_fig).name);
% image_plane_fig_r = strcat(stackdir,'\',MyFolderInfo(num_planes+plane_fig).name);
% I_fig = imread(image_plane_fig_g);
% I_fig_bn = binarize_adapt(I_fig);
% I_fig_r = imread(image_plane_fig_r);
% I_fig_bn_r = binarize_adapt(I_fig_r);
% I_fig_bn_o = I_fig_bn.*I_fig_bn_r;
% I_fig_bn_l = I_fig_bn-I_fig_bn_o;
figure;
imshow(mat2gray(I_fig_bn_l,[0,1]));
cc = bwconncomp(I_fig_bn_l);

%3D Green Channel stack;
for n=1:num_planes
    I_temp = stack_g;
    %C = [min(I_fig, [],'all') max(I_fig,[],'all')];
    imagesc(imadjust(I_fig));
end
