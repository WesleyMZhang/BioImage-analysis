function plane_low = find_plane_low(g_area,option)
%Function to identify the plane within 1-3 microns of the sample surface.
%The option is identified based on the imaging conditions

num_planes = size(g_area,1);
if option ==1
    %Flip the plane indices of SYTO9 binarized areas in the stack
    g_area = flipud(g_area);
    %Find the flipped index whose area is >1.5 bacteria atleast
    x = find(g_area>120,1);
    %Subtract the flipped index from total number of planes to obatin the
    %plane within 1-3 microns of the sample surface
    x = num_planes-x;
    plane_low=x;
else
    %Find the plane index of image with <1.5 bacteria 
    x = find(g_area<120);
    %If the index is above 15, set the minimum of x as the plane within 1-3
    %microns of the sample surface
    x = x(x>15);
    l = min(x);
    %If the sample surface is the last image in the stack or not captured
    if isempty(l)
        l =num_planes;
    end
    plane_low=l;
end