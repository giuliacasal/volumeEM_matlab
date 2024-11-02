
% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

addpath('data')
if loadData == 1
    clear all;
    endotheliumImg = tiffreadVolume('data/j7784_bv2_10nm_ECfilled.tif'); % endothelium
    pericyteImg = tiffreadVolume('data/j7784_bv2_10nm_PC.tif'); % pericyte
end

% if plotting results plt = 1, if not plt =0
plt = 1;

% threshold distance
contactDistance = 30; % nm
coveredDistance = 150; % nm
depthDistance = 150; % nm

% dimensions
dz = 70; %nm
dx = 10; %nm
dy = 10; %nm

% size of images
s  = size(endotheliumImg); xLength = s(2); yLength = s(1); zLength = s(3);

%if plot looks ridicululous i.e; x and y axis values are reversed - uncomment this
% s  = size(E); xLength = s(1); yLength = s(2); zLength = s(3);

% find edges
endotheliumPerimeter3d = bwperim(endotheliumImg,4);
pericytePerimeter3d = bwperim(pericyteImg,4);

% X and Y coordinates of all image pixels
[X,Y] = meshgrid((1:xLength),(1:yLength));
X = X.*dx;
Y = Y.*dy;

zSliceStart = 1;

zSlices = (1 : zLength) .* dz;

% Define the range of z-slices to analyze
zMin = 290; % Example minimum z-slice
zMax = 299; % Example maximum z-slice

% Ensure the range is within the available slices
zMin = max(zMin, zSliceStart);
zMax = min(zMax, zLength);

depthSliceRange = floor(depthDistance/dz);

totalConnected = 0;
totalCovered = 0;
totalEndothelium = 0;

% Loop through the specified z-slices range
for i = zMin : zMax
    tic

    neighborZSlices = max(zMin,i-depthSliceRange) : min(zMax,i+depthSliceRange);

    contactPerimeterIndices = [];
    coveragePerimeterIndices = [];

    endotheliumPerimeter2d(:,:) = endotheliumPerimeter3d(:,:,i);
    endotheliumPerimeterIndices = find(endotheliumPerimeter2d==1);

    for j = neighborZSlices      
        pericytePerimeter2d(:,:) = pericytePerimeter3d(:,:,j);  
        ip = find(pericytePerimeter2d==1);

        distanceMatrix = sqrt((X(endotheliumPerimeterIndices)-X(ip)').^2 ...
            + (Y(endotheliumPerimeterIndices)-Y(ip)').^2 ...
            + (zSlices(i)-zSlices(j)).^2);

        [contactPointIndices, ~] = find(distanceMatrix<contactDistance);
        [coveragePointIndices, ~] = find(distanceMatrix<coveredDistance);

        contactPerimeterIndices = [contactPerimeterIndices; contactPointIndices];
        coveragePerimeterIndices = [coveragePerimeterIndices; coveragePointIndices];
    end

    contactPerimeterIndices = unique(contactPerimeterIndices);
    xContactPoints = X(endotheliumPerimeterIndices(contactPerimeterIndices));
    yContactPoints = Y(endotheliumPerimeterIndices(contactPerimeterIndices));

    coveragePerimeterIndices = unique(coveragePerimeterIndices);
    xCoveragePoints = X(endotheliumPerimeterIndices(coveragePerimeterIndices));
    yCoveragePoints = Y(endotheliumPerimeterIndices(coveragePerimeterIndices));

    endotheliumPerimeter2d(:,:) = endotheliumPerimeter3d(:,:,i);
    pericytePerimeter2d(:,:) = pericytePerimeter3d(:,:,i);

    numContactPixels = length(xContactPoints);
    numCoveragePixels = length(xCoveragePoints);
    numEndotheliumPixels = length(endotheliumPerimeterIndices);

    totalConnected = totalConnected + numContactPixels;
    totalCovered = totalCovered + numCoveragePixels;
    totalEndothelium = totalEndothelium + numEndotheliumPixels;

    fprintf('Time left: %g min - frac covered: %g - frac connected: %g \n',(zLength-i) .* toc./60, totalCovered./totalEndothelium, totalConnected./totalEndothelium)

     if plt==1

        ip = find(pericytePerimeter2d==1);
        figure(1); cla; hold on;
        plot(X(endotheliumPerimeterIndices),Y(endotheliumPerimeterIndices),'.','color',[179,205,227]./255)
        plot(X(ip),Y(ip),'.','color',[253,191,111]./255)
        plot(xCoveragePoints,yCoveragePoints,'b.')
        plot(xContactPoints,yContactPoints,'r.')
        axis([0 xLength.*dx 0 yLength.*dy]); axis equal
        set(gca,'box','on','linewidth',2)
        pause(0); drawnow


%         figure(2);
%         plot3(X(ie),Y(ie),z(i).*ones(size(ie)),'.','color',[179,205,227]./255)
%         plot3(X(ip),Y(ip),z(i).*ones(size(ip)),'.','color',[253,191,111]./255)
%         plot3(xx,yy,z(i).*ones(size(xx)),'r.')
%         bp = nan;

    end

end

fracConnected = totalConnected./totalEndothelium;
fracCovered = totalCovered./totalEndothelium;

fprintf('Fraction connected: %g \n Fraction covered: %g ',fracConnected,fracCovered)

