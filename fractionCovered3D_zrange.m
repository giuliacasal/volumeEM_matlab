
% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

addpath('data')
if loadData == 1
    clear all;
    endotheliumImg = tiffreadVolume('data/J7568-EPI-BV1-5nm-EC.tif'); % endothelium
    pericyteImg = tiffreadVolume('data/J7568-EPI-BV1-5nm-PC.tif'); % pericyte
end

% if plotting results plt = 1, if not plt =0
plt = 1;

% threshold distance
contactDistance = 30; % nm
coveredDistance = 200; % nm
depthDistance = 150; % nm

% dimensions
dz = 70; %nm
dx = 5; %nm
dy = 5; %nm

% size of images
s  = size(endotheliumImg); xLength = s(2); yLength = s(1); zLength = s(3);

%if plot looks ridicululous i.e; x and y axis values are reversed - uncomment this
% s  = size(E); xLength = s(1); yLength = s(2); zLength = s(3);

% find perimeter
endotheliumPerimeter3d = bwperim(endotheliumImg,4);
pericytePerimeter3d = bwperim(pericyteImg,4);

% sample perimeter slice
samplePerimeterSlice = endotheliumPerimeter3d(:,:,round(zLength/2));
estimatedPerimeterPointsPerSlice = floor(nnz(samplePerimeterSlice) * 1.5); 

% X and Y coordinates of all image pixels
[X,Y] = meshgrid((1:xLength),(1:yLength));
X = X.*dx;
Y = Y.*dy;

zSliceStart = 1;

zSlices = (1 : zLength) .* dz;

% Define the range of z-slices to analyze
z_min = 204; % Example minimum z-slice
z_max = 354; % Example maximum z-slice

% Ensure the range is within the available slices
z_min = max(z_min, zSliceStart);
z_max = min(z_max, zLength);

n = floor(depthDistance/dz);

fCon = 0;
fCov = 0;
fTot = 0;

% Loop through the specified z-slices range
for i = z_min : z_max
    tic

    neighborZSlices = max(z_min,i-n) : min(z_max,i+n);
    endotheliumPerimeter2d(:,:) = endotheliumPerimeter3d(:,:,i);
    endotheliumPerimeterIndices = find(endotheliumPerimeter2d==1);

    contactPerimIndices = [];
    coveragePerimIndices = [];

    for j = neighborZSlices      
        pericytePerimeter2d(:,:) = pericytePerimeter3d(:,:,j);  
        pericytePerimeterIndices = find(pericytePerimeter2d==1);

        distanceMatrix = sqrt((X(endotheliumPerimeterIndices)-X(pericytePerimeterIndices)').^2 ...
            + (Y(endotheliumPerimeterIndices)-Y(pericytePerimeterIndices)').^2 ...
            + (zSlices(i)-zSlices(j)).^2);

        [contactPointIndices, ~] = find(distanceMatrix<contactDistance);
        [coveragePointIndices, ~] = find(distanceMatrix<coveredDistance);

        contactPerimIndices = [contactPerimIndices; contactPointIndices];
        coveragePerimIndices = [coveragePerimIndices; coveragePointIndices];
    end

    contactPerimIndices = unique(contactPerimIndices);
    xxCon = X(endotheliumPerimeterIndices(contactPerimIndices));
    yyCon = Y(endotheliumPerimeterIndices(contactPerimIndices));

    coveragePerimIndices = unique(coveragePerimIndices);
    xxCov = X(endotheliumPerimeterIndices(coveragePerimIndices));
    yyCov = Y(endotheliumPerimeterIndices(coveragePerimIndices));

    endotheliumPerimeter2d(:,:) = endotheliumPerimeter3d(:,:,i);
    pericytePerimeter2d(:,:) = pericytePerimeter3d(:,:,i);

    numConPixels = length(xxCon);
    numCovPixels = length(xxCov);
    numEndoPixels = length(endotheliumPerimeterIndices);

    fCon = fCon + numConPixels;
    fCov = fCov + numCovPixels;
    fTot = fTot + numEndoPixels;

    fprintf('Time left: %g min - frac covered: %g - frac connected: %g \n',(zLength-i) .* toc./60, fCov./fTot, fCon./fTot)

     if plt==1

        pericytePerimeterIndices = find(pericytePerimeter2d==1);
        figure(1); cla; hold on;
        plot(X(endotheliumPerimeterIndices),Y(endotheliumPerimeterIndices),'.','color',[179,205,227]./255)
        plot(X(pericytePerimeterIndices),Y(pericytePerimeterIndices),'.','color',[253,191,111]./255)
        plot(xxCov,yyCov,'b.')
        plot(xxCon,yyCon,'r.')
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

fracConnected = fCon./fTot;
fracCovered = fCov./fTot;

fprintf('Fraction connected: %g \n Fraction covered: %g ',fracConnected,fracCovered)

