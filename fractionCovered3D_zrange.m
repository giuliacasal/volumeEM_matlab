
% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

addpath('data')
if loadData == 1
    clear all;
    E = tiffreadVolume('data/J7784-BV1-10nm-EC.tif'); % endothelium
    P = tiffreadVolume('data/J7784-BV1-10nm-PC.tif'); % pericyte
end

% if plotting results plt = 1, if not plt =0
plt = 1;

% threshold distance
contactDistance = 30; % nm
coveredDistance = 150; % nm

% dimensions
dz = 70; %nm
dx = 10; %nm
dy = 10; %nm

% size of images
s  = size(E); xL = s(2); yL = s(1); zL = s(3);

%if plot looks ridicululous - uncomment this
% s  = size(E); xL = s(1); yL = s(2); zL = s(3);

% X and Y coordinates of all image pixels
[X,Y] = meshgrid((1:xL),(1:yL));
X = X.*dx;
Y = Y.*dy;

z0 = 1;

z = (1 : zL) .* dz;

% Define the range of z-slices to analyze
z_min = 140; % Example minimum z-slice
z_max = 290; % Example maximum z-slice

% Ensure the range is within the available slices
z_min = max(z_min, z0);
z_max = min(z_max, zL);

% find edges
ee = bwperim(E,4);
pe = bwperim(P,4);

n = floor(coveredDistance/dz);

fCon = 0;
fCov = 0;
fTot = 0;

% Loop through the specified z-slices range
for i = z_min : z_max
    tic

    J = max(z_min,i-n) : min(z_max,i+n);

    uCon = [];
    uCov = [];

    e(:,:) = ee(:,:,i);
    ie = find(e==1);

    for j = J      
        p(:,:) = pe(:,:,j);  
        ip = find(p==1);

        d = sqrt( (X(ie)-X(ip)').^2 + (Y(ie)-Y(ip)').^2 + (z(i)-z(j)).^2);
        [rowCon, ~] = find(d<contactDistance);
        [rowCov, ~] = find(d<coveredDistance);

        uCon = [uCon; rowCon];
        uCov = [uCov; rowCov];
    end

    uCon = unique(uCon);
    xxCon = X(ie(uCon));
    yyCon = Y(ie(uCon));

    uCov = unique(uCov);
    xxCov = X(ie(uCov));
    yyCov = Y(ie(uCov));

    e(:,:) = ee(:,:,i);
    p(:,:) = pe(:,:,i);

    numConPixels = length(xxCon);
    numCovPixels = length(xxCov);
    numEndoPixels = length(ie);

    fCon = fCon + numConPixels;
    fCov = fCov + numCovPixels;
    fTot = fTot + numEndoPixels;

    fprintf('Time left: %g min - frac covered: %g - frac connected: %g \n',(zL-i) .* toc./60, fCov./fTot, fCon./fTot)

     if plt==1

        ip = find(p==1);
        figure(1); cla; hold on;
        plot(X(ie),Y(ie),'.','color',[179,205,227]./255)
        plot(X(ip),Y(ip),'.','color',[253,191,111]./255)
        plot(xxCov,yyCov,'b.')
        plot(xxCon,yyCon,'r.')
        axis([0 xL.*dx 0 yL.*dy]); axis equal
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

