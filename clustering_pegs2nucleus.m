%THIS ONE WORKS
% works for J7785a_BV1, J7784_BV1,_BV2,J7568-bv-2,J7568-EPI-BV1 
%does not work for J7571_BV1, J7568-epi-bv-2
% clc
clear
% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

dataName = 'J7568-bv3';

addpath('data')
if loadData == 1
    clearvars -except dataName;
    PegImg = tiffreadVolume(sprintf('data/%s-12nm-411_PCpeg.tif',dataName)); % pericyte
    %NucleusImg = tiffreadVolume(sprintf('data/%s_10nm_PC.tif',dataName)); % pericyte
end

mkdir(sprintf('%s',dataName))

% Binary matrix representing the perimeter of the pericyte in 3D
% (generated using bwperim)
% pericytePerimeter3d = bwperim(pericyteImg, 4);

% Voxel dimensions in nanometers
dx = 12; % x-dimension (nm)
dy = 12; % y-dimension (nm)
dz = 70; % z-dimension (nm)

% dims
s = size(PegImg);

% 3D mesh
x = (1:s(1)).*dx;
y = (1:s(2)).*dy;
z = (1:s(3)).*dz;
[Y,X,Z] = meshgrid(x,y,z);

% rici's connected components
cc_peg = bwconncomp(PegImg);
% cc_nuc = bwconncomp(NucleusImg);



% number of nuclei and pegs
n_peg = cc_peg.NumObjects;
% n_nuc = cc_nuc.NumObjects;
n_nuc = 1;
% finding the central coords of the pegs
x_peg = zeros(n_peg,1);
y_peg = zeros(n_peg,1);
z_peg = zeros(n_peg,1);

rp_peg=regionprops(cc_peg);

for i = 1:n_peg

    x_peg(i) = rp_peg(i).Centroid(1).*dx;
    y_peg(i) = rp_peg(i).Centroid(2).*dy;
    z_peg(i) = rp_peg(i).Centroid(3).*dz;

end

% boundary coord of nucleus
x_nuc = [];
y_nuc = [];
z_nuc = [];

% for i = 1 : length(z)
% 
%       nuc(:,:) =  NucleusImg(:,:,i);
%        
%         b_nuc = bwboundaries(nuc);
% 
%         if isempty(b_nuc)==0
%             x_nuc = [x_nuc, x(b_nuc{1}(:,1))];
%             y_nuc = [y_nuc, y(b_nuc{1}(:,2))];
%             z_nuc = [z_nuc, z(i).*ones(size(y(b_nuc{1}(:,2))))];
%         end
% 
% end

% figure(3); cla; hold on;
% plot3(x_peg,y_peg,z_peg,'o')
% plot3(x_nuc,y_nuc,z_nuc,'.')


d_peg = zeros(n_peg,n_nuc);
i_nuc = zeros(n_peg,n_nuc);

 % going through all pegs
% for j = 1 : n_peg
% 
%     % finding distance between peg centre and all nucleus pixels
%     [dis, idx] = min( sqrt( ( x_peg(j) - y_nuc ).^2 + ( y_peg(j) - x_nuc ).^2 + ( z_peg(j) - z_nuc) .^2 ) ); 
% 
%     d_peg(j) = dis;
%     i_nuc(j) = idx;
% 
% end

% Define clustering parameters
pegCoords = [x_peg, y_peg, z_peg];

epsilon = 800; % Maximum distance for neighbourhood
minPts = 3;    % Minimum pegs required for a cluster

clusterLabels = dbscan(pegCoords, epsilon, minPts);

uniqueClusters = unique(clusterLabels);
numClusters = sum(uniqueClusters > 0); % Ignore noise (-1)

fprintf('Number of clusters: %d \n',numClusters);

% Convert DBSCAN labels to a logical clustered variable
isClustered = clusterLabels > 0; % Ignore noise (-1)

% Compute total number of pegs
totalPegs = length(clusterLabels);

% Count the number of pegs that belong to a cluster (excluding noise, -1)
clusteredPegs = sum(clusterLabels > 0);

% Compute percentage of clustered pegs
percentageClustered = (clusteredPegs / totalPegs) * 100;

% Print the result
fprintf('Percentage of pegs in a cluster: %.2f%%\n', percentageClustered);

figure(8); cla; hold on; axis equal; set(gca,'box','on')

for i = 1 : length(z)

    peg(:,:) =  PegImg(:,:,i);

    b_peg = bwboundaries(peg);

    for j = 1 : length(b_peg)
        plot3( x(b_peg{j}(:,1)), y(b_peg{j}(:,2)),z(i).*ones(size(b_peg{j}(:,2))), 'color', [0,0,0]./255)
    end

    pause(0); drawnow
end

for i = 1:n_peg
    if isClustered(i)
        plot3(pegCoords(i,2), pegCoords(i,1), pegCoords(i,3), 'o', 'MarkerEdgeColor', [0, 0.5, 0.5], 'MarkerFaceColor', [0, 0.5, 0.5], 'MarkerSize', 6) % Green for clustered
    else
        plot3(pegCoords(i,2), pegCoords(i,1), pegCoords(i,3), 'o', 'MarkerEdgeColor', [0.4, 0.2, 0], 'MarkerFaceColor', [0.4, 0.2, 0], 'MarkerSize', 6) % Red for isolated
    end
end

view(45,15)
exportgraphics(gcf,sprintf('%s/3Dplot.png',dataName),'resolution',1000)

statsTable = table(numClusters, percentageClustered, ...
    'VariableNames', {'NumClusters', 'PercentageClustered'});

writetable(statsTable, sprintf('%s/clustering_summary.csv', dataName));
