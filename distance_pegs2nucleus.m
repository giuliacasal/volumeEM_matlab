% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

dataName = 'J7785a_BV1';

addpath('data')
if loadData == 1
    clearvars -except dataName;
    PegImg = tiffreadVolume(sprintf('data/%s_10nm_PCpegs.tif',dataName)); % pericyte
    NucleusImg = tiffreadVolume(sprintf('data/%s_10nm_ECnuc.tif',dataName)); % pericyte
end

mkdir(sprintf('%s',dataName))

% Binary matrix representing the perimeter of the pericyte in 3D
% (generated using bwperim)
% pericytePerimeter3d = bwperim(pericyteImg, 4);

% Voxel dimensions in nanometers
dx = 10; % x-dimension (nm)
dy = 10; % y-dimension (nm)
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
cc_nuc = bwconncomp(NucleusImg);

% number of nuclei and pegs
n_peg = cc_peg.NumObjects;
n_nuc = cc_nuc.NumObjects;

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

for i = 1 : length(z)

      nuc(:,:) =  NucleusImg(:,:,i);
       
        b_nuc = bwboundaries(nuc);

        if isempty(b_nuc)==0
        x_nuc = [x_nuc x(b_nuc{1}(:,1))];
        y_nuc = [y_nuc y(b_nuc{1}(:,2))];
        z_nuc = [z_nuc z(i).*ones(size(y(b_nuc{1}(:,2))))];
        end

end

% figure(3); cla; hold on;
% plot3(x_peg,y_peg,z_peg,'o')
% plot3(x_nuc,y_nuc,z_nuc,'.')


d_peg = zeros(n_peg,n_nuc);

% store the min distance co-ordinate of the nucleus
x_nuc_min = [];
y_nuc_min = [];
z_nuc_min = [];

 % going through all pegs
for j = 1 : n_peg

    % finding distance between peg centre and all nucleus pixels
    dis = sqrt( ( x_peg(j) - x_nuc ).^2 + ( y_peg(j) - y_nuc ).^2 + ( z_peg(j) - z_nuc) .^2 ) ; 

    % finding the minimum of those distances
    [d_peg(j,1), min_idx] = min(dis);

    % storing the x, y, z coordinates of the nucleus with the shortest distance
    x_nuc_min(j,1) = x_nuc(min_idx);
    y_nuc_min(j,1) = y_nuc(min_idx);
    z_nuc_min(j,1) = z_nuc(min_idx);

end

T = array2table([x_peg, y_peg, z_peg, d_peg]);
T.Properties.VariableNames = {'x','y','z','d'};
writetable(T,sprintf('%s/pegPositions.csv',dataName))
%csvwrite(sprintf('%s/pegPositions.csv'),[x_peg, y_peg, z_peg, d_peg])
% n_peg, x_peg, y_peg, z_peg, d_peg,


% figure(9); cla; hold on;
% plot3(x_peg,y_peg,z_peg,'ko')


figure(8); cla; hold on; axis equal; set(gca,'box','on')
plot3(y_peg,x_peg,z_peg,'ko')
for i = 1 : length(z)

    peg(:,:) =  PegImg(:,:,i);

     nuc(:,:) =  NucleusImg(:,:,i);

    b_peg = bwboundaries(peg);
    b_nuc = bwboundaries(nuc);

    for j = 1 : length(b_peg)
        plot3( x(b_peg{j}(:,1)), y(b_peg{j}(:,2)),z(i).*ones(size(b_peg{j}(:,2))), 'color', [228,26,28]./255)
    end

    for j = 1 : length(b_nuc)
        plot3( x(b_nuc{j}(:,1)), y(b_nuc{j}(:,2)),z(i).*ones(size(b_nuc{j}(:,2))), 'color', [179,205,227]./255  )
    end

    pause(0); drawnow
end



view(45,15)
exportgraphics(gcf,sprintf('%s/3Dplot.png',dataName),'resolution',1000)