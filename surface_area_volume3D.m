
% would you like to load yuour data?
loadData = 1; % 1 for yes, 0 for no

addpath('data')
if loadData == 1
    clear all;
    pericyteImg = tiffreadVolume('data/J7568-EPI-BV1-5nm-PC.tif'); % pericyte
end

cc = bwconncomp(pericyteImg);
rp=regionprops(cc);

% dimensions
dz = 70; %nm
dx = 5; %nm
dy = 5; %nm

pericytePerimeter3d = bwperim(pericyteImg, 4);


% Step 1: Generate the isosurface
tic
fv = isosurface(pericytePerimeter3d, 0.5);
fprintf('isosurface: %g \n', toc./60);

% Step 2: Scale the vertices of the isosurface according to the voxel dimensions
tic
fv.vertices(:, 1) = fv.vertices(:, 1) * dx;  % Scale x-coordinates
fv.vertices(:, 2) = fv.vertices(:, 2) * dy;  % Scale y-coordinates
fv.vertices(:, 3) = fv.vertices(:, 3) * dz;  % Scale z-coordinates
fprintf('vertices data: %g \n', toc./60);

surfaceArea = 0;  % Initialize surface area

numberOfFaces = size(fv.faces, 1);

for i = 1:numberOfFaces
    tic
    % Get the vertex indices for the i-th triangle
    v1 = fv.vertices(fv.faces(i, 1), :);  % First vertex
    v2 = fv.vertices(fv.faces(i, 2), :);  % Second vertex
    v3 = fv.vertices(fv.faces(i, 3), :);  % Third vertex
    
    % Calculate the area of the triangle
    triangleArea = 0.5 * norm(cross(v2 - v1, v3 - v1));
    
    % Add this triangle's area to the total surface area
    surfaceArea = surfaceArea + triangleArea;
    fprintf('Faces time left: %g \n',(numberOfFaces-i) .* toc./60)
end

fprintf('Surface Area: %.2f nm^2\n', surfaceArea);

surfaceArea_um2 = surfaceArea * 1e-6;
fprintf('Surface Area: %.6f µm^2\n', surfaceArea_um2);