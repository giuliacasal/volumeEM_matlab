%split surfaces


clear all;
filename = 'J7784-BV1-10nm-ECnuclei';
addpath('data')

pericyteImg = tiffreadVolume(sprintf('data/%s.tif',filename)); % pericyte

cc = bwconncomp(pericyteImg);

n = cc.NumObjects;

%in case smaller, reduce number
threshSize = 1e2;

ii = 0;
for i = 1 : n

    N = zeros(size(pericyteImg));

    if length(cc.PixelIdxList{i})>threshSize

        ii = ii + 1;

        I = cc.PixelIdxList{i};
    
        N(I) = 1;
    
        output_file1 = sprintf('data/%s_Nuc%i.tif',filename,ii);
        for k = 1:size(N, 3)
            imwrite(N(:, :, k), output_file1, 'WriteMode', 'append',  'Compression','none');
        end

        fprintf('Number of nuclei: %i \n', ii)
    end



end

