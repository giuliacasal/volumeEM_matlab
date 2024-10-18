%split surfaces


clear all;
filename = 'j7784_bv2_10nm_ECnucleus';
addpath('data')

pericyteImg = tiffreadVolume(sprintf('data/%s.tif',filename)); % pericyte

cc = bwconncomp(pericyteImg);

n = cc.NumObjects;


for i = 1 : n

    N = zeros(size(pericyteImg));

    I = cc.PixelIdxList{i};

    N(I) = 1;

    output_file1 = sprintf('data/%s_Nuc%i.tif',filename,i);
    for k = 1:size(N, 3)
        imwrite(N(:, :, k), output_file1, 'WriteMode', 'append',  'Compression','none');
    end



end

