%% tif stacks to RGB channel
fname = './dat/MAX_EGFP+synapsin+PSD 952.tif';
info = imfinfo(fname);
num_images = numel(info);

A = zeros(512,512,3);
for k = 1:num_images
    A(:,:,k) = im2double(imread(fname, k, 'Info', info));
end

fname_out = './dat/MAX_EGFP+synapsin+PSD 952_RGB.tif';
imwrite(A,fname_out);

%% three tif to one
fname = './dat/ctrl neuron-ctrl astro-Tuj1, Syn I.lsm - Ch3-T3 - C3 Z1 T1.tif';