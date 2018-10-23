ftop = 'D:\neuro_WORK\coculture_ucdavis_WORK\GECI_test\';
runName = 'geci_nn-4001-sv4-run-0908-mean\';

load([ftop runName,'res_cor.mat']);
load([ftop runName,'res_temp.mat']);

% idNow = roi_map(338,177);
idNow = 12;

map0 = (roi_map==idNow)*1;
B = bwboundaries(map0,8,'noholes');
B1 = B{1};
B1x = B1(:,2:-1:1) - 0.5;

% map1 = map0;
% B1idx = sub2ind(size(map1),B1(:,1),B1(:,2));
% map1(B1idx) = 1:size(B1,1);

save('12.txt','B1x','-ASCII');

roi1 = (roi_map>0)*0.3;
roi1(roi_map==23) = 1;
imshow(roi1)