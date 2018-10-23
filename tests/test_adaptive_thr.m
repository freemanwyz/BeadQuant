%% read data
p = getInitParamBead();
p.crop = 1;
% p.crop_rgx = 1:500;
% p.crop_rgy = 1:2560;

p.crop_rgx = 1:3328;
p.crop_rgy = 1:500;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin))/255;
datbw = dat;

[Nx,Ny] = size(dat);
rx1 = [1 1];
rx2 = [500 3328];
rx2a = [200 3328];
ry1 = [1 1];
ry2 = [2560 500];
ry2a = [2400 150];
Nblk = 2;

% adaptive thresholding ---
for ii=1:Nblk
    crop_rgx = rx1(ii):rx2(ii);
    crop_rgy = ry1(ii):ry2(ii);
    dat0 = dat(crop_rgx,crop_rgy);
    IM = dat0;
    ws = 100;
    C = 0;
    tm = 1;
    if tm==0
        mIM=imfilter(IM,fspecial('average',ws),'replicate');
    else
        mIM=medfilt2(IM,[ws ws]);
    end
    sIM=mIM-IM+0.0005;
    bw=im2bw(sIM,0);
    bw=imcomplement(bw);
    % figure;imshow(bw);
    %     K5 = cat(3,bw,sqrt(dat0),dat0*0);
    %     figure;imshow(K5);
    % figure;imshow(dat0);
    %     K6 = cat(3,dat0>0.05,sqrt(dat0),dat0*0);
    % figure;imshow(K6);
    
    crop_rgxa = rx1(ii):rx2a(ii);
    crop_rgya = ry1(ii):ry2a(ii);
    datbw(crop_rgxa,crop_rgya) = bw(crop_rgxa,crop_rgya);
end

K5 = cat(3,datbw,dat,dat*0);
% K5 = cat(3,datbw,sqrt(dat),dat*0);
figure;imshow(K5);

imwrite(datbw,'myAdaMask1.tif');



