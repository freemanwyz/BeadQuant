%% Bead signal from grace
% Compatible with old version Matlab
p = getInitParamBead();

fin = [p.tp filesep 'dump' filesep '20150627_mask' filesep ...
    'res_bead_mask_1_thr_0p08_score_0p1_avoid_0_0p4_rad_1_median.mat'];
tmp = load(fin);
resBead = tmp.resBead;

fin = [p.tp filesep 'dump' filesep 'beadAvg16.tif'];
datAvg = double(imread(fin))/65535;

fin = [p.tp filesep 'dat' filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128' filesep ...
    '061815 Background 488 3.5%750_561 1% 800.lsm - Ch1-T1 - C1 Z1 T1.tif'];
datBg = double(imread(fin))/65535;

fin = [p.tp filesep 'dat' filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128' filesep ...
    '061815 Background 488 3.5%750_561 1% 800.lsm - Ch2-T2 - C2 Z1 T1.tif'];
sigBg = double(imread(fin))/65535;

[Nx,Ny] = size(sigBg);
nFrame = 32;

%% alignment
d0 = datAvg;
d1 = datBg;
rt0 = mean(d1(:))/mean(d0(:));
d1 = d1/rt0;

par.transform = 'euclidean'; 
par.levels = 2;
par.iterations = 20; %iterations per level
d1Warp = iat_LucasKanade(d1, d0, par); 
[d1_rigid, d1_rigid_support] = iat_inverse_warping(d1, d1Warp, par.transform, 1:Ny, 1:Nx);
[~, grayerrorECC] = iat_error2gray(d0,d1_rigid,d1_rigid_support);
figure;imshow(grayerrorECC);

figure;imshowpair(d1_rigid, d0,'Scaling','joint');title('Alignment error');

% [optimizer, metric]  = imregconfig('monomodal');
% % d1t_trans = imregtform(d1,d0,'translation',optimizer, metric);
% % d1_trans = imwarp(d1,d1t_trans,'OutputView',imref2d(size(d0)));
% 
% d1t_rigid = imregtform(d1,d0,'rigid',optimizer, metric);
% d1_rigid = imwarp(d1,d1t_rigid,'OutputView',imref2d(size(d0)));
% 
% % d1Registeredx = imwarp(d1,d1tform);
% % d1Registered = imregister(d1,d0,'rigid',optimizer, metric);
% 
% figure;imshowpair(d0, d1,'Scaling','joint');
% % figure;imshowpair(d0, d1_trans,'Scaling','joint');
% figure;imshowpair(d0, d1_rigid,'Scaling','joint');

%% ch2, signal
[sigBg_rigid, sigBg_rigid_support] = iat_inverse_warping(sigBg, d1Warp, par.transform, 1:Ny, 1:Nx);
% sigBg_rigid = imwarp(sigBg,d1t_rigid,'OutputView',imref2d(size(sigBg)));

fin = [p.tp filesep 'dat' filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128' filesep ...
    '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128.lsm - Ch2-T2 - C2 Z1 T'];
sig = zeros(Nx,Ny,nFrame+1);
sig(:,:,1) = sigBg_rigid;

for ii=1:nFrame
    fin0 = [fin num2str(ii) '.tif'];
    sig(:,:,ii+1) = double(imread(fin0))/65535;
end

%% find curves for each bead
nBeads = length(resBead);
myCurves = zeros(nBeads,nFrame+1);
SE7 = strel('square',7);
SE5 = strel('square',5);
SE3 = strel('square',3);

for ii=1:nBeads
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    pix0 = resBead{ii};
    n0 = length(pix0);
    pix0x = pix0(:,1);
    pix0y = pix0(:,2);
    rga = (min(pix0x)-3):(max(pix0x)+3);
    rgb = (min(pix0y)-3):(max(pix0y)+3);
    pix0x = pix0x - min(pix0x) + 1 + 3;
    pix0y = pix0y - min(pix0y) + 1 + 3;
    xx = zeros(length(rga),length(rgb));
    pix0 = sub2ind(size(xx),pix0x,pix0y);
    xx(pix0) = 1;
    if n0>125
        xx = imerode(xx,SE7);
%         xx = imerode(xx,SE5);
    end
    if n0>75 && n0<=125
%         xx = imerode(xx,SE3);
        xx = imerode(xx,SE5);
    end
    if n0<=75
        xx = imerode(xx,SE3);
    end
    pix1 = find(xx>0);
    sig0 = reshape(sig(rga,rgb,:),[],nFrame+1);
    sig0 = sig0(pix1,:);
    myCurves(ii,:) = mean(sig0,1);
%     for jj=1:(nFrame+1)
%         sig0 = sig(rga,rgb,ii);
%         c0(jj) = mean(sig0(pix1));
%     end
%     myCurves(ii,:) = c0;
end

%% plot F - F0
myCurvesDf = myCurves*0;
myCurvesDff = myCurves*0;
x00 = myCurves(:,1);
x00(x00<0.02) = 0.02;
for ii=1:nFrame
    myCurvesDf(:,ii+1) = myCurves(:,ii+1) - myCurves(:,1);
    myCurvesDff(:,ii+1) = (myCurves(:,ii+1) - myCurves(:,1))./x00;
end
figure;plot(myCurvesDf(:,1:33)');xlabel('Time points');ylabel('\DeltaF');
figure;plot(myCurvesDff(:,1:33)');xlabel('Time points');ylabel('\DeltaF/F');

%% analysis
xC = max(myCurvesDf,[],2);
xI = find(xC>0.7 & myCurvesDf(:,6)>0.3);

neibVec = [0, -1, 1, -Nx, Nx];
t0 = zeros(Nx,Ny);
for ii=1:length(xI)
    idx = resBead{xI(ii)};
    idx = sub2ind(size(t0),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);    
    t0(idxc) = 1;
end
K0 = cat(3,t0,datAvg,datAvg*0);
figure;imshow(K0);










