%% Bead signal from grace
p = getInitParamBead();
p.crop = 0;
% p.crop_rgx = 1501:1700;
% p.crop_rgy = 1301:1600;

% p.crop_rgx = 651:800;
% p.crop_rgy = 1051:1150;

% p.crop_rgx = 3178:3328;
% p.crop_rgy = 1:1500;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin));
% myRg = 4:13;
myRg = 4:0.5:13;
myRef = beadTemplate(myRg);

% To isolated parts ----
if p.crop
    dat0 = dat(p.crop_rgx,p.crop_rgy)/255;
else
    dat0 = dat/255;
end
[Nx,Ny] = size(dat0);
thr0 = 0.08;

% fin = [p.tp filesep 'dump' filesep 'myAdaMask1.tif'];
% myMask = double(imread(fin))/255;
% dat0ada = max(dat0,myMask);
% dat0ada(:,2497:end) = 0;
% dat0Mask = dat0ada>thr0;

dat0Mask = dat0>thr0;

dat0Mask = bwareaopen(dat0Mask, 10, 4);
dat1 = dat0.*dat0Mask;
dat1CC = bwconncomp(dat1>0,4);
nCC = dat1CC.NumObjects;
pixCC = dat1CC.PixelIdxList;

% which one -----
if p.crop==0
    idx8 = sub2ind(size(dat0),124,662);
    for ii=1:nCC
        pix0 = pixCC{ii};
        if sum( bsxfun(@eq,pix0,idx8) ) ==1
            fprintf('ii %d\n',ii);
        end
    end
end

% K10 = cat(3,dat0,dat1,dat1*0);
% figure;imshow(K10);

% dat1 = imsharpen(dat1,'Radius',2,'Amount',3);
% dat1(dat1>1) = 1;
% dat1(dat1<0) = 0;

%% detecting
resBeadx = {};
% for ii=1:100
% for ii=1023
% for ii=1:nCC
parfor ii=1:nCC
    fprintf('nCC: %d\n',ii);
    pix0 = pixCC{ii};
    dat1Ele = dat1*0;
    dat1Ele(pix0) = dat1(pix0);
    %     dat1Ele(pix0) = dat1(pix0);
    %     dat1Ele(pix0) = dat1(pix0);
    [rga, rgb] = myCropRg(pix0,Nx,Ny,20);
    d1 = dat1Ele(rga,rgb);  % as a mask
    d2 = dat0(rga,rgb);  % raw data
    resBead0 = beadMatchPairSeqMaskNoThr( d2, d1, myRef, rga(1), rgb(1), 0 );
%     resBead0 = beadMatchPairSeq( d1, myRef, rga(1), rgb(1), 0 );
    %     resBead0 = beadMatchPair( d1, myRef, rga(1), rgb(1), 0 );
    %     resBead0 = beadMatchCircle( d1, rga(1), rgb(1), 3:13, 0 );
    %     resBead0 = beadMatch( d1, rga(1), rgb(1), 5:7, 0 );
    resBeadx{ii} = resBead0;
end

resBead = {};
for ii=1:length(resBeadx)
    resBead = [resBead,resBeadx{ii}];
end

%% plot
resBorder = dat1*0;
neibVec = [0, -1, 1, -Nx, Nx];

fprintf('Dump ===== \n');
for ii=1:length(resBead)
    %     resBorder0 = dat1*0;
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBead{ii};
    idx = sub2ind(size(dat1),idx(:,1),idx(:,2));
    
%     idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
%     idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
%     idxc = setdiff(idxk,idx);
    
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);        

    resBorder(idxc) = 0.75;
    %     resBorder0(idxk) = 0.75;
    %     resBorder0(idx) = 0;
    %     resBorder = max(resBorder,resBorder0);
end

K2sqrt = cat(3,resBorder,sqrt(dat0),dat0*0);
K2 = cat(3,resBorder,dat0,dat0*0);
figure;imshow(K2);
fname = 'res_bead_mask_1_thr_0p08_score_0p1_avoid_1_0p4_rad_0_mean';
save([fname,'.mat'],'resBead','p');
imwrite(K2,[fname,'.tif']);
imwrite(K2sqrt,[fname,'_sqrtbg.tif']);



