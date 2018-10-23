%% Bead signal from grace
p = getInitParamBead();
p.crop_rgx = 651:800;
p.crop_rgy = 1051:1150;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin));
fin = [p.tp filesep 'dump' filesep 'beadCorMap.tif'];
corMap = double(imread(fin))/255;

%% To isolated parts
datOrg = dat(p.crop_rgx,p.crop_rgy)/255;
% datOrg = dat/255;
corMapOrg = corMap(p.crop_rgx,p.crop_rgy);
% corMapOrg = corMap;
dat0 = datOrg;
% dat0 = wiener2(dat0,[7 7]);
[Nx,Ny] = size(dat0);
resBead = zeros(size(dat0));
thr0 = 0.1;
dat0Mask = dat0>thr0;
dat0Mask = bwareaopen(dat0Mask, 10);
dat1 = dat0.*dat0Mask;
dat0CC = bwconncomp(dat0Mask,4);
nCC = dat0CC.NumObjects;
pixCC = dat0CC.PixelIdxList;

% idx8 = sub2ind(size(dat0),494,415);
% for ii=1:nCC
%     pix0 = pixCC{ii};
%     if sum( bsxfun(@eq,pix0,idx8) ) ==1
%         fprintf('ii %d\n',ii);
%     end
% end

%% detecting
for ii=5
    % for ii=1:nCC
    done = 0;
    % crop a part out of the large image for faster computation -----
    pix0 = pixCC{ii};
    dat1Ele = dat1*0;
    dat1Ele(pix0) = dat1(pix0);
    [ai,bi] = find(dat1Ele>0);
    rga = (min(ai)-5):(max(ai)+5);
    rgb = (min(bi)-5):(max(bi)+5);
    rga = rga(rga>0 & rga<Nx);
    rgb = rgb(rgb>0 & rgb<Ny);
    d1org = dat1Ele(rga,rgb);
    corMap1 = corMapOrg(rga,rgb);
    
    % de-blurring -----
    d1orga = d1org;
    d1orga(d1orga>0.4) = 0.4;
    PSF = fspecial('gaussian',5,5);
    V = .0001;
    INITPSF = ones(size(PSF));
    WT = zeros(size(d1org));
    WT(d1org>0.1) = 1;
    J = deconvblind(d1orga,INITPSF,20,10*sqrt(V),WT);
    J(J<0.1) = 0;
    %     figure;imshow(J);
    
    % bourndary mask ----
    baseMask0 = d1org*0;
    B = bwboundaries(d1org>0);
    for nn = 1:length(B)
        B1 = B{nn};
        idxB = sub2ind(size(d1org),B1(:,1),B1(:,2));
        baseMask0(idxB) = 1;
    end
    %     d1edge = edge(d1org,'Sobel',0.1);
    %     baseMask(d1edge>0) = 1;
    baseMask = baseMask0;
    d1orgShrink = d1org;
    d1orgShrink(baseMask>0) = 0;
    B = bwboundaries(d1orgShrink>0);
    for nn = 1:length(B)
        B1 = B{nn};
        idxB = sub2ind(size(d1org),B1(:,1),B1(:,2));
        baseMask(idxB) = 1;
    end
    
    res = struct('fig1',{},'fig2',{},'ratio',{},'beads',{},'nBeads',{});
    % watershed + circle -----
    d1 = d1org;
    %     d1(d1>0.4) = 0.4;
    resBead0 = zeros(size(d1));
    thr1rg = [0.03 0.3];
    for jj=1:length(thr1rg)
        thr0 = thr1rg(jj);
        if jj>1
            idxUsed = resBead0>0;
            d1(idxUsed>0) = 0;
        end
        if jj==1  % watershed
            resBead0 = beadEdgeWaterShed( d1,resBead0,thr0 );
        else  % draw circle
            resBead0 = beadEdgeCircle( d1,resBead0,thr0,baseMask );
        end
    end
    res(1) = beadEval( d1org,baseMask,resBead0 );
    
    % circle with deblur -----
    resBead0 = beadCircle( J,zeros(size(d1)),baseMask);
    res(2) = beadEval( d1org,baseMask,resBead0 );
    
    %     % circle seq with deblur and dilation -----
    %     resBead0 = beadCircleSeq( J,zeros(size(d1)),1);
    %     res(5) = beadEval( d1org,baseMask,resBead0 );
    
    % circle seq with deblur -----
    resBead0 = beadCircleSeq( J,zeros(size(d1)),0,baseMask);
    res(3) = beadEval( d1org,baseMask,resBead0 );
    
    % circle -----    
    resBead0 = beadCircle( d1org,zeros(size(d1)),baseMask);
    res(4) = beadEval( d1org,baseMask,resBead0 );
    
    %     % circle with multiple thresholds -----
    %     resBead0 = beadCircleSeq( d1org,zeros(size(d1)),1);
    %     res(6) = beadEval( d1org,baseMask,resBead0 );
    %
    %     % circle with multiple thresholds -----
    %     resBead0 = beadCircleSeq( d1org,zeros(size(d1)),0);
    %     res(7) = beadEval( d1org,baseMask,resBead0 );
    
    % choose best one
    %     [C,I] = max( [res.ratio] );
    dat = d1org;
    [C,I] = max( [res.ratio] - [res.nBeads]*5 );
    fprintf('nCC: %d Best: %d\n',ii,I);
    resBeadSel = res(I).beads;
    nBeadNow = max(resBead(:));
    resBeadSel(resBeadSel>0) = resBeadSel(resBeadSel>0) + nBeadNow;
    resBead(rga,rgb) = max(resBeadSel,resBead(rga,rgb));
end

% %% bead mask
% myTemplate = zeros(11,11);
% idxBead = unique(resBead);
% idxBead = idxBead(idxBead>0);
% beadCnt = 0;
% for jj=1:length(idxBead)
%     tmp = resBead == idxBead(jj);
%     stat0 = regionprops(tmp,'Centroid');
%     xy1 = round(stat0.Centroid);
%     x1 = xy1(1);
%     y1 = xy1(2);
%     if x1>7 && x1<Nx-7 && y1>7 && y1<Ny-7
%         beadCnt = beadCnt + 1;
%         myTemplate = myTemplate + datOrg( (x1-5):(x1+5), (y1-5):(y1+5) );
%     end
% end
% myTemplate = myTemplate/beadCnt;

%% plot
idxBead = unique(resBead);
idxBead = idxBead(idxBead>0);

resBorder = resBead*0;
neibVec = [0, -1, 1, -Nx, Nx];

fprintf('Dump ===== \n');
for ii=1:length(idxBead)
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = find(resBead==idxBead(ii));
    idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
    idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
    resBorder(idxk) = 0.75;
    resBorder(idx) = 0;
end

K2 = cat(3,resBorder,datOrg,dat1*0);
figure;imshow(K2);





