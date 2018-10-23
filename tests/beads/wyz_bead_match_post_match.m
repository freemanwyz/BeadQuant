%% post processing
p = getInitParamBead();
p.crop = 0;
plotAll = 1;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin))/255;
dat(dat<0.1) = 0;

fin = [p.tp filesep 'dump' filesep 'res_bead_match_20160623pm7.mat'];
resIn = load(fin);
resBeadIn = resIn.resBead;

% remove detected beads -----
datx = dat*0;
for ii=1:length(resBeadIn)
    idx = resBeadIn{ii};
    idx = sub2ind(size(dat),idx(:,1),idx(:,2));
    datx(idx) = 1;
end
datRemove = dat.*(1-datx);

% To isolated parts
thr0 = 0.1;
if p.crop
    dat0 = datRemove(p.crop_rgx,p.crop_rgy);
else
    dat0 = datRemove;
end
[Nx,Ny] = size(dat0);
dat0Mask = dat0>thr0;
dat0Mask = bwareaopen(dat0Mask, 15, 4);
dat1 = dat0.*dat0Mask;
dat1CC = bwconncomp(dat1>0,4);
nCC = dat1CC.NumObjects;
pixCC = dat1CC.PixelIdxList;

% which one -----
idx8 = sub2ind(size(dat0),636,2087);
for ii=1:nCC
    pix0 = pixCC{ii};
    if sum( bsxfun(@eq,pix0,idx8) ) ==1
        fprintf('ii %d\n',ii);
    end
end

%% detecting -----
resBeadResi = {};
resBead = {};
% for ii=1427
for ii=1:nCC
    fprintf('nCC: %d\n',ii);
    pix0 = pixCC{ii};
    datMask = dat1*0;    
    datMask(pix0) = 1;
    [rgx, rgy] = myCropRg(pix0,Nx,Ny,10);
    d1 = dat(rgx,rgy);
    d1Mask = datMask(rgx,rgy);
    resBead0 = {};
    if length(pix0) > 150
        resBead0 = beadMatchPostBig( d1, d1Mask, rgx(1), rgy(1), 9:14, 0 );        
    end
    if isempty(resBead0)        
        resBead0 = beadMatchPost( d1, d1Mask, rgx(1), rgy(1), 3:5, 0 );
    end
    if isempty(resBead0) && length(pix0) < 75
        resBead0 = beadMatchPostSmall( d1, d1Mask, rgx(1), rgy(1));
    end
    resBead = [resBead,resBead0];
end
resBeadResi = [resBeadResi, resBead];

% plot -----
resBorder = dat*0;
neibVec = [0, -1, 1, -Nx, Nx];

fprintf('Dump ===== \n');
for ii=1:length(resBeadResi)
    resBorder0 = dat1*0;
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBeadResi{ii};
    idx = sub2ind(size(dat1),idx(:,1),idx(:,2));
    idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
    idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
    resBorder0(idxk) = 0.75;
    resBorder0(idx) = 0;
    resBorder = max(resBorder,resBorder0);
end

daty = dat;
% daty(dat1==0) = daty(dat1==0)*0.4;

K2 = cat(3,resBorder,daty,dat*0);
figure;imshow(K2);

K2a = cat(3,resBorder,dat1,dat*0);
figure;imshow(K2a);

%% plot all
if plotAll
    resBeadAll = [resBeadIn,resBeadResi];
    resBorderAll = dat*0;
    neibVec = [0, -1, 1, -Nx, Nx];
    
    fprintf('Dump ===== \n');
    for ii=1:length(resBeadAll)
        resBorder0 = dat*0;
        if mod(ii,1000)==0
            fprintf('ii is %d\n',ii);
        end
        idx = resBeadAll{ii};
        idx = sub2ind(size(dat),idx(:,1),idx(:,2));
        idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
        idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
        resBorder0(idxk) = 0.75;
        resBorder0(idx) = 0;
        resBorderAll = max(resBorderAll,resBorder0);
    end
    
    K3 = cat(3,resBorderAll,dat,dat*0);
    figure;imshow(K3);    
end






