%% post processing
p = getInitParamBead();
p.crop = 0;

% load data -----
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin));

fin = [p.tp filesep 'dump' filesep 'res_bead_match_20160623pm7.mat'];
resIn = load(fin);
resBeadIn = resIn.resBead;

% preparation ------
resBeadToRemove = resBeadIn;

thr0 = 0.1;
% remove detected beads -----
datx = dat*0;
for ii=1:length(resBeadToRemove)
    idx = resBeadToRemove{ii};
    idx = sub2ind(size(dat),idx(:,1),idx(:,2));
    datx(idx) = 1;
end
% datx = imdilate(datx,strel('square',3));
datb = dat.*(1-datx);

% To isolated parts
if p.crop
    dat0 = datb(p.crop_rgx,p.crop_rgy)/255;
else
    dat0 = datb/255;
end
[Nx,Ny] = size(dat0);
dat0Mask = dat0>thr0;
dat0Mask = bwareaopen(dat0Mask, 15, 4);
dat1 = dat0.*dat0Mask;
dat1CC = bwconncomp(dat1>0,4);
nCC = dat1CC.NumObjects;
pixCC = dat1CC.PixelIdxList;

% %% which one
idx8 = sub2ind(size(dat0),1462,2030);
for ii=1:nCC
    pix0 = pixCC{ii};
    if sum( bsxfun(@eq,pix0,idx8) ) ==1
        fprintf('ii %d\n',ii);
    end
end

%% detecting -----
resBeadResi = {};
resBead = {};
useEdge = 1;
for ii=1388
% for ii=1:nCC
    fprintf('nCC: %d\n',ii);
    pix0 = pixCC{ii};
    dat1Ele = dat1*0;
    dat1Ele(pix0) = dat1(pix0);
    [rgx, rgy] = myCropRg(pix0,Nx,Ny);
    d1a = dat1Ele(rgx,rgy);
    d1 = sqrt(d1a);
    
    % select some neighbors -----
    d1bw = 1*(d1>0);
    d1bwDilate = imdilate(d1bw,p.SE);
    d1bwDif0 = d1bwDilate - d1bw;
    if useEdge
        d1edge = edge(d1,'Sobel',0.1);
        d1bwDif = max(d1bwDif0,d1edge);
        %             d1bwDif = d1bwDif0;
    else
        d1bwDif = d1bwDif0;
    end
    
    % compute the min distance to the neighbor -----
    [ai,aj] = find(d1>0);
    nCen = length(ai);
    [bi,bj] = find(d1bwDif>0);
    xneib = [bi,bj];
    dist0 = d1*0;
    for uu=1:nCen
        dist00 = sqrt((bi - ai(uu)).^2 + (bj-aj(uu)).^2);
        dist0(ai(uu),aj(uu)) = min(dist00);
    end
    
    [ai,aj] = find(d1>0);
    nCen = length(ai);
    [bi,bj] = find(d1bwDif0>0);
    xneib = [bi,bj];
    dist0a = d1*0;
    for uu=1:nCen
        dist00 = sqrt((bi - ai(uu)).^2 + (bj-aj(uu)).^2);
        dist0a(ai(uu),aj(uu)) = min(dist00);
    end
    
    % get a circle for each local maximum -----
    iCnt = 1;
    resBead00 = {};
%     resBead01 = {};
    dist0LMa = imregionalmax(dist0);
    dist0LMb = imregionalmax(dist0a);
    dist0LM = max(dist0LMa,dist0LMb);
    while sum(dist0LM(:))>0
        cIdx = find(dist0LM>0);
        [C,I] = max(dist0a(cIdx));
        [ci,cj] = ind2sub(size(dist0),cIdx(I));
        [di,dj] = find(dist0 >- 1);
        dist1 = sqrt((di - ci).^2 + (dj-cj).^2);
        if C>=3 && sum(dist1<=C)>15;
            exy = find(dist1<=C);
%             exy1 = find(dist1<=(C-1));
            resBead00{iCnt} = exy;
%             resBead01{iCnt} = exy1;
            iCnt = iCnt + 1;
        end
        dist0LM(dist1<=1) = 0;
    end
    
    % rank by homogeneity -----
    var0 = zeros(1,length(resBead00));
    for nn=1:length(resBead00)
        PixMe = resBead00{nn};
        var0(nn) = var(d1(PixMe))/length(PixMe);
%         tmp = d1*0;
%         tmp(PixMe) = 1;
%         K02 = cat(3,tmp,d1,d1*0);
%         figure(nn);imshow(K02);
    end
    
    iCnt = 1;
    resBead0 = {};
    [Cvar,Ivar] = sort(var0);    
    for nn=1:length(resBead00)
        PixMe = resBead00{Ivar(nn)};
        if iCnt==1
            if length(PixMe)>15
                resBead0{iCnt} = PixMe;
                iCnt = iCnt + 1;
            end
        else
            for mm=1:(iCnt-1)
                PixYou = resBead0{mm};
                [C,ia,ib] = intersect(PixMe,PixYou);
                PixMe(ia) = 0;
            end
            if sum(PixMe==0)/length(PixMe)<0.4 && length(PixMe)>15
                PixMe = resBead00{Ivar(nn)};                
                resBead0{iCnt} = PixMe;
                iCnt = iCnt + 1;
            end
        end
    end
    
    resBead0a = {};
    if iCnt>1
        for nn=1:length(resBead0)
            PixMe = resBead0{nn};
            [ex,ey] = ind2sub(size(d1),PixMe);
            resBead0a{nn} = [ex+rgx(1)-1,ey+rgy(1)-1];
        end
        resBead = [resBead,resBead0a];
    end
end
resBeadResi = [resBeadResi, resBead];


%% plot
resBorder = dat*0;
neibVec = [0, -1, 1, -Nx, Nx];

datx = dat*0;
for ii=1:length(resBeadIn)
    idx = resBeadIn{ii};
    idx = sub2ind(size(dat),idx(:,1),idx(:,2));
    datx(idx) = 1;
end
% datx = imdilate(datx,strel('square',3));
datb = dat.*(1-datx);

% To isolated parts
if p.crop
    dat0 = datb(p.crop_rgx,p.crop_rgy)/255;
else
    dat0 = datb/255;
end
[Nx,Ny] = size(dat0);
dat0Mask = dat0>0.1;
dat0Mask = bwareaopen(dat0Mask, 15, 4);
dat1 = dat0.*dat0Mask;

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

K2 = cat(3,resBorder,dat1,dat1*0);
figure;imshow(K2);

%% plot all
resBeadAll = [resBeadIn,resBead];
fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat1 = double(imread(fin))/256;
resBorderAll = dat1*0;
neibVec = [0, -1, 1, -Nx, Nx];

fprintf('Dump ===== \n');
for ii=1:length(resBeadAll)
    resBorder0 = dat1*0;
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBeadAll{ii};
    idx = sub2ind(size(dat1),idx(:,1),idx(:,2));
    idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
    idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
    resBorder0(idxk) = 0.75;
    resBorder0(idx) = 0;
    resBorderAll = max(resBorderAll,resBorder0);
end

K3 = cat(3,resBorderAll,dat1,dat1*0);
figure;imshow(K3);






