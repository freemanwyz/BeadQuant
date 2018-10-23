%% Bead signal from grace
% bad !!
p = getInitParamBead();
p.crop_rgx = 651:800;
p.crop_rgy = 1051:1150;

fin = [p.tp filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128' filesep ...
    '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128.lsm - Ch1-T1 - C1 Z1 T'];

% nFrame = 32;
% for ii=1:nFrame
%     fin0 = [fin num2str(ii) '.tif'];
%     dat0 = imread(fin0);
%     if ii==1
%         dat = double(dat0);
%     else
%         dat = dat + double(dat0);
%     end
% end
% dat = dat/nFrame;

fin = [p.tp filesep 'dump' filesep 'beadAvg.tif'];
dat = double(imread(fin));

%% To isolated parts
datOrg = dat(p.crop_rgx,p.crop_rgy)/255;
% datOrg = dat/65535;
dat0 = datOrg;
[Nx,Ny] = size(dat0);
resBead = zeros(size(dat0));
% thr0rg = [0.1 0.2 0.3 0.075];
thr0rg = 0.1;
% K = wiener2(dat0,[5 5]);
% for thr=[0.2,0.1]
tt = 1;
for tt=1:length(thr0rg)
    thr0 = thr0rg(tt);
    dat0Mask = dat0>thr0;
    dat0Mask = bwareaopen(dat0Mask, 10);
    dat1 = dat0.*dat0Mask;
    dat0CC = bwconncomp(dat0Mask,4);
    nCC = dat0CC.NumObjects;
    pixCC = dat0CC.PixelIdxList;
%     for ii=5
%     ii = 5;
    for ii=1:nCC
        fprintf('nCC: %d\n',ii);
        %         if ii==110 && nthr==2
        %             keyboard
        %         end
        % crop a part out of the large image for faster computation -----
        pix0 = pixCC{ii};
        dat1Ele = dat1*0;
        dat1Ele(pix0) = dat1(pix0);
        
        [ai,bi] = find(dat1Ele>0);
        rga = (min(ai)-5):(max(ai)+5);
        rgb = (min(bi)-5):(max(bi)+5);
        rga = rga(rga>0 & rga<Nx);
        rgb = rgb(rgb>0 & rgb<Ny);
        d1 = dat1Ele(rga,rgb);
        %         imshow(d1)
        %         keyboard
        
        for kk=1:100
            d1test = bwareaopen(d1, 10);
            if max(d1test(:))==0
                break
            end
            [Cmax,Imax] = max(d1(:));
            d1thr = d1*0;
            d1thr(d1>(Cmax*0.5)) = d1(d1>(Cmax*0.5));
%             imshow(d1thr);
            
            % select some neighbors -----
            d1bw = 1*(d1thr>0);
            d1bwDilate = imdilate(d1bw,p.SE);
            d1bwDif = d1bwDilate - d1bw;
            
            % compute the min distance to the neighbor -----
            [ai,aj] = find(d1thr>0);
            nCen = length(ai);
            [bi,bj] = find(d1bwDif>0);
            xneib = [bi,bj];
            dist0 = d1thr*0;
            for uu=1:nCen
                dist00 = sqrt((bi - ai(uu)).^2 + (bj-aj(uu)).^2);
                dist0(ai(uu),aj(uu)) = min(dist00);
            end
            
            % get a circle for each local maximum -----
            iCnt = 1;
            availMap = 1*(d1thr>0);
            resBead0 = zeros(size(d1thr));
            dist0LM = imregionalmax(dist0);
            while sum(dist0LM(:))>0
                cIdx = find(dist0LM>0);
                [C,I] = max(dist0(cIdx));
                [ci,cj] = ind2sub(size(dist0),cIdx(I));
                [di,dj] = find(dist0 >- 1);
                dist1 = sqrt((di - ci).^2 + (dj-cj).^2)<=C;
                if C>2 && sum(availMap(dist1))>25;
                    resBead0(dist1) = iCnt;
                    iCnt = iCnt + 1;
                end
                availMap(dist1) = 0;
                dist0LM(dist1) = 0;
            end
            
%             if (sum(resBead0(:)>0)/sum(d1bw(:)) < 0.75) && tt<4
%                 resBead0 = resBead0*0;
%             else
%                 dat0(pix0) = 0;
%             end
            dat0(pix0) = 0;
            d1(d1thr>0) = 0;
            
            nBeadNow = max(resBead(:));
            resBead0(resBead0>0) = resBead0(resBead0>0) + nBeadNow;
            resBead(rga,rgb) = max(resBead0,resBead(rga,rgb));
        end
%         keyboard
    end
end

%% plot
% K1 = cat(3,resBead,dat1,dat1*0);
% figure;imshow(K1);

idxBead = unique(resBead);
idxBead = idxBead(idxBead>0);

resBorder = resBead*0;
neibVec = [0, -1, 1, -Nx, Nx];

% idx = find(resBead>0);
% idxk = reshape(bsxfun(@plus,idx,neibVec),[],1);
% idxk = idxk(idxk>0 & idxk<=(Nx*Ny));
% resBorder(idxk) = 0.75;
% resBorder(idx) = 0;

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


% figure;imshow(resBead);
% figure;imshow(resBorder);








