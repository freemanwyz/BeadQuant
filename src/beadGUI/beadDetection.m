function [resBead,resBead1,resBeadCenter,resBeadRad,beadShow] = beadDetection(dat0,...
    smoMethod,smoSize,accumThresh,rad0,rad1,outPath)
% BEAD_DETECTION Detect bead based allow overlap
% dat0: matrix scaled to [0,1]
% resBead: cell vector with each element the coordinates of one bead

% detection -----
dat1 = sqrt(dat0);
if smoMethod==1
    if mod(smoSize,2)==0
        smoSize = smoSize + 1;
    end
    if smoSize>=3
        dat1 = wiener2(dat1,[smoSize,smoSize]);
    end
else
    if smoSize>0
        dat1 = imgaussfilt(dat1,smoSize);
    end
end
[resBeadCenter, resBeadRad, resBead, resBead1] = seqCHT(dat1,[rad0 rad1],accumThresh);

% test -----
% remove some beads for better comparison with manual counting
if 0
    [Nx,Ny] = size(dat1);
    beadValid = zeros(length(resBead),1);
    for ii=1:size(resBeadCenter,1)
        idx0 = resBead{ii};
        % remove dark beads
        idx0 = sub2ind([Nx,Ny],idx0(:,1),idx0(:,2));
        if mean(dat0(idx0))<0  % 0.02
            continue
        end
        
        % remove beads close to image border
        xx = resBeadCenter(ii,:);
        distB = 2;
        if xx(1)<distB || xx(1)>(Nx-distB+1) || xx(2)<distB || xx(2)>(Ny-distB+1)
            continue
        end
        beadValid(ii) = 1;
    end
    
    resBeadCenter = resBeadCenter(beadValid>0,:);
    resBeadRad = resBeadRad(beadValid>0,:);
    resBead = resBead(beadValid>0);
    resBead1 = resBead1(beadValid>0);
    
    KCenter = zeros(size(dat1));
    for ii=1:size(resBeadCenter,1)
        % remove beads close to image border
        xx = resBeadCenter(ii,:);
        KCenter(xx(1),xx(2)) = KCenter(xx(1),xx(2)) + 1;
    end
    K1 = cat(3,KCenter,dat0,dat0*0);
    fname = 'res_bead_center';
    imwrite(double(K1),[outPath,fname,'.tif']);
end

% plot -----
[Nx,~] = size(dat1);
resBorder = dat0*0;
neibVec = [0, -1, 1, -Nx, Nx];

fprintf('Dump ===== \n');
for ii=1:length(resBead)
    if mod(ii,1000)==0
        fprintf('ii is %d\n',ii);
    end
    idx = resBead{ii};
    idx = sub2ind(size(dat0),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);
    resBorder(idxc) = 0.75;
end

K2 = cat(3,resBorder,dat0,dat0*0);
fname = 'res_bead';
imwrite(double(K2),[outPath,fname,'.tif']);
% imwrite(K2,[outPath,filesep,fname,'.tif']);
beadShow = double(K2);

end

