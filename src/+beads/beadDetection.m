function [resBead,resBeadCenter,resBeadRad] = beadDetection(dat0,thr0,thrMulti,rad0,rad1,usePara)
% BEAD_DETECTION Detect bead based allow overlap
% dat0: matrix scaled to [0,1]
% resBead: cell vector with each element the coordinates of one bead

myRg = rad0:0.5:rad1;
% myRg = 4:0.5:13;
myRef = beadTemplate(myRg);

[Nx,Ny] = size(dat0);
dat0Mask = dat0>thr0;

dat0Mask = bwareaopen(dat0Mask, 10, 4);
dat1 = dat0.*dat0Mask;
dat1CC = bwconncomp(dat1>0,4);
nCC = dat1CC.NumObjects;
pixCC = dat1CC.PixelIdxList;

rMin = min(myRg);
rMax = max(myRg);
rPhi = exp(1i*2*pi*(log(myRg)-log(rMin))/(log(rMax)-log(rMin)));

if 0
    thr0 = 0.14;
    dat0Mask = dat0>thr0;    
    dat0Mask = bwareaopen(dat0Mask, 10, 4);
    dat1 = dat0.*dat0Mask;
    dat1CC = bwconncomp(dat1>0,4);
    nCC = dat1CC.NumObjects;
    pixCC = dat1CC.PixelIdxList;    
    cc=cellfun(@length,pixCC);
    [xC,xI]= max(cc);
    idx000 = pixCC{xI};
    tmp = dat0*0;
    tmp(idx000) =1;
    K0 = cat(3,tmp,dat0,dat0*0);
    figure;imshow(K0)
end

% detecting -----
resBeadx = {};
resBeadCenterx = {};
resBeadRadx = {};
if usePara
%     h = msgbox('Detecting beads...');
    parfor ii=1:nCC
        fprintf('nCC: %d\n',ii);
        pix0 = pixCC{ii};
        dat1Ele = dat1*0;
        dat1Ele(pix0) = dat1(pix0);
        [rga, rgb] = myCropRg(pix0,Nx,Ny,20);
        d1 = dat1Ele(rga,rgb);  % as a mask
        d2 = dat0(rga,rgb);  % raw data
        [resBead0,resBeadCenter0,resBeadRad0] = beadMatchMulti( d2, d1, myRef, rga(1), rgb(1), 0, thrMulti, rPhi);
        resBeadx{ii} = resBead0;
        resBeadCenterx{ii} = resBeadCenter0;
        resBeadRadx{ii} = resBeadRad0;
    end
%     close(h);
else
    % for ii=1:nCC
    % for ii=299
%     h = waitbar(0,'Detecting beads...');
    for ii=1:nCC
        fprintf('nCC: %d\n',ii);
        pix0 = pixCC{ii};
        dat1Ele = dat1*0;
        dat1Ele(pix0) = dat1(pix0);
        [rga, rgb] = myCropRg(pix0,Nx,Ny,20);
        d1 = dat1Ele(rga,rgb);  % as a mask
        d2 = dat0(rga,rgb);  % raw data
        %     resBead0 = beadMatchPairSeqMaskNoThr( d2, d1, myRef, rga(1), rgb(1), 0 );
        [resBead0,resBeadCenter0,resBeadRad0] = beadMatchMulti( d2, d1, myRef, rga(1), rgb(1), 0, thrMulti, rPhi);
        resBeadx{ii} = resBead0;
        resBeadCenterx{ii} = resBeadCenter0;
        resBeadRadx{ii} = resBeadRad0;
%         waitbar(ii/nCC);
    end
%     close(h)
end

resBead = {};
resBeadCenter = {};
resBeadRad = {};
for ii=1:length(resBeadx)
    resBead = [resBead,resBeadx{ii}];
    resBeadCenter = [resBeadCenter,resBeadCenterx{ii}];
    resBeadRad = [resBeadRad,resBeadRadx{ii}];
end

% % plot -----
% resBorder = dat1*0;
% neibVec = [0, -1, 1, -Nx, Nx];
% 
% fprintf('Dump ===== \n');
% for ii=1:length(resBead)
%     if mod(ii,1000)==0
%         fprintf('ii is %d\n',ii);
%     end
%     idx = resBead{ii};
%     idx = sub2ind(size(dat1),idx(:,1),idx(:,2));
%     idxk = bsxfun(@plus,idx,neibVec);
%     idxa = ismember(idxk,idx);
%     idxSel = sum(idxa,2)<5;
%     idxc = idx(idxSel);
%     resBorder(idxc) = 0.75;
% end
% 
% K2 = cat(3,resBorder,dat0,dat0*0);
% fname = 'res_bead';
% imwrite(K2,[outPath,filesep,fname,'.tif']);
% beadShow = K2;

end

