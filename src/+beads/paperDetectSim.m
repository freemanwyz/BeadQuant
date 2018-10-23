function [ resPrecision,resRecall,resMSE ] = paperDetectSim( pathTop,pathLbTop,...
    pathOutTop,thrMulti,thrDetect,radRg0,radRg1,Nimg,plotme,noiseAmp )
%beadPaperDetectSim run a single simulation

smoMethod = 1;
resPrecision = zeros(1,Nimg);
resRecall = zeros(1,Nimg);
resMSE = zeros(1,Nimg);

for ii=1:Nimg
    % cropped orignal image
    fname = [num2str(ii),'.tif'];
    fpname = [pathTop,fname];
    X0 = tiffread(fpname);
    bgMean = double(X0(1).data)./65535;
    
    if noiseAmp>0
        bgMean = sqrt(bgMean) + randn(size(bgMean))*noiseAmp;
        bgMean(bgMean<0) = 0;
        bgMean(bgMean>1) = 1;
        bgMean = bgMean.^2;
    end
    
    % labelling
    fname = [num2str(ii),'_labeled.tif'];
    fpname = [pathLbTop,filesep,fname];
    XLb = imread(fpname);
    bgCenter = double(XLb(:,:,3))./255;
    bgCenter = (bgCenter>0.5)*1;
    bwcc = bwconncomp(bgCenter,4);
    nBeadLb = bwcc.NumObjects;
    
    % detection
    outPath = [pathOutTop num2str(ii) '_'];
    [resBead,~,resBeadCenter,~,~] = beadDetection(bgMean,...
        smoMethod,thrMulti,thrDetect,radRg0,radRg1,outPath);
    nBeadDetect = length(resBead);
    resCenter1 = bgMean*0;
    idx1 = sub2ind(size(bgMean),resBeadCenter(:,1),resBeadCenter(:,2));
    resCenter1(idx1) = 1;
    
    % distance, sensitivity and precision
    beadHit = zeros(nBeadDetect,1);
    distHit = zeros(nBeadDetect,1);
    for jj=1:nBeadDetect
        tmpMask = bgMean*0;
        c0 = resBeadCenter(jj,:);
        idxSub = resBead{jj};
        idxInd = sub2ind(size(bgMean),idxSub(:,1),idxSub(:,2));
        tmpMask(idxInd) = 1;
        [cx,cy] = find(tmpMask.*bgCenter > 0);
        if ~isempty(cx)
            beadHit(jj) = 1;
            if length(cx)>1
                dist01 = (cx-c0(1)).^2 + (cy-c0(2)).^2;
                [~,I] = min(dist01);
                cx = cx(I);
                cy = cy(I);
            end
            distHit(jj) = (cx-c0(1)).^2 + (cy-c0(2)).^2;
        end
    end
    resPrecision(ii) = sum(beadHit)/nBeadDetect;
    resRecall(ii) = sum(beadHit)/nBeadLb;
    if resRecall(ii)>1
        resRecall(ii) = 1;
    end
    resMSE(ii) = mean(distHit(beadHit>0));
    
    if plotme
        K0 = cat(3,resCenter1,sqrt(bgMean),bgCenter);
        figure;imshow(K0);
    end    
end

end

