function [centers, r_estimated, resBead, resBead1, centersAll, r_estimatedAll, metricAll, accumMatrixAll] = seqCHT(A,radiusRange,accumThresh)
%seqCHT Find circles using Circular Hough Transform iteratively

%% setup 
[mA,nA] = size(A);
mask = ones(size(A),'single');
maskCenter = ones(size(A),'single');
% accumThresh = 0.05;
rRg = radiusRange(1):0.5:radiusRange(2);
sqSize = round((radiusRange(2)+5)/2)*2+1;
sq1 = (sqSize-1)/2 + 1;
rRg1 = radiusRange(1):0.5:(radiusRange(2)+1);
cirs = beadTemplate( rRg1, sqSize );

dbg = 0;
removeMode = 0;  % 0: all, 1: boundary only
nTry = 3;
minMetric = 0;
centersAll = cell(1,nTry);
r_estimatedAll = cell(1,nTry);
metricAll = cell(1,nTry);
accumMatrixAll = cell(1,nTry);

%% detection
h = waitbar(0,'Detecting beads...');
for ii=1:nTry
    fprintf('nTry: %d\n',ii);
    % Compute the accumulator array
    [accumMatrix, ~] = chaccum(A, radiusRange, mask);
    
    % dbg >>>>>
    if ii==1 && dbg==1
        % [672,1846] looks strange in chcenters
        keyboard
        a0 = accumMatrix;
        K0 = cat(3,a0/max(a0(:)),A*0.5,A*0);
        figure;imshow(K0);
    end
    % dbg <<<<<
    
    % Estimate the centers
    [centers, metric] = chcenters(accumMatrix, accumThresh);
    
    % remove centers in detected circles
    centersIdx = sub2ind(size(A),round(centers(:,2)),round(centers(:,1)));
    centersVal = maskCenter(centersIdx)==1;
    centers = centers(centersVal,:);
    metric = metric(centersVal);
    
    % dbg >>>>>
    if ii==1 && dbg==1
        keyboard
        tmp = A*0;
        tmp(sub2ind(size(A),round(centers(:,2)),round(centers(:,1)))) = metric/max(metric(:));
        K1 = cat(3,tmp,A,A*0);
        figure;imshow(K1);
    end
    % dbg <<<<<
    
    % Retain circles with metric value greater than threshold corresponding to AccumulatorThreshold
    fprintf('%f\n',mean(metric));
    thr1 = max(mean(metric)/2,minMetric);
    idx2Keep = find(metric >= thr1);
    centers = centers(idx2Keep,:);
    metric = metric(idx2Keep);
    
    % Estimate radii
    r_estimated = chradiiphcode(centers, accumMatrix, radiusRange);
    
    % dbg >>>>>
    if dbg==1
        K0 = cat(3,A*0,A,A*0);
        figure;imshow(K0);
        viscircles(centers, r_estimated,'EdgeColor','b');
        title(['nTry',num2str(ii)]);
    end
    % dbg <<<<<
    
    % update mask for bead center
    % center with higher metric will remove centers within it
    centerVal = ones(length(metric),1);
    for jj=1:length(r_estimated)
        c1 = round(centers(jj,:));        
        if maskCenter(c1(2),c1(1))==1
            r1 = round(r_estimated(jj)*2)/2;
            if r1>=rRg(end)*0.7 && ii>nTry/2
                centerVal(jj) = 0;
            else
                idx1 = find(rRg==r1);                
                t1 = cirs(idx1(1)).maskCir;
                x1 = reshape(t1(:,1) - sq1 + c1(2),[],1);
                y1 = reshape(t1(:,2) - sq1 + c1(1),[],1);                
                idxVal = x1 > 0 & x1 < mA & y1 >0 & y1 < nA;
                x1 = x1(idxVal);
                y1 = y1(idxVal);
                maskCenter(sub2ind(size(A),x1,y1)) = 0;
            end
        else
            centerVal(jj) = 0;
        end      
    end
    
    centers = centers(centerVal>0,:);
    r_estimated = r_estimated(centerVal>0);
    metric = metric(centerVal>0);
    
    % update mask for gradient map
    if removeMode==0
        mask = maskCenter;
    else
        for jj=1:length(r_estimated)
            c1 = round(centers(jj,:));
            r1 = round(r_estimated(jj)*2)/2;
            idx1 = find(rRg==r1);
            
            tx = cirs(idx1(1)).maskX;
            ty = cirs(idx1(1)).maskY;
            x1 = reshape(tx - sq1 + c1(2),[],1);
            y1 = reshape(ty - sq1 + c1(1),[],1);
            
            idxVal = x1 > 0 & x1 < mA & y1 >0 & y1 < nA;
            x1 = x1(idxVal);
            y1 = y1(idxVal);
            mask(sub2ind(size(A),x1,y1)) = 0;
        end
    end        
    
    % gather results
    centersAll{ii} = centers;
    r_estimatedAll{ii} = r_estimated;
    metricAll{ii} = metric;
    accumMatrixAll{ii} = accumMatrix;
    
    h = waitbar(ii/nTry);
end

if isvalid(h)
    close(h)
end

%% gather
nBeads = 0;
for ii=1:length(r_estimatedAll)
    nBeads = nBeads + length(r_estimatedAll{ii});
end

centers = zeros(nBeads,2);
r_estimated = zeros(nBeads,1);
nCnt = 1;
for ii=1:length(r_estimatedAll)
    c0 = centersAll{ii};
    r0 = r_estimatedAll{ii};
    centers(nCnt:(nCnt+length(r0)-1),:) = c0;
    r_estimated(nCnt:(nCnt+length(r0)-1)) = r0;
    nCnt = nCnt + length(r0);
end

resBead = cell(length(r_estimated),1);
for jj=1:nBeads
    c1 = round(centers(jj,:));
    r1 = round(r_estimated(jj)*2)/2;
    idx1 = find(rRg==r1);
    t1 = cirs(idx1(1)).maskCir;
%     t1 = cirs(idx1(1)).maskCir0;
    x1 = reshape(t1(:,1) - sq1 + c1(2),[],1);
    y1 = reshape(t1(:,2) - sq1 + c1(1),[],1);
    idxVal = x1 > 0 & x1 < mA & y1 >0 & y1 < nA;
    x1 = x1(idxVal);
    y1 = y1(idxVal);
    resBead{jj} = [x1,y1];
end

% enlarge by one -----
resBead1 = cell(length(r_estimated),1);
for jj=1:nBeads
    c1 = round(centers(jj,:));
    r1 = round(r_estimated(jj)*2)/2;
    idx1 = find(rRg==r1);
    t1 = cirs(idx1(1)+1).maskCir;
    x1 = reshape(t1(:,1) - sq1 + c1(2),[],1);
    y1 = reshape(t1(:,2) - sq1 + c1(1),[],1);
    idxVal = x1 > 0 & x1 < mA & y1 >0 & y1 < nA;
    x1 = x1(idxVal);
    y1 = y1(idxVal);
    resBead1{jj} = [x1,y1];
end

centers = round([centers(:,2),centers(:,1)]);
r_estimated = round(r_estimated);

end





