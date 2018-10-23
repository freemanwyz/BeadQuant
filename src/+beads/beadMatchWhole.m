function resBead = beadMatchWhole( datIn, myRef, plotme )
% BEADMATCHPAIRSEQ Compute border pair score and sequential removal

thrx = 0.1;
[Nx,Ny] = size(datIn);
gap = 20;
datInAvoid = zeros(Nx,Ny)-1;
datInAvoid( (gap+1):(Nx-gap), (gap+1):(Ny-gap)) = 0;  % only use center part
radRg = [myRef.radius];
% ofst0 = sub2ind(size(datIn),gap,gap);

resBead = {};
beadCnt = 1;
nTry = 1;
while 1
    fprintf('Ntry: %d\n',nTry);
    nTry = nTry + 1;
    % scan radius -----
    [canX,canY] = find(datInAvoid==0);
    canX = reshape(canX,1,[]);
    canY = reshape(canY,1,[]);
    canXY = sub2ind(size(datIn),canX,canY);
    nCan = length(canX);
    matchMap = zeros(Nx,Ny,length(radRg));
    
    parfor ii=1:length(radRg)
        refX0 = myRef(ii).maskX;
        refY0 = myRef(ii).maskY;
        Cs = zeros(Nx,Ny);
        canX0 = canX-gap-1;
        canY0 = canY-gap-1;
        idx1x = bsxfun(@plus,refX0(:,1),canX0);
        idx1y = bsxfun(@plus,refY0(:,1),canY0);
        idx1 = reshape(sub2ind(size(datIn),idx1x(:),idx1y(:)),[],nCan);
        idx2x = bsxfun(@plus,refX0(:,2),canX0);
        idx2y = bsxfun(@plus,refY0(:,2),canY0);
        idx2 = reshape(sub2ind(size(datIn),idx2x(:),idx2y(:)),[],nCan);        
        idx3x = bsxfun(@plus,refX0(:,3),canX0);
        idx3y = bsxfun(@plus,refY0(:,3),canY0);
        idx3 = reshape(sub2ind(size(datIn),idx3x(:),idx3y(:)),[],nCan);      
        
        dif23 = datIn(idx2) - datIn(idx3);
        idxAvoid = datInAvoid(idx1)==1;
        dif23(idxAvoid) = NaN;
%         score0 = zeros(1,nCan);
%         for jj=1:nCan
%             xx = dif23(:,jj);
%             xx = xx(~isnan(xx));
%             score0(jj) = median(xx);            
%         end        
        score0 = nanmedian(dif23,1);
        nAvoid = sum(idxAvoid,1);
        score0( nAvoid>=size(refX0,1)*0.3 ) = 0;
        Cs(canXY) = score0;
        Cs(Cs<0) = 0;
        matchMap(:,:,ii) = Cs;
    end
    
    % find best one -----
    [matchMapBest,radBest] = max(matchMap,[],3);
    matchCenter = (imregionalmax(matchMapBest)>0 & matchMapBest>thrx) .* matchMapBest;
    [X,Y] = meshgrid(1:Nx,1:Ny);
    X = X.';
    Y = Y.';
    
    % dump -----
    if plotme
        figure;contour(flipud(matchMapBest),'ShowText','on');title(num2str(nTry));
        K0 = cat(3,matchCenter>0,datIn,datIn*0);
        figure;imshow(K0);
        K1 = cat(3,radBest/10.*(matchCenter>0),datIn,datIn*0);
        figure;imshow(K1);
    end
    
    if max(matchCenter(:))>thrx
        % if having the same score, choose the one with largest radius ----
        [xC,xI] = max(matchCenter(:)); 
        xI = find(matchCenter(:)==xC);
        rad0 = radRg(radBest(xI));
        [zC,zI] = max(rad0);
        xI = xI(zI);
        % get the pixels in the bead and update available map -----
        [mx,my] = ind2sub(size(datIn),xI);
        rad0 = radRg(radBest(xI));
        dist0 = sqrt((X-mx).^2 + (Y-my).^2);
        datInAvoid(dist0 <= (rad0+1-0.01)) = 1;
        exy = find(dist0<=(rad0-0.01));
        resBead{beadCnt} = exy;
        beadCnt = beadCnt + 1;
    else
        break
    end
end

end







