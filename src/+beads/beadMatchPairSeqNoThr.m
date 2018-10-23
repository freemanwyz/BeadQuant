function resBead = beadMatchPairSeqNoThr( datIn, myRef, plotme )
% BEADMATCHPAIRSEQ Compute border pair score and sequential removal

thrx = 0.05;
[Nx0,Ny0] = size(datIn);
resBead = {};
beadCnt = 1;
datInAvoid = zeros(Nx0,Ny0);
nTry = 1;

while 1
    fprintf('Ntry: %d\n',nTry);
    nTry = nTry + 1;
    datEnlarge = zeros(Nx0+40,Ny0+40);
    datEnlarge(21:(Nx0+20),21:(Ny0+20)) = datIn;
    dat = datEnlarge;
    
    datAvoidEnlarge = zeros(Nx0+40,Ny0+40)-1;
    datAvoidEnlarge(21:(Nx0+20),21:(Ny0+20)) = datInAvoid;
    
    % scan radius -----
    radRg = [myRef.radius];
    [Nx,Ny] = size(dat);
    [canX,canY] = find(datAvoidEnlarge==0);
    nCan = length(canX);
    matchMap = zeros(Nx,Ny,length(radRg));
    
    parfor ii=1:length(radRg)
        ref0 = myRef(ii).mask;
        Cs = zeros(Nx,Ny);
        for jj=1:nCan
            canX0 = canX(jj);
            canY0 = canY(jj);
            rgX0 = (canX0-20):(canX0+20);
            rgY0 = (canY0-20):(canY0+20);
            imgCrop = dat(rgX0,rgY0);
            dif0 = imgCrop(ref0(:,2)) - imgCrop(ref0(:,3));
            avoidCrop = datAvoidEnlarge(rgX0,rgY0);
            idxAvoid = avoidCrop(ref0(:,1))==1;
            dif0 = dif0(~idxAvoid);
            nAvoid = sum(idxAvoid);
            if nAvoid < size(ref0,1)*0.3
                Cs(canX0,canY0) = median(dif0);
            else
                Cs(canX0,canY0) = 0;
            end
            
%             if ii==1 && canX0==122 && canY0==42
%                 tmp1 = imgCrop*0;
%                 xx = (imgCrop(ref0(:,1)) + imgCrop(ref0(:,2)))/2 - imgCrop(ref0(:,3));
%                 xx(xx>0) = 1;
%                 xx(xx<0) = 0.25;
%                 xx(xx==0) = 0.5;
%                 tmp1(ref0(:,1)) = xx;
%                 K3 = cat(3,tmp1,imgCrop,imgCrop*0);
%                 figure;imshow(K3)
%                 keyboard
%             end
            
        end
        Cs(Cs<0) = 0;
        matchMap(:,:,ii) = Cs;
        %     figure;contour(flipud(Cs),'ShowText','on');title(num2str(ii));
    end
    
    % find best one -----
    [matchMapBest,radBest] = max(matchMap,[],3);
    matchMapBest = matchMapBest(21:(Nx0+20),21:(Ny0+20));
    radBest = radBest(21:(Nx0+20),21:(Ny0+20));
    matchCenter = (imregionalmax(matchMapBest)>0 & matchMapBest>thrx) .* matchMapBest;
    [X,Y] = meshgrid(1:Nx0,1:Ny0);
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
        [xC,xI] = max(matchCenter(:));        
        xI = find(matchCenter(:)==xC);
        rad0 = radRg(radBest(xI));
        [zC,zI] = max(rad0);
        xI = xI(zI);
        [mx,my] = ind2sub(size(datIn),xI);
        rad0 = radRg(radBest(xI));
        dist0 = sqrt((X-mx).^2 + (Y-my).^2);
        [ex,ey] = find(dist0<=(rad0-0.01));
        datInAvoid(dist0 <= (rad0+1-0.01)) = 1;
        resBead{beadCnt} = [ex,ey];
        beadCnt = beadCnt + 1;
    else
        break
    end
end

end







