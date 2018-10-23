function [resBead,resBeadCenter,resBeadRad] = beadMatchMulti( datIn, datMask, myRef, ofstx, ofsty, plotme, thrMulti )
% BEADMATCHPAIRSEQ Compute border pair score and sequential multiple removal

datIn = sqrt(datIn);
thrx = 0.1;
[Nx0,Ny0] = size(datIn);
resBead = {};
resBeadCenter = {};
resBeadRad = {};

beadCnt = 1;
datInAvoid = zeros(Nx0,Ny0)-1;
datInAvoid(datMask>0) = 0;
nTry = 1;

while 1
%     fprintf('Ntry: %d\n',nTry);
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
    
    for ii=1:length(radRg)
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
            if nAvoid < size(ref0,1)*0.4
                Cs(canX0,canY0) = median(dif0);
            else
                Cs(canX0,canY0) = 0;
            end
        end
        Cs(Cs<0) = 0;
        matchMap(:,:,ii) = Cs;
    end
    
    % find local maximum -----
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
        xCRef = max(matchCenter(:));
        while 1
            xC = max(matchCenter(:));
            if xC < xCRef*thrMulti  && thrMulti<1
                break
            end
            xI = find(matchCenter(:)==xC);
            rad0 = radRg(radBest(xI));
            [~,zI] = max(rad0);
            xI = xI(zI);
            [mx,my] = ind2sub(size(datIn),xI);                     
            rad0 = radRg(radBest(xI));
            dist0 = sqrt((X-mx).^2 + (Y-my).^2);
            [ex,ey] = find(dist0<=(rad0-0.01));
            datInAvoid(dist0 <= (rad0+1-0.01)) = 1;
            matchCenter(dist0 <= (rad0+1-0.01)) = 0;    % for remove multiple
            resBead{beadCnt} = [ex+ofstx-1,ey+ofsty-1];
            resBeadCenter{beadCnt} = [mx+ofstx-1,my+ofsty-1];
            resBeadRad{beadCnt} = rad0;
            beadCnt = beadCnt + 1;
            if  thrMulti>=1
                break
            end
        end        
    else
        break
    end
end

end







