function resBead = beadMatchPair( datIn, myRef, ofstx, ofsty, plotme )
% BEADMATCHPAIR Compute border pair score
% not tested !!!

% thrx = 0.15;
thrx = 0;

[Nx0,Ny0] = size(datIn);
datEnlarge = zeros(Nx0+40,Ny0+40);
datEnlarge(21:(Nx0+20),21:(Ny0+20)) = datIn;
dat = datEnlarge;

% scan radius -----
radRg = [myRef.radius];
[Nx,Ny] = size(dat);
d1 = sqrt(dat);
d1bw = bwareaopen(d1>0,20);
[canX,canY] = find(d1bw>0);
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
        imgCrop = d1(rgX0,rgY0);
        Cs(canX0,canY0) = median( (imgCrop(ref0(:,1)) + imgCrop(ref0(:,2)))/2 - imgCrop(ref0(:,3)));
        
        %         if ii==9 && canX0==75 && canY0==51
        %             tmp1 = imgCrop*0;
        %             xx = (imgCrop(ref0(:,1)) + imgCrop(ref0(:,2)))/2 - imgCrop(ref0(:,3));
        %             xx(xx>0) = 1;
        %             xx(xx<=0) = 0.5;
        %             tmp1(ref0(:,1)) = xx;
        %             K3 = cat(3,tmp1,imgCrop,imgCrop*0);
        %             figure;imshow(K3)
        %             keyboard
        %         end
        
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
[X,Y] = meshgrid(1:Nx,1:Ny);
X = X.';
Y = Y.';

% dump -----
if plotme
    figure;contour(flipud(matchMapBest),'ShowText','on');title(num2str(ii));
    K0 = cat(3,matchCenter>0,datIn,datIn*0);
    figure;imshow(K0);
    K1 = cat(3,radBest/10.*(matchCenter>0),datIn,datIn*0);
    figure;imshow(K1);
end

resBead = {};
beadCnt = 1;
while max(matchCenter(:))>thrx
    [~,xI] = max(matchCenter(:));
    [mx,my] = ind2sub(size(d1),xI);
    rad0 = radRg(radBest(xI));
    dist0 = sqrt((X-mx).^2 + (Y-my).^2);
    matchCenter(dist0<=(rad0)) = 0;
    if min(datIn(dist0 <= (rad0/2)))>0.1
        [ex,ey] = find(dist0<=rad0);
        resBead{beadCnt} = [ex+ofstx-1,ey+ofsty-1];
        beadCnt = beadCnt + 1;
    end
end


end





