function resBead = beadMatchPostBig( dat, datMask, ofstx, ofsty, radRg, plotme )

thrx = 0.5;

[Nx,Ny] = size(dat);
d1 = sqrt(dat);
% d1 = dat;
d1bw = 1*(d1>0);
d1bw = bwareaopen(d1bw,20);

% extract background -----
d1CC = bwconncomp(d1bw,4);
nCC = d1CC.NumObjects;
pixCC = d1CC.PixelIdxList;
for ii=1:nCC
    pixMe = pixCC{ii};
    if sum(datMask(pixMe))>10
        break
    end
end
d1Mask = d1*0;
d1Mask(pixMe) = 1;
d1 = d1.*d1Mask;

% template matching -----
Gmag = imgradient(d1);
Gmag(d1bw==0) = 0;  % !! or -1?
matchMap = zeros(Nx,Ny,length(radRg));

for ii=1:length(radRg)
    myRef = zeros(21,21);
    [ax,ay] = find(myRef>-100);
    C0 = radRg(ii);
    aDist = sqrt((ax-11).^2 + (ay-11).^2);
    aOut = aDist>C0-2 & aDist<=C0+1;
    aOut0 = aDist>C0-1 & aDist<=C0;
    aOut1 = aDist>C0-2 & aDist<=C0-1;
    aOut2 = aDist>C0 & aDist<=C0+1;
    aIn = aDist<C0-2;
    myRef(aOut0) = 1;
    myRef(aOut1) = 1/sqrt(2);
    myRef(aOut2) = 1/sqrt(2);
    myRef(aOut) = myRef(aOut)/sum(myRef(aOut));
    myRef(aIn) = -1;  % -1?
    myRef(aIn) = myRef(aIn)/abs(sum(myRef(aIn)));
%     myRef(aIn) = 0;
    
    Cs = conv2(Gmag,myRef,'same');
    Cs(Cs<0) = 0;
    Cs = Cs.*d1bw;
    %     Csf = flipud(Cs);
    matchMap(:,:,ii) = Cs;
    %     figure;contour(flipud(Cs),'ShowText','on');title(num2str(ii));
end

radMin = 9;
dist0 = bwdist(1-d1bw);
dist1 = bwdist(1-datMask);
[matchMapBest,radBest] = max(matchMap,[],3);
matchCenter = (imregionalmax(matchMapBest)>0 & matchMapBest>thrx & dist0>=radMin & (dist1>=radMin | matchMapBest>1)  & datMask>0) .* matchMapBest;
myScore = matchCenter;
% myScore = matchCenter + radBest/5 + dist0/5;

if plotme
    figure;mesh(flipud(Gmag));
    figure;contour(flipud(matchMapBest),'ShowText','on');title(num2str(ii));
    K0 = cat(3,matchCenter,dat,dat*0);
    figure;imshow(K0);
end

[X,Y] = meshgrid(1:Nx,1:Ny);
X = X.';
Y = Y.';
resBead = {};
beadCnt = 1;
while max(myScore(:))>thrx
    [~,xI] = max(myScore(:));
    [mx,my] = ind2sub(size(d1),xI);
    rad0 = radRg(radBest(xI));
    dist0 = sqrt((X-mx).^2 + (Y-my).^2);    
    if min(dat(dist0 <= (rad0/2)))>0.1 && datMask(mx,my)>0
        xx = dat(dist0<=rad0);
        vm = var(xx);
        if vm<0.05
            [ex,ey] = find(dist0<=rad0);
            resBead{beadCnt} = [ex+ofstx-1,ey+ofsty-1];
            beadCnt = beadCnt + 1;
            myScore(dist0<=(rad0)) = 0;
        else
            myScore(dist0<=2) = 0;
        end
    else
        myScore(dist0<=2) = 0;
    end
end

end





