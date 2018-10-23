function resBead = beadMatchPostSmall( dat, datMask, ofstx, ofsty)

resBead = {};
statSel = regionprops(datMask,'Area','Perimeter');
rtSel = statSel.Perimeter^2/statSel.Area/4/pi;
if rtSel>1.5
    return
end

[Nx,Ny] = size(dat);
d1 = sqrt(dat);

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

[X,Y] = meshgrid(1:Nx,1:Ny);
X = X.';
Y = Y.';
dist1 = bwdist(1-datMask);
[xC,xI] = max(dist1(:));
[mx,my] = ind2sub(size(d1),xI);
dist0 = sqrt((X-mx).^2 + (Y-my).^2);
% notTooMuchOut = sum(dat(dist0 <= (xC))==0)<2;
notTooMuchOut = 1;
if xC>=2 && notTooMuchOut && sum(dist0(:)<=xC)>15
    [ex,ey] = find(dist0<=xC);
    resBead{1} = [ex+ofstx-1,ey+ofsty-1];
end

end





