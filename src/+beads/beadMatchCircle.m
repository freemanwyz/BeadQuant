function resBead = beadMatchCircle( dat, ofstx, ofsty, radRg, plotme )

thrx = 0.15;
% thrx = 0;

[Nx,Ny] = size(dat);
d1 = sqrt(dat);
% d1 = dat;
d1bw = 1*(d1>0);
d1bw = bwareaopen(d1bw,20);
matchMap = zeros(Nx,Ny,length(radRg));

for ii=1:length(radRg)
    myRef = zeros(41,41);
    [ax,ay] = find(myRef>-100);
    C0 = radRg(ii);
    aDist = sqrt((ax-21).^2 + (ay-21).^2);
    aOut0 = aDist>C0-2 & aDist<=C0;
%     aOut0 = aDist>C0-1 & aDist<=C0;
%     aOut0 = aDist<=C0;
    aOut2 = aDist>C0 & aDist<=C0+1;
    myRef(aOut0) = 1/sum(aOut0);
    myRef(aOut2) = -1/sum(aOut2);    
    Cs = conv2(d1,myRef,'same');
    Cs(Cs<0) = 0;
    Cs = Cs.*d1bw;
    matchMap(:,:,ii) = Cs;
%     figure;contour(flipud(Cs),'ShowText','on');title(num2str(ii));
end

[matchMapBest,radBest] = max(matchMap,[],3);
matchCenter = (imregionalmax(matchMapBest)>0 & matchMapBest>thrx) .* matchMapBest;
[X,Y] = meshgrid(1:Nx,1:Ny);
X = X.';
Y = Y.';

if plotme
    figure;contour(flipud(matchMapBest),'ShowText','on');title(num2str(ii));
    K0 = cat(3,matchCenter>0,dat,dat*0);
    figure;imshow(K0);
end

resBead = {};
beadCnt = 1;
while max(matchCenter(:))>thrx
    [~,xI] = max(matchCenter(:));
    [mx,my] = ind2sub(size(d1),xI);
    rad0 = radRg(radBest(xI));
    dist0 = sqrt((X-mx).^2 + (Y-my).^2);
    matchCenter(dist0<=(rad0)) = 0;
    if min(dat(dist0 <= (rad0/2)))>0.1
        [ex,ey] = find(dist0<=rad0);
        resBead{beadCnt} = [ex+ofstx-1,ey+ofsty-1];
        beadCnt = beadCnt + 1;
    end
end


end





