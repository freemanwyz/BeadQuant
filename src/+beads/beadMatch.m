function resBead = beadMatch( dat, varargin )

if nargin<3
    ofstx = 0;
    ofsty = 0;
else
    ofstx = varargin{1};
    ofsty = varargin{2};
end
if nargin<4
    radRg = 5:9;
else
    radRg = varargin{3};
end
if nargin<5
    plotme = 0;
else
    plotme = varargin{4};
end

thrx = 0.4;

[Nx,Ny] = size(dat);
d1 = sqrt(dat);
% d1 = dat;
d1bw = 1*(d1>0);
d1bw = bwareaopen(d1bw,20);
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

[matchMapBest,radBest] = max(matchMap,[],3);
matchCenter = (imregionalmax(matchMapBest)>0 & matchMapBest>thrx) .* matchMapBest;
[X,Y] = meshgrid(1:Nx,1:Ny);
X = X.';
Y = Y.';

if plotme
    figure;mesh(flipud(Gmag));
    figure;contour(flipud(matchMapBest),'ShowText','on');title(num2str(ii));
    K0 = cat(3,matchCenter,dat,dat*0);
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





