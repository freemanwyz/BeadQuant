function myCurves = beadExtract(sigch2, resBeadCenter, resBeadRad, rad0, rad1, outPath, maxVal, plotme)
% sigch2 is the (combined) signal channel

[Nx,Ny] = size(sigch2(:,:,1));
nFrame = size(sigch2,3);

% circle template -----
if rad0>=4
    rad0 = rad0 - 1;
end
myRg = rad0:0.5:(rad1+2);
Nx0 = 41;
myTemp = zeros(Nx0,Nx0,length(myRg));
myTempXY = cell(length(myRg),1);
ofstCenter = (Nx0-1)/2+1;
for ii = 1:length(myRg)
    tmp = zeros(Nx0,Nx0);
    [ax,ay] = find(tmp>-100);
    C0 = myRg(ii)-0.01;  % to make the circle look better
    aDist = sqrt((ax-ofstCenter).^2 + (ay-ofstCenter).^2);
    aOut0 = aDist>C0-1 & aDist<=C0;
    tmp(aOut0) = 1/sum(aOut0(:)>0);  % normalize the circle
    myTemp(:,:,ii) = tmp;
    [x,y] = find(tmp>0);
    myTempXY{ii} = [x,y]-ofstCenter;
end

% mean for time lapse channel 2 -----
sigch2Mean = zeros(Nx,Ny);
for ii=1:nFrame
    sigch2Mean = sigch2Mean + double(sigch2(:,:,ii))/maxVal;
end
sigch2Mean = sigch2Mean/nFrame;

% convolution with different circle size -----
bestSizeIdx = ones(Nx,Ny);
sigch2CircleBest = zeros(Nx,Ny);
for ii=1:length(myRg)
    sigch2CircleEle = conv2(sigch2Mean,myTemp(:,:,ii),'same');
    sIdx = sigch2CircleEle > sigch2CircleBest;
    bestSizeIdx(sIdx) = ii;
    sigch2CircleBest(sIdx) = sigch2CircleEle(sIdx);
%     sigch2Circle(:,:,ii) = conv2(sigch2Mean,myTemp(:,:,ii),'same');
end
sIdx = bestSizeIdx;
% [~,sIdx] = max(sigch2Circle,[],3);

% find curves for each bead -----
t0 = zeros(Nx,Ny);
sigs = reshape(sigch2,[],nFrame);
myCurves = zeros(size(resBeadCenter,1),nFrame);
for ii=1:size(resBeadCenter,1)
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    m0 = resBeadCenter(ii,:);  % detected center
    r0 = resBeadRad(ii);  % detected radius
    r1Idx = sIdx(m0(1),m0(2));
    r1 = myRg(r1Idx);
    if r1>r0+1
        r1 = r0;
    end
    % shift the circle templates to the bead center
    r1Idx = find(myRg==r1);
    xy0 = myTempXY{r1Idx};
    xy1 = [m0(1)+xy0(:,1),m0(2)+xy0(:,2)];
    rgIn = 1:1;
    for jj=rgIn  % use two inside circles
        if r1Idx<=jj
            break
        end
        xy0 = myTempXY{r1Idx-jj};
        xy1a = [m0(1)+xy0(:,1),m0(2)+xy0(:,2)];
        xy1 = [xy1;xy1a];
    end
    rgIn = 1:1;
    for jj=rgIn  % use one outside circles
        xy0 = myTempXY{r1Idx+jj};
        xy1a = [m0(1)+xy0(:,1),m0(2)+xy0(:,2)];
        xy1 = [xy1;xy1a];
    end
    % remove out of image ones -----
    idxy = xy1(:,1)>0 & xy1(:,1)<=Nx & xy1(:,2)>0 & xy1(:,2)<=Ny;
    xy1 = xy1(idxy,:);
    xy2 = sub2ind([Nx,Ny],xy1(:,1),xy1(:,2));    
    t0(xy2) = 1;
    % extract curves -----
    sig0 = double(sigs(xy2,:))/maxVal;
    myCurves(ii,:) = mean(sig0,1);
end

if plotme
    K0 = cat(3,t0,sigch2Mean*0,sigch2Mean);
    imwrite(double(K0),[outPath,filesep,'extracted_ring.tif'])
end


end








