%Bead selection pipeline
% debug only

%% load files
if 1
    warning('off', 'Images:initSize:adjustingMag');
    [fname,stub] = uigetfile('*.lsm','Select the background .lsm file');
    path0 = fullfile(stub,fname);
    [fname,stub] = uigetfile(fullfile(stub,'*.lsm'),'Select the .lsm file you are analyzing');
    path1 = fullfile(stub,fname);
    
    % background
    dat0 = tiffread(path0);
    maxVal = 2^dat0.bits-1;
    bgch1 = double(dat0.data{1})/maxVal;
    bgch2 = double(dat0.data{2})/maxVal;
    % [Nx,Ny] = size(bgch1);
    
    % time lapse
    dat1 = tiffread(path1);
    nFrame = length(dat1);
    % sigch2 = zeros(Nx,Ny,nFrame);
    for ii=1:nFrame
        tmp = dat1(ii);
        datEle = double(tmp.data{1})/maxVal;
        if ii==1
            sigch1Mean = datEle;
        else
            sigch1Mean = sigch1Mean + datEle;
        end
        %     datEle = double(tmp.data{2})/maxVal;
        %     sigch2(:,:,ii) = datEle;
    end
    sigch1Mean = sigch1Mean/nFrame;
end

%%
cropImg = 0;
addNoise = 0;
if cropImg
    imga = sigch1Mean(300:500,300:500);
else
    imga = sigch1Mean;
end
image_intensity0 = 0:0.01:1;
varMax = 0.3;
var0 = 0:varMax/100:varMax;
if addNoise
    imgb = imnoise(imga,'localvar',image_intensity0,var0);
else
    imgb = imga;
end
d0 = sqrt(imgb);

% d0 = wiener2(d0,[7 7]);
% d0 = dd.bgMean;  % Please first load the session file from beadGUI for this option
% d0 = sqrt(sigch1Mean);

% d0a = d0;
% d0a(d0a<0.25) = 0;
% xx = bwareaopen(d0a>0,50,4);

%% bead detection using Hough transform
d0a = d0;
% d0a = imgaussfilt(d0a,2);
% d0a = wiener2(d0a,[5 5]);
% d0a(d0a<0.25) = 0;
% xx = bwareaopen(d0a>0,50,4);
% d0a(xx==0) = 0;
[centersBright, radiiBright, metric, accumMatrix] = imfindcircles0(d0a,[3 14],-1,0.98);
K0 = cat(3,d0*0,imgb,imgb*0);
figure;imshow(K0);
viscircles(centersBright, radiiBright,'EdgeColor','b');

a0 = abs(accumMatrix);
% a0(xx==0) = 0;
K0 = cat(3,a0/max(a0(:)),imga*0.5,imga*0);
figure;imshow(K0);
% figure;contour(flipud(a0)/max(a0(:)))
% figure;mesh(a0)

%% our method
if 1
    d0a = d0;
%     d0a = imgaussfilt(d0a,2);
    d0a(d0a<0.25) = 0;
    xx = bwareaopen(d0a>0,50,4);
    d0a(xx==0) = 0;
    [resBead,resBeadCenter,resBeadRad] = beadDetection(d0a,0.25,0.8,5,10,0);
    K0 = cat(3,d0*0,imgb,imgb*0);
    figure;imshow(K0);
    x0 = cell2mat(resBeadCenter');
    x0 = [x0(:,2),x0(:,1)];
    viscircles(x0, cell2mat(resBeadRad'), 'EdgeColor','b');
end

%% our method, score only
if 1
    myRg = 4:0.5:13;
    myRef = beadTemplate(myRg);
    d0a = d0;
%     d0a(d0a<0.25) = 0;
%     xx = bwareaopen(d0a>0,50,4);
%     d0a(xx==0) = 0;
    datIn = d0a;
    [Nx0,Ny0] = size(datIn);
    
    datInAvoid = zeros(Nx0,Ny0)-1;
    datInAvoid(datIn>0.1) = 0;
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
    matchMapPhase = zeros(Nx,Ny,length(radRg));
    
    rMin = min(myRg);
    rMax = max(myRg);
    rPhi = exp(1i*2*pi*(log(myRg)-log(rMin))/(log(rMax)-log(rMin)))-pi;

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
%             Cs(canX0,canY0) = mean(dif0);
            Cs(canX0,canY0) = median(dif0);
        end
        Cs(Cs<0) = 0;
        matchMap(:,:,ii) = Cs;
        Cs = Cs*rPhi(ii);
        matchMapPhase(:,:,ii) = Cs;
    end
    
    matchMapPhaseBest = abs(sum(matchMapPhase,3));
    matchMapPhaseBest = matchMapPhaseBest(21:(Nx0+20),21:(Ny0+20));
    [matchMapBest,radBest] = max(matchMap,[],3);
    matchMapBest = matchMapBest(21:(Nx0+20),21:(Ny0+20));
    
    K1 = cat(3,matchMapBest/max(matchMapBest(:)),imga*0.5,imga*0);
    figure;imshow(K1);title('ours, max');
    K1p = cat(3,matchMapPhaseBest/max(matchMapPhaseBest(:)),imga*0.5,imga*0);
    figure;imshow(K1p);title('ours, phase');
    
    figure;contour(flipud(matchMapPhaseBest/max(matchMapPhaseBest(:))))
    figure;contour(flipud(matchMapBest/max(matchMapBest(:))))
    figure;mesh(matchMapPhaseBest)
    
end

%%
if 0
    x0 = a0/max(a0(:));
    x1 = matchMapBest/max(matchMapBest(:));
    nr = 134;
    figure;
    plot(x0(nr,:))
    hold on
    plot(x1(nr,:))
    legend('Hough','Ours');
end

%% remove wrong circles
if 0
    % start from low score beads
    D = squareform(pdist(centersBright));
    for ii=length(metric):-1:1
        if mod(ii,1000)==0
            fprintf('%d\n',ii);
        end
        idx1 = D(ii,:) < radiiBright(ii);
        idx2 = D(ii,:) < radiiBright(ii)/2;
        idx1(ii) = 0;
        idx2(ii) = 0;
        if sum(idx1)>0 && radiiBright(ii)>8
            if metric(ii) < max(metric(idx))
                %             D(ii,:) = 0;
                D(:,ii) = 0;
            end
        elseif sum(idx2)>0 && radiiBright(ii)<=8
            if metric(ii) < max(metric(idx2))
                D(:,ii) = 0;
            end
        end
    end
    
    idx = mean(D,1)>0;
    c1 = centersBright(idx,:);
    r1 = radiiBright(idx);
    figure;imshow(d0);
    viscircles(c1, r1,'EdgeColor','b');
end

%%
if 0
    resBorder = zeros(Nx,Ny);
    resCenter = zeros(Nx,Ny);
    myRg = 4:15;
    myRef = beadTemplate(myRg);
    for ii=1:size(centersBright,1)
        if mod(ii,1000)==0
            fprintf('ii is %d\n',ii);
        end
        idx = round(centersBright(ii,:));
        resCenter(idx(2),idx(1)) = 1;
        rad0 = round(radiiBright(ii))+1;
        if rad0<4
            rad0 = 4;
        end
        if rad0>15
            rad0 = 15;
        end
        nx = idx(2) + myRef(rad0-3).maskX(:,1) - 21;
        ny = idx(1) + myRef(rad0-3).maskY(:,1) - 21;
        nxy = sub2ind(size(sigch1Mean),nx,ny);
        resBorder(nxy) = 1;
    end
    K2 = cat(3,resBorder,sigch1Mean,resCenter);
    figure;imshow(K2)
    imwrite(K2,'061815_Hough.tiff');
end


