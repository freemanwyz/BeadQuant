%% Bead signal from grace
p = getInitParamBead();
figId = 1;
meanRg = 0:0.01:1;
N1 = length(meanRg)-1;
p.crop = 0;

%% data
% path1 = [p.tp filesep 'dat' filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128.lsm'];
[fname,stub] = uigetfile('*.lsm','Select the .lsm file you are analyzing');
path1 = fullfile(stub,fname);

dat1 = tiffread(path1);
nFrame = length(dat1);
maxVal = 2^dat1(1).bits-1;
if iscell(dat1(1).data)
    [Nx,Ny] = size(dat1(1).data{1});
else
    [Nx,Ny] = size(dat1(1).data);
end

sigch1 = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
    tmp = dat1(ii);
    if iscell(tmp.data)
        datEle = double(tmp.data{1})/maxVal;
    else
        datEle = double(tmp.data)/maxVal;
    end
    if ii==1
        sigch1Mean = datEle;
    else
        sigch1Mean = sigch1Mean + datEle;
    end
    sigch1(:,:,ii) = datEle;
end
sigch1Mean = sigch1Mean/nFrame;

if p.crop
    dat = sigch1(500:end,500:end,:);
    datRef = sigch1Mean(500:end,500:end);
else
    dat = sigch1;
    datRef = sigch1Mean;
end
dats = reshape(dat,[],nFrame);

%% foi's method
if 0
    dat0 = dat(:,:,figId);
    get_var_foi(dat0);
end

%% watershed
if 0
    fprintf('Watershed dist ----- \n');
    dat0 = dat(:,:,figId);
    
    bw = wiener2(dat0,[7,7])>0.05;
    D = bwdist(~bw);
    D = -D;
    D(~bw) = -Inf;
    % D = -dat0;
    
    L = watershed(D);
    
    tmp = 1*(L==0);
    tmp = imdilate(tmp,strel('square',4));
    L(tmp>0) = 0;
    
    rgb0 = label2rgb(L,'jet',[.5 .5 .5]);
    figure;imshow(rgb0,'InitialMagnification','fit')
    
    dat0(L==0) = Inf;
    xVar4 = zeros(1,N1);
    nPix4 = zeros(1,N1);
    
    parfor ii=1:N1
        fprintf('%d\n',ii);
        x1 = meanRg(ii);
        x2 = meanRg(ii+1);
        idx = find(dat0>=x1 & dat0<x2);
        nPix4(ii) = length(idx);
        [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
        xVar4(ii) = varOut;
    end
end

%% do not use edge
if 1
    fprintf('No edge ----- \n');
    dat0 = dat(:,:,figId);
    % dat0 = datRef;
    [Gmag,Gdir] = imgradient(dat0);
    dat0(Gmag>0.5) = Inf;
    xVar3 = zeros(1,N1);
    nPix3 = zeros(1,N1);
    
    parfor ii=1:N1
        fprintf('%d\n',ii);
        x1 = meanRg(ii);
        x2 = meanRg(ii+1);
        idx = find(dat0>=x1 & dat0<x2);
        nPix3(ii) = length(idx);
        [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
        xVar3(ii) = varOut;
    end
end

%% Variance and mean by time series
fprintf('Multi ----- \n');
dats1 = dats;
datMean = mean(dats1,2);
% xMean = zeros(1,N1);
xVar = zeros(1,N1);
xVarVar = zeros(1,N1);
xVarSd = zeros(1,N1);
nPix = zeros(1,N1);
% for ii=11
for ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(datMean>=x1 & datMean<x2);
    nPix(ii) = length(idx);
    datSel = dats1(idx,:);
    varSel = var(datSel,0,2);
    varSel = varSel(varSel<0.02);
    xVar(ii) = mean(varSel);
    xVarVar(ii) = var(varSel);
    xVarSd(ii) = std(varSel);
    
    %     datB = datRef*0;
    % %     idx1 = idx(varSel>0.02);
    %     idx1 = idx;
    %     datB(idx1) = 1;
    %     K0 = cat(3,datB,datRef.*(1-datB)*0.3,datRef*0);
    %     figure;imshow(K0);
    
end

meanRg1 = meanRg(1:(end-1));
% figure;errorbar(meanRg1,xVar,xVarSd);
% figure;plot(meanRg1,xVar);
% figure;semilogy(meanRg1,nPix);

%% Variance and mean by single image
if 1
    fprintf('Mean ----- \n');
    dat0 = dat(:,:,figId);
    % dat0 = datRef;
    xVar1 = zeros(1,N1);
    nPix1 = zeros(1,N1);
    
    parfor ii=1:N1
        fprintf('%d\n',ii);
        x1 = meanRg(ii);
        x2 = meanRg(ii+1);
        idx = find(dat0>=x1 & dat0<x2);
        nPix1(ii) = length(idx);
        [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
        xVar1(ii) = varOut;
    end
end

%% Variance and median by single image
if 1
    fprintf('Median ----- \n');
    dat0 = dat(:,:,figId);
    % dat0 = datRef;
    xVar1a = zeros(1,N1);
    nPix1a = zeros(1,N1);
    
    parfor ii=1:N1
        fprintf('%d\n',ii);
        x1 = meanRg(ii);
        x2 = meanRg(ii+1);
        idx = find(dat0>=x1 & dat0<x2);
        nPix1a(ii) = length(idx);
        [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'median',1 );
        xVar1a(ii) = varOut;
    end
end

%% Gaussian noise
if 0
    dat0 = dat(:,:,figId);
    dat0 = randn(size(dat0))/5+0.5;
    meanRg = 0:0.01:1;
    N1 = length(meanRg)-1;
    xVar2 = zeros(1,N1);
    nPix2 = zeros(1,N1);
    
    parfor ii=1:N1
        fprintf('%d\n',ii);
        x1 = meanRg(ii);
        x2 = meanRg(ii+1);
        idx = find(dat0>=x1 & dat0<x2);
        nPix2(ii) = length(idx);
        [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'median',1 );
        xVar2(ii) = varOut;
    end
    
    idx = 1:length(dat0(:));
    [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',2 );
end

%% plot all
xMax = max(xVar);
% xMax = max(sqrt(xVar));
if 1
    meanRg1 = meanRg(1:(end-1));
    % figure;plot(meanRg1,xVar1);
    % figure;semilogy(meanRg1,nPix1);
    
    figure;
    plot(meanRg1,xVar,'r');ylim([0,xMax*1.5]);
    hold on
    plot(meanRg1,xVar1/2,'b');
    %     plot(meanRg1,xVar4/2,'g');
    plot(meanRg1,xVar3/2,'k');
    % plot(meanRg1,xVar2,'k');ylim([0,0.02]);
    hold off
    ylabel('Variance');
    legend('Multiple','Single','Single without edge','Location','NorthWest');
    
    figure;
    plot(meanRg1,log10(nPix+1),'r');
    hold on
    plot(meanRg1,log10(nPix1+1),'b');
    %     plot(meanRg1,log10(nPix4+1),'g');
    plot(meanRg1,log10(nPix3+1),'k');
    hold off
    ylabel('Pixel number');
    legend('Multiple','Single','Single without edge');
end

xMax = max(sqrt(xVar));
if 0
    meanRg1 = meanRg(1:(end-1));
    % figure;plot(meanRg1,xVar1);
    % figure;semilogy(meanRg1,nPix1);
    
    figure;
    plot(meanRg1,sqrt(xVar),'r');ylim([0,xMax*1.5]);
    hold on
    plot(meanRg1,sqrt(xVar1/2),'b');
    %     plot(meanRg1,xVar4/2,'g');
    plot(meanRg1,sqrt(xVar3/2),'k');
    % plot(meanRg1,xVar2,'k');ylim([0,0.02]);
    hold off
    ylabel('sigma');
    legend('Multiple','Single','Single without edge','Location','NorthWest');
end

% figure;plot(meanRg1,sqrt(xVar),'r');
% hold on;plot(meanRg1,sqrt(xVar3),'b');hold off;
% legend('Multiple','Single');






