%% Bead signal from grace
p = getInitParamBead();
p.tp = 'D:\neuro_WORK\glia_ucdavis_WORK\spontanous activity\';
fin = [p.tp filesep 'control.lsm'];
datIn = tiffread(fin);
nFrame = length(datIn);
[Nx,Ny] = size(datIn(1).data);
datAvg = zeros(Nx,Ny);
nMax = 2^datIn(1).bits-1;

for ii=1:nFrame
    datAvg = datAvg + double(datIn(ii).data)/nMax;
end
datAvg = datAvg/nFrame;

% datAvg = datAvg - min(datAvg(:));
% datAvg = datAvg/max(datAvg(:));

dat = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
    fin0 = [fin num2str(ii) '.tif'];
    dat(:,:,ii) = double(datIn(ii).data)/nMax;
end
dats = reshape(dat,[],nFrame);
datRef = datAvg;

%% parameters
figId = 1;
meanRg = 0.05:0.01:1;
N1 = length(meanRg)-1;
meanRg1 = meanRg(1:(end-1));

%% Variance and mean by time series
dats1 = dats;
datMean = mean(dats1,2);
xVar = zeros(1,N1);
xVarVar = zeros(1,N1);
xVarSd = zeros(1,N1);
nPix = zeros(1,N1);
% for ii=11
parfor ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(datMean>=x1 & datMean<x2);
    nPix(ii) = length(idx);
    datSel = dats1(idx,:);
    varSel = var(datSel,0,2);
    varSel = varSel(varSel<0.05);
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
figure;errorbar(meanRg1,xVar,xVarSd);
figure;plot(meanRg1,xVar);
figure;semilogy(meanRg1,nPix);

%% do not use edge
dat0 = dat(:,:,figId);
% dat0 = datRef;
[Gmag,Gdir] = imgradient(dat0);
dat0(Gmag>0.8) = Inf;
K0 = cat(3,Gmag>0.8,datRef,datRef*0);
figure;imshow(K0);

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

%% Variance and mean by single image
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

%% Variance and median by single image
dat0 = dat(:,:,figId);
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

%% plot all
figure;
plot(meanRg1,xVar,'r');ylim([0,0.03]);
hold on
plot(meanRg1,xVar1/2,'b');
plot(meanRg1,xVar1a/2,'g');
plot(meanRg1,xVar3/2,'k');
hold off
ylabel('Variance');
legend('Multiple','Single','Single median','Single without edge','Location','NorthWest');

figure;
plot(meanRg1,log10(nPix+1),'r');
hold on
plot(meanRg1,log10(nPix1+1),'b');
plot(meanRg1,log10(nPix3+1),'k');
hold off
ylabel('log10(Pixel number)');
legend('Multiple','Single or Single median','Single without edge');










