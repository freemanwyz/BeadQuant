%% Bead signal from grace
p = getInitParamBead();
nFrame = 32;

fin = [p.tp filesep 'dump' filesep 'beadAvg16.tif'];
datAvg = double(imread(fin))/65535;
[Nx,Ny] = size(datAvg);

%% data
fin = [p.tp filesep 'dat' filesep '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128' filesep ...
    '061815 BCP 50nM 488 3.5%750_561 1% 800 0-128.lsm - Ch1-T1 - C1 Z1 T'];
dat = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
    fin0 = [fin num2str(ii) '.tif'];
    dat(:,:,ii) = double(imread(fin0))/65535;
end
dat = dat(500:end,500:end,:);
dats = reshape(dat,[],nFrame);
datRef = datAvg(500:end,500:end);

%% Variance and mean by time series
dats1 = dats;
datMean = mean(dats1,2);
meanRg = 0:0.01:1;
N1 = length(meanRg)-1;
% xMean = zeros(1,N1);
xVar = zeros(1,N1);
xVarVar = zeros(1,N1);
xVarSd = zeros(1,N1);
nPix = zeros(1,N1);
% for ii=11
for ii=1:N1
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
figure;errorbar(meanRg1,xVar,xVarSd);
figure;plot(meanRg1,xVar);
figure;semilogy(meanRg1,nPix);

%% Variance and mean by single image
dat0 = dat(:,:,1);
dat0 = randn(size(dat0))/5+0.5;

meanRg = 0:0.01:1;
N1 = length(meanRg)-1;
xVar1 = zeros(1,N1);
nPix1 = zeros(1,N1);

for ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(dat0>=x1 & dat0<x2);
    nPix1(ii) = length(idx);
    [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'median',1 );    
    xVar1(ii) = varOut;
end

%% Gaussian noise
dat0 = dat(:,:,1);
dat0 = randn(size(dat0))/5+0.5;
meanRg = 0:0.01:1;
N1 = length(meanRg)-1;
xVar2 = zeros(1,N1);
nPix2 = zeros(1,N1);

for ii=1:N1
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

%% plot all
meanRg1 = meanRg(1:(end-1));
% figure;plot(meanRg1,xVar1);
% figure;semilogy(meanRg1,nPix1);

figure;
plot(meanRg1,xVar,'r');ylim([0,0.02]);
hold on
plot(meanRg1,xVar1,'b');ylim([0,0.02]);
% plot(meanRg1,xVar2,'k');ylim([0,0.02]);
hold off
legend('Multiple','Single','Gaussian');













