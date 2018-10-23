% %% truncated normal
% NN = 1e6;
% xx0 = randn(NN,1);
% rtxRg = 0.1:0.1:2;
% xxVar = zeros(1,length(rtxRg));
% xxSelVar = xxVar;
% for ii=1:length(rtxRg)
%     rtx = rtxRg(ii);
%     xx = xx0*rtx;
%     [xC,xI] = sort(abs(xx),'descend');
%     xxSel = xx(xI((round(NN*0.15)+1):NN));
%     xxVar(ii) = var(xx);
%     xxSelVar(ii) = var(xxSel);
% end
% mean(xxVar./xxSelVar)

%% Bead signal from grace
Nx = 1000;
Ny = 1000;
datAvg = randi(255,Nx,Ny);
nFrame = 20;
rt0 = 25;

dat = zeros(Nx,Ny,nFrame);
parfor ii=1:nFrame
    datNoisy0 = datAvg + randn(Nx,Ny).*sqrt(datAvg*rt0);
    dat(:,:,ii) = datNoisy0;
end
dats = reshape(dat,[],nFrame);
datRef = datAvg;

% parameters -----
figId = 10;
meanRg = 0:255;
N1 = length(meanRg)-1;
meanRg1 = meanRg(1:(end-1));

%% Variance and mean by time series
dats1 = dats;
datMean = mean(dats1,2);
xVar = zeros(1,N1);
xVarVar = zeros(1,N1);
xVarSd = zeros(1,N1);
nPix = zeros(1,N1);
parfor ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(datMean>=x1 & datMean<x2);
    nPix(ii) = length(idx);
    datSel = dats1(idx,:);
    varSel = var(datSel,0,2);
    xVar(ii) = mean(varSel);
    xVarVar(ii) = var(varSel);
    xVarSd(ii) = std(varSel);    
end
figure;errorbar(meanRg1,xVar,xVarSd);
figure;plot(meanRg1,xVar);
figure;semilogy(meanRg1,nPix);

%% Variance and mean by single image
dat0 = dat(:,:,figId);
xVar1 = zeros(1,N1);
nPix1 = zeros(1,N1);
for ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(dat0>=x1 & dat0<x2);
    trueVar = datAvg(idx)*rt0;
    xVar1(ii) = mean(trueVar);
    nPix1(ii) = length(idx);
%     [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
%     xVar1(ii) = varOut;
end

%%
figure;
plot(meanRg1,xVar,'b');
hold on
plot(meanRg1,xVar1,'r');
hold off
legend('Average','Single','Location','NorthWest')








