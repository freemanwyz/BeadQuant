%% simulate small patches
Nx = 1024;
Ny = 1024;
datAvg = zeros(Nx,Ny);
idx0 = 1:256;
[idx0x,idx0y] = ind2sub([16 16],idx0);
nn = 0;
gap0 = Nx/16;
for ii=1:256
    id1 = ((idx0x(ii)-1)*gap0+1):(idx0x(ii)*gap0);
    id2 = ((idx0y(ii)-1)*gap0+1):(idx0y(ii)*gap0);
    if ii>=86 && ii<=116
%         datAvg(id1,id2) = 0;
        datAvg(id1,id2) = nn;
    else
        datAvg(id1,id2) = nn;
    end    
    nn = nn + 1;
end

nFrame = 20;
rt0 = 1;
dat = zeros(Nx,Ny,nFrame);
parfor ii=1:nFrame
    datNoisy0 = datAvg + randn(Nx,Ny).*sqrt(datAvg*rt0);
%     datNoisy0(datNoisy0<0) = 0;
%     datNoisy0(datNoisy0>255) = 255;
    dat(:,:,ii) = datNoisy0;
end
dats = reshape(dat,[],nFrame);
datRef = datAvg;

%% estimation
figId = 10;
meanRg = 0:255;
meanRg1 = meanRg(1:(end-1));

% Variance and mean by time series
[xVar1,xVar1Var,xVar1Sd,nPix1] = getVarCurveMulti(dats,meanRg);

% Variance and mean by single image
dat0 = dat(:,:,figId);
[xVar2a,nPix2a,nPairs2a] = getVarCurve(dat0,meanRg,1,'mean',2);
[xVar2b,nPix2b,nPairs2b] = getVarCurve(dat0,meanRg,1,'meanAvg',2);
[xVar2c,nPix2c,nPairs2c] = getVarCurve(dat0,meanRg,1,'meanCorrect',2);
[xVar2d,nPix2d,nPairs2d] = getVarCurve(dat0,meanRg,1,'meanSingleDirection',2);

%% plot image
tmp = dat0*0;
xx = 100;
xx1 = xx + 1;
tmp(dat0>=xx & dat0<(xx1)) = 1;
idx000 = find(datAvg==100);
tmp1 = dat0/255.*(1-tmp);
tmp2 = dat0*0;
tmp2(idx000) = 1;
K0 = cat(3,tmp,tmp1,tmp2);
figure;imshow(K0);title(num2str(xx));

idx00 = find(dat0>=xx & dat0<xx1);
x1 = datAvg(idx00);
figure;hist(x1);

%% analysis 1
ii0 = 100;
x1 = meanRg(ii0);
x2 = meanRg(ii0+1);
maskIn = find(dat0>=x1 & dat0<x2);
trueVar = datAvg(maskIn);
[Nx,Ny] = size(dat0);
d = zeros(Nx+2,Ny+2)+Inf;
xrg = 2:(Nx+1);
yrg = 2:(Ny+1);
d(xrg,yrg) = dat0;
zx = zeros(length(maskIn),8);
xvec = [0,1,1,-1,0,-1,1,-1];
yvec = [1,0,1,1,-1,0,-1,-1];
for ii=1:8
    tmp0 = dat0 - d(xrg+xvec(ii),yrg+yvec(ii));
    zx(:,ii) = tmp0(maskIn);
end
zx(isinf(zx) | isnan(zx)) = NaN;
noise = nanmean(zx,2);
noise = noise(~isinf(noise) & ~isnan(noise));

% noise2 = noise.^2*8/9;
noise2 = zx(:,2).^2/2;

x00 = min(trueVar);
x01 = max(trueVar);
varTrueVec = zeros(1,x01-x00+1);
tt = 1;
for jj=x00:x01
    varTrueVec(tt) = nanmean(noise2(trueVar==jj));
    tt = tt + 1;
end

figure;scatter(trueVar,noise2,10,'filled');
hold on
plot(x00:x01,varTrueVec,'r');
plot(x00:x01,varTrueVec*0+ii0,'g');
hold off;
xlabel('True intensity');
ylabel('Estimated variance');
legend('Pair','Variance','Desired variance');

% sn = sum(2*trueVar.^2);
% lya = sum((abs(noise2-trueVar)).^3)/sn^3;

%% analysis 2
kt = 1;
pdf0 = @(murg,muii0,kt) exp(-(murg-muii0).^2/2./murg/kt)./sqrt(murg*kt)/sqrt(2*pi);
% K = round(gap0^2*(255:(-1):1)/255+100);
% K = round(gap0^2*exp(-(1:255)/50));
K = gap0^2*ones(1,255);
% K(86:116) = 0;
murg = 1:255;
varPred = zeros(1,255);
for ii0=1:255
    Nele = pdf0(murg,ii0,kt);
    N0 = sum(K.*Nele);
    N1 = 1/2/N0*sum(((murg-ii0).^2+kt*ii0).*Nele.*K);
    varPred(ii0) = N1;
end

%% plot all
figure;
plot(meanRg1,xVar1,'b');
hold on
plot(meanRg1,xVar2a,'r');
plot(meanRg1,xVar2b,'k');
plot(meanRg1,xVar2c,'g');
plot(meanRg1,xVar2d,'c');
plot(meanRg1,varPred,'m');
hold off
xlabel('Intensity');
ylabel('Variance');
t=title('Simulated patches');
set(t,'Interpreter','none');
legend('Multiple','Single mean','Single mean average',...
    'Single mean correction','Single one direction',...
    'Predicted for one direction','Location','SouthEast')


figure;errorbar(meanRg1,xVar1,xVar1Sd);
xlabel('Intensity');
ylabel('Variance with standard deviation');
t=title('Simulated patches');
set(t,'Interpreter','none');



