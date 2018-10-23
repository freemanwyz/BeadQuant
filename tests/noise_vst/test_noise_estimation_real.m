%% Bead signal from grace
p = getInitParamBead();
figId = 10;
meanRg = (0:255)/255;
% meanRg = 0:255;
meanRg1 = meanRg(2:end);
N1 = length(meanRg)-1;
p.crop = 0;

%% data
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
% dat = dat*255;
% datRef = round(datRef*255);
dats = reshape(dat,[],nFrame);

%% foi's method
if 1
    dat0 = dat(:,:,figId);
    fitparams = get_var_foi(dat0);
    %     fitparams = get_var_foi(dat0/255);
    xVar3 = meanRg1*fitparams(1) + fitparams(2);
end

%% Variance and mean by time series
fprintf('Multi ----- \n');
[xVar1,xVar1Var,xVar1Sd,nPix1] = getVarCurveMulti(dats,meanRg);

%% Variance and mean by single image
fprintf('Mean ----- \n');
dat0 = dat(:,:,figId);
% dat0 = wiener2(dat0,[5,5]);
[xVar2a,nPix2a,nPairs2a,xMean2a] = getVarCurve(dat0,meanRg,1,'mean',2);
[xVar2b,nPix2b,nPairs2b,xMean2b] = getVarCurve(dat0,meanRg,1,'meanAvg',2);

if 0
    dat0 = dat(:,:,figId);
    [Gmag,Gdir] = imgradient(dat0);
    dat0(Gmag>0.4) = Inf;
    % dat0(Gmag>(0.4*255)) = Inf;
    [xVar2c,nPix2c,nPairs2c] = getVarCurve(dat0,meanRg,1,'mean',2);
    % [xVar2d,nPix2d,nPairs2d] = getVarCurve(dat0,meanRg,1,'meanSingleDirection',2);
end

%% analytical
if 0
    K = nPix1;
    kt = 6.2;  % study from multiple
    pdf0 = @(murg,muii0,kt) exp(-(murg-muii0).^2/2./murg/kt)./sqrt(murg*kt)/sqrt(2*pi);
    murg = 1:255;
    varPred = zeros(1,255);
    for ii0=1:255
        Nele = pdf0(murg,ii0,kt);
        N0 = sum(K.*Nele);
        N1 = 1/2/N0*sum(((murg-ii0).^2+kt*ii0).*Nele.*K);
        varPred(ii0) = N1;
    end
end

%% analysis 1
ii0 = 30;
x1 = meanRg(ii0);
x2 = meanRg(ii0+1);
maskIn = find(dat0>=x1 & dat0<x2);
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
noiseMean = nanmean(zx,2);
noiseMean = noiseMean(~isinf(noiseMean) & ~isnan(noiseMean));

noiseAll = zx(:);
noiseAll = noiseAll(~isinf(noiseAll) & ~isnan(noiseAll));

var(noiseMean)/9*8
var(noiseAll)/2
figure;hist(noiseMean,50);
figure;hist(noiseAll,50);

%% plot all
yMax = max(xVar1);
figure;
plot(meanRg1,xVar1,'b');ylim([-yMax*0.1, yMax*1.2]);xlim([0,1.1])
hold on
plot(meanRg1,xVar3,'r');
% plot(meanRg1,xVar2a,'g');
% plot(meanRg1,xVar2b,'k');
% plot(meanRg1,xMean2a/10,'m');
% plot(meanRg1,xVar2c,'g');
% plot(meanRg1,xVar2d,'g');
% plot(meanRg1,varPred,'m');
hold off
xlabel('Intensity');
ylabel('Variance');
% t=title(fname);
% set(t,'Interpreter','none');
% legend('Multiple','Foi','Single Mean','Single Mean Avg',...
%     'Location','NorthWest')
% legend('Multiple','Single mean','Single mean average','Single mean no edge',...
%     'Predicted for one direction',...
%     'Location','NorthWest')

%%
% idxk0 = find(dat0>=0.5 & dat0<0.6);
% tmp = dat0*0;
% tmp(idxk0) = 1;
% K0 = cat(3,tmp,dat0,dat0*0);
% figure;imshow(K0)

figure;plot(meanRg1,xMean2a,'r');
hold on
plot(meanRg1,xMean2b,'b');

figure;
plot(meanRg1,log10(nPairs2a))


