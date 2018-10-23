%% Bead signal from grace
figId = 1;
meanRg = (0:255)/255;
meanRg1 = meanRg(2:end);
N1 = length(meanRg)-1;
[fname,stub] = uigetfile('*.lsm','Select the .lsm file you are analyzing');
path1 = fullfile(stub,fname);
[dat,datRef,dats] = readStack(path1,0);
dat0 = dat(:,:,figId);

%% Reference
% foi's method
[fitparams,cx,cy] = get_var_foi(dat0);
xVar3 = meanRg1*fitparams(1) + fitparams(2);

% Variance and mean by time series
[xVar1,xVar1Var,xVar1Sd,nPix1] = getVarCurveMulti(dats,meanRg);

% yMax = max(xVar1);
% figure;
% plot(meanRg1,xVar1,'b','LineWidth',2);ylim([-yMax*0.1, yMax*1.2]);xlim([0,1.1])
% hold on
% % plot(meanRg1,xVar3,'r');
% plot(cx{4},cy{4}.^2,'r','LineWidth',2)
% hold off
% set(gca, 'FontSize', 12)
% xlabel('Intensity','FontSize',14);
% ylabel('Variance','FontSize',14);
% h_legend = legend('100 images','Single image','Location','northwest');
% set(h_legend,'FontSize',14);

%% my implementation of Foi
edgeTau = 1;
nLevels = 600;
delta0 = 1/nLevels;

% wavelet and scaling functions
psi1 = [0.035 0.085 -0.135 -0.460 0.807 -0.333];
phi1 = [0.025 -0.06 -0.095 0.325 0.571 0.235];
phi2 = (conv2(phi1.',phi1)).^2;
phi22 = sum(phi2(:));
[zwapp,cH,cV,zwdet] = dwt2(dat0,phi1,psi1);

% smoothing and edge detection
hw = fspecial('average', 7);
hl = fspecial('laplacian',1);
hs = fspecial('sobel');
zsmo = imfilter(zwapp,hw);
s = imfilter(sqrt(pi/2)*abs(zwdet),hw);
zwappMedian = medfilt2(zwapp,[3,3]);
zwappL = imfilter(zwappMedian,hl);
zwappS = imfilter(zwappL,hs);
xsmo = (abs(zwappS) + abs(zwappL)) < edgeTau*s;
K_xsmo = cat(3,xsmo,zwapp,zwapp*0);
% figure;imshow(K_xsmo);

% level sets
levelSets = {};
zsmo1 = zsmo;
zsmo1(xsmo==0) = -1;
zsmo1Max = max(zsmo1(:));
for ii=1:nLevels
    z0 = (ii-1)*delta0;
    z1 = z0 + delta0;
    if z0>zsmo1Max
        break
    end
    zIdx = find(zsmo1>=z0 & zsmo1<z1);
    if length(zIdx)>20
        levelSets{ii} = zIdx;
    end
end
nsi = cellfun(@length,levelSets);
% figure;semilogy(nsi);

% variance mean pairs
nPair = sum(nsi>0);
xMean = zeros(1,nPair);
xStd = zeros(1,nPair);
xVar = zeros(1,nPair);
xKai = zeros(1,nPair);
xCi = zeros(1,nPair);
xDi = zeros(1,nPair);
xNi = zeros(1,nPair);
n0 = 1;
for ii=1:length(levelSets)
    pix0 = levelSets{ii};
    if ~isempty(pix0)
       xMean(n0) = mean(zwapp(pix0));
       zwdet0 = zwdet(pix0);
       zwdet0Mean = mean(zwdet0);
       ni = length(zwdet0);
       xNi(n0) = ni;
       kai = 1-1/4/ni-7/32/ni^2;
       xKai(n0) = kai;
       xStd(n0) = sqrt(sum((zwdet0 - zwdet0Mean).^2)/(ni-1))/kai;
       xVar(n0) = sum((zwdet0 - zwdet0Mean).^2)/(ni-1);
       xCi(n0) = phi22/ni;
       xDi(n0) = (1-kai^2)/(kai^2);
       n0 = n0 + 1; 
    end
end
% figure;scatter(xMean,xStd)
% figure;scatter(xMean,xVar)

% LS Estimation
A = zeros(nPair,2);
A(:,1) = xMean;
A(:,2) = 1;
b = reshape(xVar,[],1);
pEst = (A.'*A)\(A.'*b);

% ML Estimation
pEst = fminsearch(@(x) pgL(x,xMean,xStd,xCi,xDi), pEst);
xVar4 = meanRg1*pEst(1) + pEst(2);

%% variance stablization
% Anscombe, generalized Anscombe or square root
% pEst0 = fitparams;
pEst0 = pEst;
useAnscombe = 1;
alpha = pEst0(1);
sigma2 = max(pEst0(2),1e-6);
datsa = dats*0;
if useAnscombe
    t0 = 2/alpha*sqrt(alpha*0+3/8*alpha^2+sigma2);
    t1 = 2/alpha*sqrt(alpha*1+3/8*alpha^2+sigma2);
else
    ofst = 0;
    pwrr = 0.5;
    cha1 = 1/pEst0(1);
    t0 = 2*(ofst)^pwrr;
    t1 = 2*(ofst+1*cha1)^pwrr;
end
for ii=1:size(dats,2)
    if useAnscombe
        datsa(:,ii) = 2/alpha*sqrt(alpha*dats(:,ii)+3/8*alpha^2+sigma2);
    else
        datsa(:,ii) = 2*(dats(:,ii)*cha1+ofst).^pwrr;
    end
end
datsan = datsa;
datsan = (datsan-t0)/(t1-t0);
meanRg2 = meanRg;
% meanRg2 = t0:(t1-t0)/255:t1;
[xVar1a,xVar1Vara,xVar1Sda,nPix1a] = getVarCurveMulti(datsan,meanRg2);

%% plot all
yMax = max(xVar1);
figure;
plot(meanRg1,xVar1,'b');ylim([-yMax*0.1, yMax*1.2]);xlim([0,1.1])
hold on
% plot(meanRg1,xVar3,'r');
plot(cx{4},cy{4}.^2,'r')
plot(meanRg1,xVar4,'k')
scatter(xMean,xVar,5,'m')
plot(meanRg1,xVar1a,'g')
hold off
xlabel('Intensity');
ylabel('Variance');
h_legend = legend('Multiple images','Single image (Foi)',...
    'Single image (Foi, our implement)','Mean-variance Pairs','After VST','Location','northwest');













