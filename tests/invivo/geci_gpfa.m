%% load data
USE_SIM = 1;

if USE_SIM==1
    xCurve = zeros(1,2000);
%     xCurve(150:160) = 0.8*sin(pi*(0:0.1:1)); % slow
%     xCurve(150:154) = 0.8*sin(pi*(0:0.25:1)); % fast
    xCurve  = xCurve + 0.2;
    datRoiSel = repmat(xCurve,44,1);
    datRoiSel = datRoiSel + randn(size(datRoiSel))*0.3;
    datRoiSel(datRoiSel>1) = 1;
    datRoiSel(datRoiSel<0) = 0;
else
    ftop = 'D:\neuro_WORK\coculture_ucdavis_WORK\';
    fname = 'nn-4001-sv4';
    runName = ['geci_',fname,'-run-0908-mean\'];
    
    % -- load ROIs
    % x = load([ftop runName,'res_cor.mat']);
    x = load([ftop 'GECI_test\' runName,'res_temp.mat']);
    roiMap = x.roi_map;
    
    % -- load data
    fTiffName = [ftop,'GECI_tiff\',fname,'.tif'];
    h = imfinfo(fTiffName);
    nW = h(1).Height;
    nH = h(1).Width;
    nT = length(h);
    maxVal = 2^h(1).BitsPerSample - 1;
    dat = zeros(nW,nH,nT);
    for ii=1:length(h)
        temp = double(imread(fTiffName,ii,'Info',h));
        temp = temp/maxVal;
        dat(:,:,ii) = temp;
    end
    
    % -- curves for ROI
    idNow = 103;
    [iX,iY] = find(roiMap == idNow);
    nPix = length(iX);
    datRoi = zeros(nPix,nT);
    for ii=1:nPix
        datRoi(ii,:) = dat(iX(ii),iY(ii),:);
    end
%     datRoiSel = datRoi(:,201:600);
    % datRoiSel = sqrt(datRoi);
    datRoiSel = datRoi;
end

%% GFPA
fname = 'output/roi_sim';
xDim = 1;
binWidth = 1000;
startTau = 100000;
startEps = 1e-3;
extraOpts = {'kernSDList', 30};
seqTrain = struct('trialId',{},'T',{},'y',{});
seqTest = seqTrain;
seqTrain(1).T = size(datRoiSel,2);
seqTrain(1).trialId = 1;
seqTrain(1).y = datRoiSel;
gpfaEngine(seqTrain, seqTest, fname,...
    'xDim', xDim, 'binWidth', binWidth, 'startTau', startTau, 'startEps', startEps, extraOpts{:});

%% plot
load(fname);
x0 = seqTrain.xsm';
x1 = mean(datRoiSel,1);
x0 = x0 - mean(x0(:));
x0 = x0/max(x0(:));
x1 = x1 - mean(x1);
x1 = x1/max(x1);
figure;
subplot(2,1,1);plot(x0);title('GPFA');ylim([-1 1])
% legend('1','2');
subplot(2,1,2);plot(x1);title('Mean');ylim([-1 1])

currentParams
estParams




