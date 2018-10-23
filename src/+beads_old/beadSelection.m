%Bead selection pipeline
% debug only
addpath(genpath('./src/'))
addpath(genpath('./output/'))
addpath(genpath('./tools/'))
warning('off', 'Images:initSize:adjustingMag');

run_detect = 1;
run_fitting = 1;
datMode = 0; % 0: deltaF, 1: deltaF/F

%% load files
[fname,stub] = uigetfile('*.lsm','Select the background .lsm file');
path0 = fullfile(stub,fname);
[fname,stub] = uigetfile(fullfile(stub,'*.lsm'),'Select the .lsm file you are analyzing');
path1 = fullfile(stub,fname);

% background
dat0 = tiffread(path0);
maxVal = 2^dat0.bits-1;
bgch1 = double(dat0.data{1})/maxVal;
bgch2 = double(dat0.data{2})/maxVal;
[Nx,Ny] = size(bgch1);

% time lapse
dat1 = tiffread(path1);
nFrame = length(dat1);
sigch2 = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
    tmp = dat1(ii);
    datEle = double(tmp.data{1})/maxVal;
    if ii==1
        sigch1Mean = datEle;
    else
        sigch1Mean = sigch1Mean + datEle;
    end
    datEle = double(tmp.data{2})/maxVal;
    sigch2(:,:,ii) = datEle;
end
sigch1Mean = sigch1Mean/nFrame;

%% bead detection
% crop_rgx = 1801:2400;
% crop_rgy = 1701:2200;
% sigch1MeanCrop = sigch1Mean(crop_rgx,crop_rgy);
sigch1MeanCrop = sigch1Mean;
resBead = beadDetection(sigch1MeanCrop,0.1);

%% post processing
% Alignment and curve extaction
myCurves = beadRegistrationCurve(resBead, sigch1Mean, sigch2, bgch1, bgch2, nFrame);

% Compute deltaF-F0 and fit the curve
[dff_fit,tau_fit,dff_mat,time_vect] = beadProcessor(myCurves,path1,datMode);
goodBeadIdx = beadChoose(dff_fit,tau_fit,dff_mat,time_vect,datMode);

% plot them
beadsDraw(resBead,sigch1Mean,goodBeadIdx)





