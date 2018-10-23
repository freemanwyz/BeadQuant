% files parameters
pp.path0 = [];
pp.path1 = [];
pp.fname0 = [];
pp.fname1 = [];
pp.ftype0 = [];
pp.ftype1 = [];
pp.outPath = 'output';

% image parameters
pp.imgNum = [];      % image number, if 2, use first one as f0, if 1, use first time point as f0
pp.timeRg0 = 1;
pp.timeRg1 = 1;
pp.timeStep = 240;
pp.nChan = [];       % number of channels
pp.bgch = 1;        % back grouand channel in image
pp.sigch = 2;       % signal channel in image
pp.maxVal = 2^16-1;

% detection parameters
pp.thrDetect = 0.05;    % detection threshold
pp.thrMulti = 3;    % smoothing window size
pp.radRg0 = 4;
pp.radRg1 = 13;

% general parameters
pp.useDff = 0;      % 0: deltaF, 1: deltaF/F
pp.fitModel = 1;    % 1: exponential, 2: linear (2 points)
pp.doneStage = 0;
pp.usePara = 1;
pp.smoMethod = 1;   % 1: Wiener, 2: Gaussian

% data and results
dd.img0 = [];
dd.img1 = [];
dd.myCurves = [];
dd.meanCurve = [];
dd.resBead = [];
dd.bgMean = [];
dd.sig = [];

dd.dff_fit = [];
dd.tau_fit = [];
dd.dff_mat = [];
dd.dff_fit_mean = [];
dd.tau_fit_mean = [];
dd.time_vect = [];
dd.goodBeadIdx = [];
dd.controlBeadIdx = [];
dd.controlBeadCenters = [];

save('./cfg/cfgDefault.mat','pp','dd')