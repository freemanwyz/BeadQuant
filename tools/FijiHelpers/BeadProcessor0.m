%Bead ProcessingScript
%It will then automatically fit time series for each bead with a decaying
%exponential, returning the asymptotic value as a max DFF estimate and tau
%as an estimation of speed with which bead changes fluorescence.
%
%User then selects beads based on desirable DFF and tau combinations.
%Output is series of graphs with fits for all slected beads.

addpath(genpath('./data/'))
addpath(genpath('./output/'))
addpath(genpath('./FijiHelpers/'))
warning('off', 'Images:initSize:adjustingMag');

%% Read in the .lsm file/data in FIJI results window
[fname,stub] = uigetfile('*.lsm','Select the .lsm file you are analyzing');
path = fullfile(stub,fname);
[lsmA,lsmB,lsmC] = lsminfo(path);
%Use this if the timeStamp aquisition reflects 
% time_prep = lsmA.TimeStamps.TimeStamps+ lsmA.TimeStamps.AvgStep;
% time_vect = [0;time_prep];

fin = './output/resCurve.mat';
tmp = load(fin);
M = transpose(tmp.myCurves);

%Making the time vector using size of the M matrix
time_prep = 0:size(M,1)-1;
time_vect = time_prep*lsmA.TimeStamps.AvgStep;

%Calculate and plot DFF for all beads
datMode = 0;
dff_mat = calcDF(M,1,datMode,1);
figure
plot(time_vect,dff_mat);
xlabel('time (sec)');
if datMode
    ylabel('\DeltaF/F');
else
    ylabel('\DeltaF');
end

%% Fit all of the data. Return maxDFF, tau, and R^2
dff_fit = zeros(1,size(dff_mat,2));
tau_fit = zeros(1,size(dff_mat,2));
parfor i = 1:size(dff_mat,2)
    try
        [cf,gof] = DecayingExponentialFit(time_vect,dff_mat(:,i));
        %Inputting the fit evaluated at the last time point for max
        %activation
        fitVal = feval(cf,time_vect(end));
        dff_fit(i) = fitVal;
        %Using the tau fit for speed of activation
        tau_fit(i) = cf.tau;
        fprintf('%d\n',i);
    catch
        dff_fit(i) = 0;
        tau_fit(i) = 0;
    end
end

scatter_x = 1:size(dff_mat,2);

%% Visualize Data
%Plot ROIs with bad fits
badROI = find(dff_fit==0&tau_fit==0);

figure;
plot(time_vect,dff_mat(:,badROI));
xlabel('time (sec)');
if datMode
    ylabel('\DeltaF/F');
else
    ylabel('\DeltaF');
end
title('Plots of ROIs with failed fits');

%% Select region of points for analysis
%Select the region of plot where interesting beads are
f = figure;
scatter(dff_fit,tau_fit);
ylabel('\tau')
xlabel('\DeltaF/F')
h = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
title('Zoom and center plot on points of interest, then hit "Continue"')

uiwait(f); 
axLim = [get(gca,'xlim'),get(gca,'ylim')];
close(f);

%Select beads of interest by enclosing within a rectangle
FIGDisp = 1;%loop until user happy with rectangle
while FIGDisp==1
    f = figure;
    scatter(dff_fit,tau_fit);
    axis(axLim)
    ylabel('\tau')
    xlabel('\DeltaF/F')
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    plot(x,y,'r','linewidth',2)          % draw box around selected region
    prompt = 'Input 1 to select a different box, 0 to continue';
    FIGDisp = str2num(cell2mat(inputdlg(prompt)));
    close(f);
end

%% Get and make plots for selected beads
%Get selected beads
goodDFF = dff_fit>min(x) & dff_fit<max(x);
goodTau = tau_fit>min(y) & tau_fit<max(y);
goodBeads = goodDFF & goodTau;
goodBeadIdx = find(goodBeads);

%make plots for selected beads
for ii = 1:sum(goodBeads)
    figure
    DecayingExponentialFitwGraph(time_vect,dff_mat(:,goodBeadIdx(ii)));
    xlabel('time (sec)')
    ylabel('\DeltaF/F')
    titStr = sprintf('Plot for Bead %d',goodBeadIdx(ii));
    title(titStr)
end

%% Draw the selected beads
load('./output/resBeads.mat');
datAvg = double(imread('./output/beadAvg16.tif'))/65535;
[Nx,Ny] = size(datAvg);
xI = goodBeadIdx;

neibVec = [0, -1, 1, -Nx, Nx];
t0 = zeros(Nx,Ny);
for ii=1:length(xI)
    idx = resBead{xI(ii)};
    idx = sub2ind(size(t0),idx(:,1),idx(:,2));
    idxk = bsxfun(@plus,idx,neibVec);
    idxa = ismember(idxk,idx);
    idxSel = sum(idxa,2)<5;
    idxc = idx(idxSel);    
    t0(idxc) = 1;
end
K0 = cat(3,t0,datAvg,datAvg*0);
figure;imshow(K0);






