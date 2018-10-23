function handles = beadDetectAlignExtract( handles )

pp = handles.pp;
dd = handles.dd;

%% prepare data
datBase = [];   % f0
if isempty(dd.img0)
    dat = dd.img1;
elseif isempty(dd.img1)
    dat = dd.img0;
else
    t0 = size(dd.img0,3);
%     t1 = size(dd.img1,3);
    if t0==1
        dat = dd.img1;
        datBase = single(dd.img0)/pp.maxVal;
    else
        dat = dd.img0;
        datBase = single(dd.img1)/pp.maxVal;
    end
end

[Nx,Ny] = size(dat(:,:,1,1));
bgMean = zeros(Nx,Ny);
% nFrame = size(dat,3);
rgx = pp.timeRg0:pp.timeRg1;
nFrame = length(rgx);
meanCurve = zeros(nFrame,1);
for nn=1:nFrame
    ii = rgx(nn);
    bgMean = bgMean + single(dat(:,:,ii,pp.bgch))/pp.maxVal;
    datEle = single(dat(:,:,ii,pp.sigch))/pp.maxVal;
    meanCurve(nn) = mean(datEle(:));
end
bgMean = bgMean/nFrame;
bg1 = single(dat(:,:,1,pp.bgch))/pp.maxVal;
bgn = single(dat(:,:,nFrame,pp.bgch))/pp.maxVal;

%% ROI detection
% [resBead,resBeadCenter,resBeadRad,beadShow] = beadDetection(bgMean,pp.thrDetect,...
%     pp.thrMulti,pp.radRg0,pp.radRg1,pp.usePara,pp.outPath);
[resBead,resBead1,resBeadCenter,resBeadRad,beadShow] = beadDetection(bgMean,...
    pp.smoMethod,pp.thrMulti,pp.thrDetect,pp.radRg0,pp.radRg1,pp.outPath);

%% select control beads
% f0 = figure('name','Click to select the control beads centers','NumberTitle','off');
set(handles.myinfo,'String','Left click to pick control beads, right click to remove; click done to continue');
imageHandle = imshow(beadShow);
set(imageHandle,'ButtonDownFcn',{@controlBeadClick});
% set(f0,'CloseRequestFcn',@controlBeadClose);
btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
    'Position', [20 20 60 25], 'Callback', @controlBeadClose);
uiwait;
controlBeadCenters = getappdata(imageHandle,'test');
delete(btn);
set(handles.myinfo,'String','Please continue');
pause(0.25)
fprintf('done\n');

%% Alignment and return registered background image
sig = dat(:,:,rgx,pp.sigch);
% sig = dat(:,:,:,pp.sigch);
if ~isempty(datBase)
    % alignment
    baseBg = datBase(:,:,1,pp.bgch);
    baseSig = datBase(:,:,1,pp.sigch);
    meanCurve = [meanCurve;mean(baseSig(:))];
    h = msgbox('Aligning images...');
    [baseBg_rigid1,baseSig_rigid1,qual1] = beadAlignment(bgMean, baseBg, baseSig, pp.outPath, 1);
    [baseBg_rigid0,baseSig_rigid0,qual0] = beadAlignment(bgMean, baseBg, baseSig, pp.outPath, 0);
    if qual1>qual0
        baseBg_rigid = baseBg_rigid1;
        baseSig_rigid = baseSig_rigid1;
    else
        baseBg_rigid = baseBg_rigid0;
        baseSig_rigid = baseSig_rigid0;
    end
    sig = cat(3,uint16(round(baseSig_rigid*pp.maxVal)),sig);
    close(h);
    
    % remove beads only exist in base or sig
    beadVal = ones(length(resBead),1);
    for ii=1:length(resBead)
        p0 = resBead{ii};
        p0 = sub2ind(size(bgMean),p0(:,1),p0(:,2));
        beadVal(ii) = median(bgMean(p0));
%         beadVal(ii) = mean(bgMean(p0));
    end
    beadMin = max(min(beadVal),0.05);
    
    idxSel = ones(length(resBead),1);
    baseBg_rigid_smo = wiener2(baseBg_rigid,[5,5]);
    bg1_smo = wiener2(bg1,[5,5]);
    bgn_smo = wiener2(bgn,[5,5]);
    for ii=1:length(resBead)
        p0 = resBead{ii};
        p0 = sub2ind(size(bgMean),p0(:,1),p0(:,2));
        dif01 = sum((baseBg_rigid_smo(p0)-bg1_smo(p0)).^2)/sum((baseBg_rigid_smo(p0)).^2);
        dif0n = sum((baseBg_rigid_smo(p0)-bgn_smo(p0)).^2)/sum((baseBg_rigid_smo(p0)).^2);
%         dif0 = sum((baseBg_rigid_smo(p0)-bgMean(p0)).^2)/sum((bgMean(p0)).^2);
%         corx0 = corrcoef(baseBg_rigid_smo(p0),bgMean(p0));
%         baseVal = mean(baseBg_rigid(p0));
%         timeVal = mean(bgMean(p0));
        baseVal = median(baseBg_rigid_smo(p0));
%         timeVal = median(bgMean(p0));
%         rt0 = baseVal/(timeVal+1e-8);
%         rt1 = timeVal/(baseVal+1e-8);
        if baseVal<beadMin || dif01>0.1 || dif0n>0.1
            idxSel(ii) = 0;
        end
%         if rt0 > 2 || rt1 >2
%             idxSel(ii) = 0;
%         end
%         if baseVal<beadMin || rt0 > 1.1 || rt1 > 1.1
%             idxSel(ii) = 0;
%         end
    end
    resBeadx = resBead(idxSel>0);
    resBead1x = resBead1(idxSel>0);
    resBeadCenterx = resBeadCenter(idxSel>0,:);
    resBeadRadx = resBeadRad(idxSel>0);
    % dump again
    resBorder = bgMean*0;
    neibVec = [0, -1, 1, -Nx, Nx];    
    fprintf('Dump valid beads ===== \n');
    for ii=1:length(resBeadx)
        if mod(ii,1000)==0
            fprintf('ii is %d\n',ii);
        end
        idx = resBeadx{ii};
        idx = sub2ind(size(bgMean),idx(:,1),idx(:,2));
        idxk = bsxfun(@plus,idx,neibVec);
        idxa = ismember(idxk,idx);
        idxSel = sum(idxa,2)<5;
        idxc = idx(idxSel);
        resBorder(idxc) = 0.75;
    end
    K2 = cat(3,resBorder,bgMean,bgMean*0);
    fname = 'res_bead_valid';
    imwrite(double(K2),[pp.outPath,filesep,fname,'.tif']);    
    Ka = cat(3,baseBg_rigid,bg1,resBorder);
    fname = 'alignment_circled_t1';
    imwrite(double(Ka),[pp.outPath,filesep,fname,'.tif']);
    Kb = cat(3,baseBg_rigid,bgn,resBorder);
    fname = 'alignment_circled_tn';
    imwrite(double(Kb),[pp.outPath,filesep,fname,'.tif']);
end

%% curve extaction
h = msgbox('Extracting curves...');
if pp.extractMethod==1
    fprintf('Extract curve with ring\n');
    myCurves = beadExtract(sig, resBeadCenterx, resBeadRadx, pp.radRg0,pp.radRg1,pp.outPath,pp.maxVal,1);
elseif pp.extractMethod==2
    fprintf('Extract curve with mean\n');
    myCurves = beadExtractOtsu(sig, resBeadx, pp.outPath, pp.maxVal);
else
    warning('Unsupported method\n');
    myCurves = beadExtractNMF(sig, resBeadx, pp.outPath, pp.maxVal);
end

% myCurves = beadExtractFANorm(sig, resBeadx, pp.outPath, pp.maxVal, 1, 1);
% myCurves = beadExtractNMF(sig, resBeadx, pp.outPath, pp.maxVal);
% myCurves = beadExtractIter(sig, resBeadx, pp.outPath, pp.maxVal);
if isvalid(h)
    close(h);
end

% save results
handles.dd.myCurves = myCurves;
handles.dd.meanCurve = meanCurve;
handles.dd.resBead = resBeadx;
handles.dd.bgMean = bgMean;
handles.dd.sig = sig;
handles.dd.controlBeadCenters = controlBeadCenters;

end

