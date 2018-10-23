%% bead detection for selected areas
pathTop = 'D:\neuro_WORK\peptide_ucdavis_WORK\dat\082715 bead detection performance\';
pathLbTop = 'D:\neuro_WORK\peptide_ucdavis_WORK\dat\083115_grace_labelling';
pathOutTop = 'D:\neuro_WORK\peptide_ucdavis_WORK\dump\paper_test\';

Nimg = 3;  % 3 or 10 images to test
plotme = 0;

% reference paramters
thrMulti = 3;  % 3
thrDetect = 0.05;
radRg0 = 4;
radRg1 = 15;
noiseAmp = 0;

%% score map
if 0
    fp1 = [pathOutTop,'score\img1_smo3_cht_res.mat'];
    tmp = load(fp1);
    K1 = tmp.K2 * 0.3;
    K1(:,:,2) = K1(:,:,2)*2;
    for tt=1:6
        fp0 = [pathOutTop,'score\img1_smo3_cht_iter',num2str(tt),'_accum.mat'];        
        tmp = load(fp0);
        accu1 = tmp.accu1;
        accu1 = accu1./max(accu1(:));
        K1(:,:,3) = accu1*2;
        figure;imshow(K1);
    end
end

%% impact of smoothing
smoRg = [0 1 3 5 7];
% smoRg = 3;
Ncond = length(smoRg);
resPrecisionA = zeros(Ncond,Nimg);
resRecallA = zeros(Ncond,Nimg);
resMSEA = zeros(Ncond,Nimg);
for kk=1:Ncond
    smo1 = smoRg(kk);
    [ resPrecision,resRecall,resMSE ] = beads.paperDetectSim( pathTop,pathLbTop,...
        pathOutTop,smo1,thrDetect,radRg0,radRg1,Nimg,plotme,noiseAmp);
    resPrecisionA(kk,:) = resPrecision;
    resRecallA(kk,:) = resRecall;
    resMSEA(kk,:) = resMSE;
end
resFsA = 2*resPrecisionA.*resRecallA./(resPrecisionA+resRecallA);

figure;plot(smoRg,resFsA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('Smoothing radius');
ylabel('F_1 score');

%% impact of thresholds
thrDetectRg = [0.02 0.05 0.1 0.15 0.2];
Ncond = length(thrDetectRg);
resPrecisionA = zeros(Ncond,Nimg);
resRecallA = zeros(Ncond,Nimg);
resMSEA = zeros(Ncond,Nimg);
for kk=1:Ncond
    thr1 = thrDetectRg(kk);
    [ resPrecision,resRecall,resMSE ] = beads.paperDetectSim( pathTop,pathLbTop,...
        pathOutTop,thrMulti,thr1,radRg0,radRg1,Nimg,plotme,noiseAmp );
    resPrecisionA(kk,:) = resPrecision;
    resRecallA(kk,:) = resRecall;
    resMSEA(kk,:) = resMSE;
end
resFsA = 2*resPrecisionA.*resRecallA./(resPrecisionA+resRecallA);

figure;plot(thrDetectRg,resFsA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('Detection threshold');
ylabel('F_1 score');

figure;plot(thrDetectRg,resPrecisionA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('Detection threshold');
ylabel('Precision');

figure;plot(thrDetectRg,resRecallA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('Detection threshold');
ylabel('Recall');

%% impact of radius range
rg1Rg = 2:6;
rg2Rg = 8:20;
Ncond1 = length(rg1Rg);
Ncond2 = length(rg2Rg);
resPrecisionA = zeros(Ncond1,Ncond2);
resRecallA = zeros(Ncond1,Ncond2);
resMSEA = zeros(Ncond1,Ncond2);
for k1=1:Ncond1
    fprintf('%d\n',k1);
    for k2=1:Ncond2
        r1 = rg1Rg(k1);
        r2 = rg2Rg(k2);
        [ resPrecision,resRecall,resMSE ] = beads.paperDetectSim( pathTop,pathLbTop,...
            pathOutTop,thrMulti,thrDetect,r1,r2,1,plotme,noiseAmp );
        resPrecisionA(k1,k2) = resPrecision;
        resRecallA(k1,k2) = resRecall;
        resMSEA(k1,k2) = resMSE;
    end
end
resFsA = 2*resPrecisionA.*resRecallA./(resPrecisionA+resRecallA);

[gridx,gridy] = meshgrid(rg2Rg,rg1Rg);
figure;surf(gridx,gridy,resFsA);
xlabel('Max radius');
ylabel('Min radius');
zlabel('F_1 score');

figure;surf(gridx,gridy,resPrecisionA);
xlabel('Max radius');
ylabel('Min radius');
zlabel('Precision');

figure;surf(gridx,gridy,resRecallA);
xlabel('Max radius');
ylabel('Min radius');
zlabel('Recall');

%% noise levels estimation
a = zeros(Nimg,1);
b = zeros(Nimg,1);
for ii=1:Nimg
    fname = [num2str(ii),'.tif'];
    fpname = [pathTop,fname];
    X0 = tiffread(fpname);
    bgMean = double(X0(1).data)./65535;
    fitparam = noise_vst.get_var_foi(bgMean);
    a(ii) = fitparam(1);
    b(ii) = fitparam(2);
end

%% impact of poisson noise levels
noiseRg = [0 0.01 0.05 0.1 0.2];
noisedbRg = -10*log10(noiseRg.^2);
Ncond = length(noisedbRg);
resPrecisionA = zeros(Ncond,Nimg);
resRecallA = zeros(Ncond,Nimg);
resMSEA = zeros(Ncond,Nimg);
for kk=1:Ncond
    noise0 = noiseRg(kk);
    [ resPrecision,resRecall,resMSE ] = beads.paperDetectSim( pathTop,pathLbTop,...
        pathOutTop,thrMulti,thrDetect,radRg0,radRg1,Nimg,plotme,noise0 );
    resPrecisionA(kk,:) = resPrecision;
    resRecallA(kk,:) = resRecall;
    resMSEA(kk,:) = resMSE;
end
resFsA = 2*resPrecisionA.*resRecallA./(resPrecisionA+resRecallA);

figure;plot(noisedbRg,resFsA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('-10log10(noise variance)');
ylabel('F_1 score');

figure;plot(noisedbRg,resPrecisionA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('-10log10(noise variance)');
ylabel('Precision');

figure;plot(noisedbRg,resRecallA)
legend('Image 1','Image 2','Image 3','location','southeast');
xlabel('-10log10(noise variance)');
ylabel('Recall');


%% impact of gaussian noise levels










