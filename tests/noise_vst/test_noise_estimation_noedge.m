p = getInitParam('preset','yinxue_nucleus');
fp0 = [p.tp filesep 'data' filesep 'dat_test' filesep 'C2-4-DS4-cherry-s100b-3.tif'];
% p.ch = 2;
% fp0 = [p.tp filesep 'data' filesep '20150519_Tuj1_SynI\neuron only\' filesep '2-weak-z_Maximum intensity projection.lsm'];
[Img, ImgCrop, q] = loadData(fp0, 'tif', p);

meanRg = 0:0.01:1;
N1 = length(meanRg)-1;

%% 1 mean, noedge
dat0 = Img/255;
% dat0 = datRef;
[Gmag,Gdir] = imgradient(dat0);
dat0(Gmag>0.5) = Inf;
% K0 = cat(3,Gmag>1,dat0,dat0*0);
% figure;imshow(K0)
xVar3 = zeros(1,N1);
nPix3 = zeros(1,N1);
parfor ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(dat0>=x1 & dat0<x2);
    nPix3(ii) = length(idx);
    [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
    xVar3(ii) = varOut;
end

%% 2 mean
dat0 = Img/255;
xVar1 = zeros(1,N1);
nPix1 = zeros(1,N1);
parfor ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(dat0>=x1 & dat0<x2);
    nPix1(ii) = length(idx);
    [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'mean',1 );
    xVar1(ii) = varOut;
end

%% 3 median
dat0 = Img/255;
xVar2 = zeros(1,N1);
nPix2 = zeros(1,N1);
parfor ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(dat0>=x1 & dat0<x2);
    nPix2(ii) = length(idx);
    [ varOut, nPairs ] = getVarNeib( dat0,idx,0,'median',1 );
    xVar2(ii) = varOut;
end

%% plot
meanRg1 = meanRg(1:(end-1));
figure;
plot(meanRg1,xVar1/2,'r');ylim([0,0.05]);
hold on
plot(meanRg1,xVar2/2,'b');ylim([0,0.05]);
plot(meanRg1,xVar3/2,'k');ylim([0,0.05]);
hold off
legend('Single','Single median','Single no edge','Location','NorthWest');


