%% preprocessing
% outDir = 'mtsyn_trim_small';
% mkdir('dump',outDir);

%% load data
[ flst, ~ ] = getFlst( './src/f1.txt' );
ii = 1;
fin = ['./data/',flst{ii},'.lsm'];
fprintf('%s  >>>>>>>>>>\n',fin);
[I, Iarg] = lsmread(fin);
Nx = Iarg.dimX;
Ny = Iarg.dimY;
neibVec = [0, -1, 1, -Nx, Nx];
Lx = Iarg.voxSizeX*1e6;
Lx2 = Lx^2;
Pix_per_synapse = round(1/Lx2);
Isyn = I(1,2,1,:,:);
Isyn = reshape(Isyn,[Nx,Ny]);
G = double(Isyn);
% Gt = G.^0.5;
% Gt = log(G+2);
% lambda = 0.1;
% Gt = (G.^lambda-1)/lambda;
Gt = G;
figure;hist(Gt(:));

% idx  = find(G<0.05);
% figure;hist(Gt(idx));
% figure;hist(G(idx));

% xx = 1:256;
% figure;plot(xx.^0.5);
% figure;plot(log(xx+1));

% K = Gt;
% K(K<=-5) = -5;
% K = K + 5;
% K = K/max(K(:));
% imshow(K);
figure;imshow(Gt/max(Gt(:)));

%% test
[Nx,Ny] = size(G);
d = zeros(Nx+2,Ny+2)+Inf;
xrg = 2:(Nx+1);
yrg = 2:(Ny+1);
d(xrg,yrg) = G;

l1 = G - d(xrg,yrg+1);
l2 = G - d(xrg+1,yrg);
l3 = G - d(xrg+1,yrg+1);
l4 = G - d(xrg-1,yrg+1);

thr = 0:1:180;
thr1 = [thr,256];
ii = 30;
maskIn = find(G>=thr1(ii) & G<thr1(ii+1));
noise = [l1(maskIn),l2(maskIn),l3(maskIn),l4(maskIn)];
noise = noise(~isinf(noise));
noise2 = (noise.^2)/2;
nPairs = length(noise);
noiseSel = noise2;
varOutMean = mean(noiseSel);
varOutMode = median(noiseSel);

a = unique(noise);
outa = [a,histc(noise(:),a)];
figure;hist(noise(:),a)

b = unique(noise2);
outb = [b,histc(noise2(:),b)];
% stem(outb(:,1),outb(:,2));
figure;hist(noise2(:),b)

%% compute noise level
maskOnly = 0;
% maskIn = 1:(Nx*Ny);
% maskIn0 = (find(G<0.05)).';
% maskIn1 = (find(G>=0.05)).';
% % maskIn1 = (find(G>=0.7 & G<0.75)).';
% v = getVarNeib(Gt,maskIn,maskOnly);
% v0 = getVarNeib(Gt,maskIn0,maskOnly);
% v1 = getVarNeib(Gt,maskIn1,maskOnly);

thr = 0:5:180;
thr1 = [thr,256];
vx = thr*0;
for ii=1:length(thr)
    mask0 = find(G>=thr1(ii) & G<thr1(ii+1));
%     mask0 = find(G>=thr1(ii));
    vx(ii) = getVarNeib(Gt,mask0,maskOnly);
end
figure;plot(thr,vx);

% %% analysis
% % histogram -----
% NN = 20;
% [~,centers] = hist(G(:),NN-1);
% centers = [min(G(:)),centers,max(G(:))];
% intensityMean = zeros(1,NN);
% noiseVar = zeros(1,NN);
% nPairs = zeros(1,NN);
% for jj=1:NN
%     idxx = (find(G>=centers(jj) & G<centers(jj+1) )).';
%     [noiseVar(jj),nPairs(jj)] = getVarNeib(G,idxx);
%     intensityMean(jj) = median(G(idxx));
% end
% figure;plot(intensityMean,noiseVar);

% % quantile -----
% centers = quantile(G(:),9);
% centers = [min(G(:)),centers,max(G(:))];
% intensityMean = zeros(1,10);
% noiseVar = zeros(1,10);
% for jj=1:10
%     idxx = (find(G>=centers(jj) & G<centers(jj+1) )).';
%     noiseVar(jj) = getVarNeib(G,idxx);
%     intensityMean(jj) = median(G(idxx));
% end
% figure;plot(intensityMean,noiseVar);

%% direct correction
thr = 0:5:180;
thr1 = [thr,256];
Gx = G;
for ii=1:length(thr)
    mask0 = find(G>=thr1(ii) & G<thr1(ii+1));
    if ii<100
        Gx(mask0) = G(mask0)/sqrt(vx(ii));
    else
        Gx(mask0) = G(mask0)/sqrt(vx(20));
    end
end
Gx1 = Gx./max(Gx(:));
figure;imshow(Gx1);


thr = 0:5:180;
thr1 = [thr,256];
vx1 = thr*0;
for ii=1:length(thr)
    mask0 = find(G>=thr1(ii) & G<thr1(ii+1));
%     mask0 = find(G>=thr1(ii));
    vx1(ii) = getVarNeib(Gx,mask0,maskOnly);
end
figure;plot(thr,vx1);










