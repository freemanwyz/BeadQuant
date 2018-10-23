%%
p.crop = 1;
p.trans = 'vst2'; % na, log, sqrt, vst1, vst2
p.test = 'order'; % t or order
p.correct = 'bonf'; % bonf
p.thrSig = 0.05;
p.thr0 = 10;
p.thr1 = 180;
p.thrg = 5;
p.thrZ = norminv(1-p.thrSig);
p.min_size = 4;  % minimun synapse size
p.max_ratio = 2;  % maximum synapse size (ratio to 1um^2)
p.conn_direc = 8;  % define neighbor, 8 to use corner
p.SE3 = strel('square',3);
tbl = load('./src/ostbl.mat');
p.mu = tbl.mu;
p.sigma = tbl.sigma;
[ flst, ~ ] = getFlst( './src/f1.txt' );

%% load data
q = [];
fin = ['./data/',flst{1},'.lsm'];
[I, Iarg] = lsmread(fin);
Nx = Iarg.dimX;
Ny = Iarg.dimY;
Isyn = reshape(I(1,2,1,:,:),[Nx,Ny]);
G0 = double(Isyn);
figure;hist(G0(:),256);

fin1 = strrep(fin,'weak','strong');
[I, Iarg] = lsmread(fin1);
Nx = Iarg.dimX;
Ny = Iarg.dimY;
Isyn = reshape(I(1,1,1,:,:),[Nx,Ny]);
G1 = double(Isyn);

kMask = G0>5 & G0<230 & G1<230;
G0Mean = mean(G0(kMask));
G1Mean = mean(G1(kMask));
ratioG1toG0 = G1Mean/G0Mean;

nPix = Nx*Ny;
noiseAdd = (rand(1,nPix)-0.5)*1e-6;
[G0s, I0] = sort(G0(:));
[G1s, I1] = sort(G1(:));

%% histogram matching -----
G1x = zeros(size(G1));
G0p1 = G0 + 1;
G1p1 = G1 + 1;
rtx = zeros(1,256);
for ii=1:256
    rt0 = mean(G1p1(G0p1==ii))/mean(G0p1(G0p1==ii));
    G1x(G0p1==ii) = G1p1(G0p1==ii)/rt0;
    rtx(ii) = rt0;
end
figure;plot(1:256,rtx);

%% weak as template
% var est multi -----
% G1a = G1./ratioG1toG0;
G0 = G0p1;
G1a = G1x;

Gdif = G0 - G1a;
% GdifSel = Gdif(kMask);

thr1 = 1:120;
vx = zeros(1,length(thr1));
mx = zeros(1,length(thr1));
nPixx = vx;
for ii=thr1
    x0 = Gdif(G0==ii);
    vx(ii) = mean(x0.^2);
    mx(ii) = mean(x0);
    nPixx(ii) = length(x0);
end
% figure;plot(thr1,vx);
figure;plot(thr1,mx);
figure;plot(thr1,log10(nPixx+1));

% var est single -----
vxb = zeros(1,length(thr1));
for jj=1:length(thr1)
    mask0 = find(G0==thr1(jj));
    vx0 = getVarNeib(G0,mask0,0,'mean');
    vxb(jj) = vx0;
end
% figure;plot(thr1,vxb);

% plot -----
figure;plot(thr1,[vx;vxb]);
ylabel('\sigma^2');
legend('two images','single image','Location','northwest');
title('Weak one as reference');
figure;plot(thr1,[log10(vx+1);log10(vxb+1)]);
ylabel('log(\sigma^2+1)');
legend('two images','single image','Location','northwest');
title('Weak one as reference');

%% strong as template
if 0
    % var est multi
    G0a = G0.*ratioG1toG0;
    G0a(G0a>255) = 255;
    
    Gdif = G1 - G0a;
    % GdifSel = Gdif(kMask);
    
    thr1 = 1:240;
    mx = zeros(1,length(thr1));
    for ii=thr1
        x0 = Gdif(G1==ii);
        vx(ii) = mean(x0.^2)/2;
        mx(ii) = mean(x0);
    end
    % figure;plot(thr1,vx);
    figure;plot(thr1,mx);
    
    % var est single
    vxb = zeros(1,length(thr1));
    for jj=1:length(thr1)
        mask0 = find(G1==thr1(jj));
        vx0 = getVarNeib(G1,mask0,0);
        vxb(jj) = vx0;
    end
    % figure;plot(thr1,vxb);
    
    % plot
    figure;plot(thr1,[vx;vxb]);
    ylabel('\sigma^2');
    legend('two images','single image','Location','northwest');
    title('Strong one as reference');
    figure;plot(thr1,[log10(vx+1);log10(vxb+1)]);
    ylabel('log(\sigma^2+1)');
    legend('two images','single image','Location','northwest');
    title('Strong one as reference');
end

