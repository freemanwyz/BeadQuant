%% simulation for VST

p = getInitParam('preset','yinxue_nucleus');
fp0 = [p.tp filesep 'data' filesep 'dat_test' filesep 'C2-4-DS4-cherry-s100b-3.tif'];
[Img, ImgCrop, q] = loadData(fp0, 'tif', p);
p.trans = 'vst2';

% I = Img;
I = uint8(Img);
J = imnoise(I,'poisson');

%% VST
Gtrue = double(I);
Gorg = double(J);
G = double(J);

thr1 = 4:1:60;

vx = Gorg*0;
vxVec = zeros(1,length(thr1)-1);
for jj=1:(length(thr1)-1)
    mask0 = find(Gorg>=thr1(jj) & Gorg<thr1(jj+1));
    dif = double(Gorg(mask0) - Gtrue(mask0));
    vx0 = var(dif);
    vxVec(jj) = vx0;
    vx(mask0) = vx0;
end
figure;plot(thr1(1:(end-1)),vxVec);

vx1 = Gorg*0;
vxVec1 = zeros(1,length(thr1)-1);
for jj=1:(length(thr1)-1)
    mask0 = find(Gorg>=thr1(jj) & Gorg<thr1(jj+1));
    vx0 = getVarNeib(Gorg,mask0,0,'mean',1);
    vxVec1(jj) = vx0;
    vx1(mask0) = vx0;
end
figure;plot(thr1(1:(end-1)),vxVec1);

idxSel = find(vx>0);
y = vx(idxSel);
if strcmp(p.trans, 'vst2')
    x = Gorg(idxSel).^2;
else
    x = Gorg(idxSel);
end
corxy = corrcoef(x,y);
s2eta = corxy(1,2).*std(y)./std(x);
s2eps = var(y - s2eta.*x);
a = mean(Gorg(Gorg<4));
c = s2eps/s2eta;
yall = G(:);
yallRef = Gtrue(:);
Gt = zeros(size(G));
GtRef = zeros(size(G));
if strcmp(p.trans, 'vst2')
    Gt(:) = log( yall - a + sqrt((yall-a).^2+c) );
    GtRef(:) = log( yallRef - a + sqrt((yallRef-a).^2+c) );
else
    Gt(:) = sqrt( yall - a + c);
end

% vxVec = zeros(1,length(thr1)-1);
% for jj=1:(length(thr1)-1)
%     mask0 = find(G>=thr1(jj) & G<thr1(jj+1));
%     vx0 = getVarNeib(Gt,mask0,0,'mean',1);
%     vxVec(jj) = vx0;
% end
% figure;plot(thr1(1:(end-1)),vxVec);


vxVec = zeros(1,length(thr1)-1);
for jj=1:(length(thr1)-1)
    mask0 = find(Gorg>=thr1(jj) & Gorg<thr1(jj+1));
    dif = double(Gt(mask0) - GtRef(mask0));
    vx0 = var(dif);
    vxVec(jj) = vx0;
end
figure;plot(thr1(1:(end-1)),vxVec);


