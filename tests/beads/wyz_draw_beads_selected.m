%% plot selected beads

xI = [1 2 3 4 5];

load('resBeads.mat');
datAvg = double(imread('beadAvg16.tif'))/65535;
[Nx,Ny] = size(datAvg);

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
