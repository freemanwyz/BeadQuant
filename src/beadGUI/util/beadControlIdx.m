function controlBeadIdx = beadControlIdx(resBead,sigch1Mean,idx0)
% idx0 = evalin('base','controlBeadCenters');
[Nx,Ny] = size(sigch1Mean);
nB = length(resBead);
kMap = zeros(Nx,Ny);
for ii=1:nB
    pix0 = resBead{ii};
    pix0L = sub2ind(size(sigch1Mean),pix0(:,1),pix0(:,2));
    kMap(pix0L) = ii;
end

nCtrl = size(idx0,1);
controlBeadIdx = zeros(1,nCtrl);
for ii=1:nCtrl
    controlBeadIdx(ii) = kMap(idx0(ii,1),idx0(ii,2));
end

controlBeadIdx = controlBeadIdx(controlBeadIdx>0);  % control beads might be drifted

end