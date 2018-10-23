function res = beadEval( d1org,baseMask,resBead0 )
K0a = cat(3,resBead0>0,d1org,baseMask);
resBead0x = resBead0>0;
resMask0 = d1org*0;
B = bwboundaries(resBead0x);
for nn = 1:length(B)
    B1 = B{nn};
    idxB = sub2ind(size(d1org),B1(:,1),B1(:,2));
    resMask0(idxB) = 1;
end
cc = resMask0>0 & baseMask;
rt0 = sum(cc(:));
%     rt0 = sum(resBead0>0)/sum(d1org>0);
K0 = cat(3,resMask0,d1org,baseMask);
res.fig1 = K0;
res.fig2 = K0a;
res.ratio = rt0;
res.beads = resBead0;
res.nBeads = length(unique(resBead0(:)))-1;
end

