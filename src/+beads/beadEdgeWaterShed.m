function resBead0 = beadEdgeWaterShed( d1,resBead0,thr0 )
d1bw = 1*(d1>0);
d1bw = bwareaopen(d1bw,20);
d1 = d1.*d1bw;
d1edge = edge(d1,'Sobel',thr0);
d1edge = bwareaopen(d1edge, 2);
d1bw(d1edge)=0;
D = bwdist(~d1bw);
D(D>4) = 4;  % avoid two nearby local maximum for watershed
D = -D;
D(~d1bw) = -Inf;
L = watershed(D);
%             Lrgb = label2rgb(L,'jet',[.5 .5 .5]);
%             figure;imshow(Lrgb,'InitialMagnification','fit')
idxLabel = unique(L(:));
idxLabel = idxLabel(idxLabel>0);
synCnt = max(resBead0(:));
for tt = 1:length(idxLabel)
    idxSel = find(L==idxLabel(tt));
    meanSel = mean(d1(idxSel));
    if meanSel>0.1 && length(idxSel)>25
        tmp = zeros(size(d1));
        tmp(idxSel) = 1;
        statSel = regionprops(tmp,'Area','Perimeter');
        rtSel = statSel.Perimeter^2/statSel.Area/4/pi;
        if rtSel<=1
            synCnt = synCnt + 1;
            tmp = bwmorph(tmp,'spur');
            tmp = imdilate(tmp,strel('square',4));
            idxSel = tmp>0;
            resBead0(idxSel) = synCnt;
        end
    end
end

end

