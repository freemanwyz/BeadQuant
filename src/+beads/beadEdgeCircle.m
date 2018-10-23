function resBead0 = beadEdgeCircle( d1,resBead0, thr0, baseMask )
[Nx,Ny] = size(d1);
neibVec = [0, -1, 1, -Nx, Nx, -Nx+1, -Nx-1, Nx+1, Nx-1];
d1bw = 1*(d1>0);
d1bw = bwareaopen(d1bw,20);
d1 = d1.*d1bw;
d1edge = edge(d1,'Sobel',thr0);
d1edge = bwareaopen(d1edge, 2);
d1bw(d1edge)=0;
D = bwdist(~d1bw);
if max(D(:))>2
    iCnt = max(resBead0(:));
    availMap = 1*(d1>0);
    dist0LM = imregionalmax(D);
    while sum(dist0LM(:))>0
        cIdx = find(dist0LM>0);
        [C,I] = max(D(cIdx));
        [ci,cj] = ind2sub(size(D),cIdx(I));
        [di,dj] = find(D > - 1);
        dist0 = sqrt((di - ci).^2 + (dj-cj).^2);
        nOut = sum(D(dist0<=C)<2); % avoid inner circles       
        if C>2 && sum(availMap(dist0<=C))>15 && nOut>5;
            idx0 = find(dist0<=C & dist0>(C-1));
            resxx = zeros(1,8);
            % alignment -----
            for tt = 1:8
                idx0a = idx0 + neibVec(tt);
                idx0a = idx0a(idx0a>0);
                idx0a = idx0a(idx0a<=Nx*Ny);
                resxx(tt) = sum(baseMask(idx0a));
            end         
            [~,I1] = max(resxx);
            idx0 = find(dist0<=C) + neibVec(I1);
            idx0 = idx0(idx0>0);
            idx0 = idx0(idx0<=Nx*Ny);
            iCnt = iCnt + 1;
%             idxOverlap = resBead0(:)>0 & dist1>0;  % overlap part set to 0
            resBead0(idx0) = iCnt;
%             resBead0(idxOverlap) = 0;
        end
        idx0 = find(dist0<=C/2);        
%         idx0 = idx0(idx0>0);
%         idx0 = idx0(idx0<=Nx*Ny);
        availMap(idx0) = 0;
        dist0LM(idx0) = 0;
    end
end

end

