function resBead0 = beadCircleSeq( d1,resBead0,useDilate,baseMask)
done = 0;
iCnt = max(resBead0(:));
[Nx,Ny] = size(d1);
neibVec = [0, -1, 1, -Nx, Nx, -Nx+1, -Nx-1, Nx+1, Nx-1];
while ~done
    d1bw = 1*(d1>0);
    d1bw = bwareaopen(d1bw,20);
    d1 = d1.*d1bw;
    D = bwdist(~d1bw);
    if max(D(:))>2
        dist0LM = imregionalmax(D);
        cIdx = find(dist0LM>0);
        [xC,xI] = sort(D(cIdx));
        if length(xC)>1
            C = xC(2);
            I = xI(2);
        else
            C = xC(1);
            I = xI(1);
        end
        [ci,cj] = ind2sub(size(D),cIdx(I));
        [di,dj] = find(D > - 1);
        dist0 = sqrt((di - ci).^2 + (dj-cj).^2);
        if C>2 && sum(dist0<=C)>20 && C<9;
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
            resBead0(idx0) = iCnt;
            tmp = zeros(size(d1));
            tmp(idx0) = 1;
        else
            tmp = zeros(size(d1));
            tmp(dist0<=C) = 1;
        end
        if useDilate
            tmp = imdilate(tmp,strel('square',3));
        end        
        d1(tmp>0) = 0;
    else
        done = 1;
    end
end


end

