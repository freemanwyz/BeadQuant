function Gt  = preTransform( Img, ImgCrop, p )
%PRETRANSFORM Transformation to stablize variance

Gorg = Img;
G = ImgCrop;

if strcmp(p.trans, 'na')
    Gt = G;
end
if strcmp(p.trans, 'log')
    Gt = log(G+1);
    Gt = Gt - min(Gt(:));
    Gt = Gt/max(Gt(:))*255;
end
if strcmp(p.trans, 'sqrt')
    Gt = sqrt(G);
    Gt = Gt - min(Gt(:));
    Gt = Gt/max(Gt(:))*255;
end
if strcmp(p.trans, 'vst2') || strcmp(p.trans, 'vst1')
    thr1 = 4:1:80;
    vx = Gorg*0;
    vxVec = zeros(1,length(thr1)-1);
    for jj=1:(length(thr1)-1)
        mask0 = find(Gorg>=thr1(jj) & Gorg<thr1(jj+1));
        vx0 = getVarNeib(Gorg,mask0,0,'mean',1);
        vxVec(jj) = vx0;
        vx(mask0) = vx0;
    end
%     figure;plot(thr1(1:(end-1)),vxVec);
%     xlabel('Intensity');
%     ylabel('Variance');
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
    Gt = zeros(size(G));
    if strcmp(p.trans, 'vst2')
        Gt(:) = log( yall - a + sqrt((yall-a).^2+c) );
    else
        Gt(:) = sqrt( yall - a + c);
    end
    
    vxVec = zeros(1,length(thr1)-1);
    for jj=1:(length(thr1)-1)
        mask0 = find(G>=thr1(jj) & G<thr1(jj+1));
        vx0 = getVarNeib(Gt,mask0,0,'mean',1);
        vxVec(jj) = vx0;
    end
%     figure;plot(thr1(1:(end-1)),vxVec);
%     xlabel('Intensity');
%     ylabel('Variance');
    
end

end

