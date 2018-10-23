function myCurves = beadExtractIter(sigch2, resBead, outPath, maxVal)
% sigch2 is the (combined) signal channel

nFrame = size(sigch2,3);
myCurves = zeros(length(resBead),nFrame);
[Nx,Ny] = size(sigch2(:,:,1));

sig1 = reshape(sigch2,[],nFrame);
x0 = zeros(Nx,Ny);

for ii=1:length(resBead)
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    d = double(tmp)/maxVal;
    
    W = mean(d,1);
    for tt=1:10
        kkInv = 1./sum(d.^2,2);
        kb = d*reshape(W,[],1);
        H = kkInv.*kb;
        d1 = d.*repmat(H,1,nFrame);
        rf1 = repmat(reshape(W,1,[]),length(beadIdxEle),1);
        dif1 = d1 - rf1;
%         err1 = mean(dif1.^2,2).^2;
        err1 = mean(dif1.^2,2);
        err1Weight = 1./err1;
        d2 = d.*repmat(err1Weight,1,nFrame);
        W = sum(d2,1)/sum(err1Weight); 
    end
    myCurves(ii,:) = W';   
    m1 = err1Weight./max(err1Weight);
    x0(beadIdxEle) = m1;   
end

sigch2Mean = zeros(Nx,Ny);
for ii=1:nFrame
    sigch2Mean = sigch2Mean + double(sigch2(:,:,ii))/maxVal;
end
sigch2Mean = sigch2Mean/nFrame;

K0 = cat(3,x0,x0*0,sigch2Mean*0.5);
imwrite(double(K0),[outPath,filesep,'extracted_weighted.tif']);

end








