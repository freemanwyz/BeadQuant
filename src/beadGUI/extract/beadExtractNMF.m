function myCurves = beadExtractNMF(sigch2, resBead, outPath, maxVal)
% sigch2 is the (combined) signal channel

nFrame = size(sigch2,3);
myCurves = zeros(length(resBead),nFrame);
[Nx,Ny] = size(sigch2(:,:,1));

sig1 = reshape(sigch2,[],nFrame);
x0 = zeros(Nx,Ny);

USE_WEIGHT = 1;
USE_MNF = 0;

for ii=1:length(resBead)
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    d = double(tmp)/maxVal;
    if USE_NMF
        dt = d.';
        [W,H] = nnmf(dt,1);
    else
        [H0, psi, W, mu, llh] = fa(d,1);
        H = 1./psi.*H0;
    end    
    
    if USE_WEIGHT
        % weight based on scale
        dt1 = dt.*repmat(H,nFrame,1);
        myCurves(ii,:) = sum(dt1,2)/sum(H);
        m1 = H(1,:)/max(H(1,:));
        x0(beadIdxEle) = m1;
    else
        % weight based on error
        dt1 = dt./repmat(H,nFrame,1);
        rf1 = repmat(W,1,length(beadIdxEle));
        dif1 = dt1 - rf1;
        err1 = mean(dif1.^2,1);
        err1inv = 1./err1;
        dt2 = dt.*repmat(err1inv,nFrame,1);
        myCurves(ii,:) = sum(dt2,2)/sum(err1inv);
        m1 = err1inv./max(err1inv);
        x0(beadIdxEle) = m1;
    end   
    
end

sigch2Mean = zeros(Nx,Ny);
for ii=1:nFrame
    sigch2Mean = sigch2Mean + double(sigch2(:,:,ii))/maxVal;
end
sigch2Mean = sigch2Mean/nFrame;

K0 = cat(3,x0,sigch2Mean*0.5,x0*0);
imwrite(K0,[outPath,filesep,'extracted_weighted.tif']);

end








