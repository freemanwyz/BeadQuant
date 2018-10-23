function [myCurves,H0,psi,Ez] = beadExtractFANorm(sigch2, resBead, outPath, maxVal, vst, plotme)
% sigch2 is the (combined) signal channel
% USE delta F over F (DFF)

nFrame = size(sigch2,3);
myCurves = zeros(length(resBead),nFrame);
[Nx,Ny] = size(sigch2(:,:,1));

sig1 = reshape(sigch2,[],nFrame);
x0 = zeros(Nx,Ny);
x0a = zeros(Nx,Ny);
x0noise = zeros(Nx,Ny);

nBead = length(resBead);
x0idx = cell(nBead,1);
x0aidx = cell(nBead,1);
x0noiseidx = cell(nBead,1);

iirg = 1:nBead;
% parfor ii=iirg
for ii=iirg
    if mod(ii,100)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    d = double(tmp)/maxVal;
    if vst
        d = sqrt(d);
    end
    df0 = mean(d(:,1));
%     df0a = max(df0,0.01);
%     d1 = d;
    d1 = d - df0;
%     d1 = (d - df0)/df0a;
    
    % factor analysis or NMF or PCA -----
    [H0, psi, Ez, mu, llh] = fa(d1,1);
    H0(H0<0) = 0;   
    
    Ez1 = Ez - min(Ez);
%     Ez1 = (Ez - min(Ez))*df0a;
%     Ez1 = (Ez - min(Ez))*df0a*10 + df0;
    
%     H = 1./psi.*H0;
%     d2 = bsxfun(@times,d1,H);
%     Ez1 = sum(d2,1)./sum(H);
    
%     Ez1(Ez1<0) = 0;
    
    if vst
        Ez1 = Ez1.^2;
    end    
    myCurves(ii,:) = Ez1;  % from DFF to F
    
    m0 = H0./max(H0);
    x0idx{ii} = m0;
    
    H = 1./psi.*H0;
    m0a = H./max(H);
    x0aidx{ii} = m0a;
    
    x0noiseidx{ii} = psi;
end

% Ez1 = Ez - min(Ez);
% dRec = H0 * Ez1 * df0 + df0;
% figure;plot(dRec');ylim([0 1]);

for ii=iirg
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx Ny],tmpIdx(:,1),tmpIdx(:,2));
    x0(beadIdxEle) = x0idx{ii};
    x0a(beadIdxEle) = x0aidx{ii};
    x0noise(beadIdxEle) = x0noiseidx{ii};
end

sigch2Mean = zeros(Nx,Ny);
for ii=1:nFrame
    sigch2Mean = sigch2Mean + double(sigch2(:,:,ii))/maxVal;
end
sigch2Mean = sigch2Mean/nFrame;

if plotme
    K0 = cat(3,x0,x0*0,sigch2Mean*0);
    imwrite(K0,[outPath,filesep,'extracted_weighted_fa_dff_s.tif']);
    K0a = cat(3,x0a,x0*0,sigch2Mean*0);
    imwrite(K0a,[outPath,filesep,'extracted_weighted_fa_dff_se.tif']);
    K0noise = cat(3,x0noise,x0*0,sigch2Mean*0);
    imwrite(K0noise,[outPath,filesep,'extracted_weighted_fa_dff_e.tif']);
end


end








