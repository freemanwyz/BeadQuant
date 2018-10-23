function myCurves = beadExtractFA(sigch2, resBead, outPath, maxVal)
% sigch2 is the (combined) signal channel
% USE delta F over F (DFF)

nFrame = size(sigch2,3);
myCurves = zeros(length(resBead),nFrame);
[Nx,Ny] = size(sigch2(:,:,1));

sig1 = reshape(sigch2,[],nFrame);
x0 = zeros(Nx,Ny);

nBead = length(resBead);
x0idx = cell(nBead,1);
% x0aidx = cell(nBead,1);

iirg = 1:nBead;
parfor ii=iirg
    % for ii=iirg
    if mod(ii,100)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    dorg = double(tmp)/maxVal;
    dorgf0 = dorg(:,1);
    dorgdf = bsxfun(@minus,dorg,dorgf0);
    dorgdff = bsxfun(@rdivide,dorgdf,max(dorgf0, 0.01));
    
    d = sqrt(dorg);
    df0 = d(:,1);
    d = bsxfun(@minus,d,df0);
    d = bsxfun(@rdivide,d,df0+0.01);
    
    % factor analysis or NMF or PCA -----
    [H0, psi] = fa(d,1);
    %     [H0, psi, Ez, mu, llh] = fa(d,1);
    H0(H0<0) = 0;
    H = 1./psi.*H0;
    
    dt1 = dorgdff'.*repmat(reshape(H,1,[]),nFrame,1);
    resdff = sum(dt1,2)/sum(H);
    myCurves(ii,:) = resdff*(max(dorgf0(1), 0.01)) + dorgf0(1);  % from DFF to F
    
    m0 = H0./max(H0);
    x0idx{ii} = m0;
    %     m0a = H./max(H);
    %     x0aidx{ii} = m0a;
end

for ii=iirg
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};
    beadIdxEle = sub2ind([Nx Ny],tmpIdx(:,1),tmpIdx(:,2));
    x0(beadIdxEle) = x0idx{ii};
    %     x0a(beadIdxEle) = x0aidx{ii};
end

sigch2Mean = zeros(Nx,Ny);
for ii=1:nFrame
    sigch2Mean = sigch2Mean + double(sigch2(:,:,ii))/maxVal;
end
sigch2Mean = sigch2Mean/nFrame;

K0 = cat(3,x0,x0*0,sigch2Mean*0.5);
imwrite(K0,[outPath,filesep,'extracted_weighted_fa_dff.tif']);

end








