%% beads
% load('D:\neuro_WORK\peptide_ucdavis_WORK\dump\20150811_ring\stp3.session.mat');
sig = dd.sig;
[Nx,Ny,nFrame] = size(sig);
sig1 = reshape(sig,[],nFrame);
beadIdx = dd.resBead;
nBead = length(beadIdx);

x0 = zeros(size(dd.bgMean));
x0a = zeros(size(dd.bgMean));
xc = zeros(length(beadIdx),nFrame);

clear x0idx x0aidx
x0idx = cell(nBead,1);
x0aidx = cell(nBead,1);

pp.maxVal = 65535;

% iirg = 3;
% iirg = 2143;
% iirg = 157;
% iirg = 11397;
iirg = 1:nBead;
% iirg = 1:10;

ii = iirg;
parfor (ii=iirg,10)
% for ii=iirg
    if mod(ii,10)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = beadIdx{ii};
    beadIdxEle = sub2ind(size(dd.bgMean),tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    dorg = double(tmp)/pp.maxVal;
    dorgf0 = dorg(:,1);
    dorgdf = bsxfun(@minus,dorg,dorgf0);
    dorgdff = bsxfun(@rdivide,dorgdf,dorgf0+0.01);   
    
    d = sqrt(dorg);
    df0 = d(:,1);
    d = bsxfun(@minus,d,df0);
    d = bsxfun(@rdivide,d,df0+0.01);

    % factor analysis or NMF or PCA -----
    [H0, psi, Ez, mu, llh] = fa(d,1);
    H0(H0<0) = 0;
    H = 1./psi.*H0;
    
    dt1 = dorgdff'.*repmat(reshape(H,1,[]),nFrame,1);
    xc(ii,:) = sum(dt1,2)/sum(H);

    m0 = H0./max(H0);
    m0a = H./max(H);
    x0idx{ii} = m0;
    x0aidx{ii} = m0a;
end

for ii=iirg
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = beadIdx{ii};
    beadIdxEle = sub2ind(size(dd.bgMean),tmpIdx(:,1),tmpIdx(:,2));
    x0(beadIdxEle) = x0idx{ii};
    x0a(beadIdxEle) = x0aidx{ii};
end

if 1
    K0 = cat(3,x0,dd.bgMean*0.25,x0*0);
    figure;imshow(K0);title('Weight')
    K0a = cat(3,x0a,dd.bgMean*0.25,x0a*0);
    figure;imshow(K0a);title('Weight/error');
    % K0a = cat(3,x0a,dd.bgMean*0.25,x0a*0);
    % figure;imshow(K0a);    
end

if 1
    figure;plot(xc.');
else
    figure;plot(dorgdff');
    figure;plot(d');
    xc1 = xc(ii,:);
    figure;plot(xc1');
end





