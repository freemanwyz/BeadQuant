%% simulation
if 0
    x0 = -0.5:0.05:0.45;
    X = [repmat(x0,50,1);zeros(50,20)];
    Xn = X + randn(100,20)*0.1;
    % [Lambda,Psi,T,stats,F] = factoran(Xn,1);
    Xnt = transpose(Xn);
    [coeff,score,latent,tsquared,explained,mu] = pca(Xnt,'NumComponents',1);
end

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

x0idx = cell(nBead,1);
x0aidx = cell(nBead,1);

pp.maxVal = 65535;

iirg = 3;
% iirg = 2143;
% iirg = 157;
% iirg = 11397;
% iirg = 1:nBead;
% iirg = 1:10;

ii = iirg;
% parfor (ii=iirg,10)
for ii=iirg
    fprintf('%d\n',ii);
%     if mod(ii,1000)==0
%         fprintf('%d\n',ii);
%     end
    tmpIdx = beadIdx{ii};
    beadIdxEle = sub2ind(size(dd.bgMean),tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    d = double(tmp)/pp.maxVal;    
    
    % saturated points can introduce problems -----
%     idxdis1 = find(d>0.95);
%     d1pert = randn(length(idxdis1),1)*0.1;
%     d(idxdis1) = d(idxdis1) + d1pert;
    
    %     d = sqrt(d);
    dt = d.';

    % factor analysis or NMF or PCA -----
    [H0, psi, Ez, mu, llh] = fa(d,1);
    H0(H0<0) = 0;
    H = 1./psi.*H0;
    
    %     sigEle = mean(d,1);
    %     % figure;plot(sigEle);
    %
    %     % PCA
    %     [coeff,score,latent,tsquared,explained,mu] = pca(dt,'NumComponents',1);
    %     % figure;plot(coeff);
    %     % figure;plot(score);
    %
    %     % figure;plot(d(5,:));
    %     % figure;plot(d(19,:));
    %     % figure;plot(d(112,:));
    
%     % NMF
%     [Ez,W] = nnmf(dt,1);
%     % figure;plot(W);
%     % figure;plot(H);
%     % figure;plot(sort(H));
%     
%     % weight based on error
%     dt1 = dt./repmat(H,nFrame,1);
%     rf1 = repmat(W,1,length(beadIdxEle));
%     dif1 = dt1 - rf1;
%     err1 = mean(dif1.^2,1);
%     err1inv = 1./err1;
%     dt2 = dt.*repmat(err1inv,nFrame,1);
%     
%     xc(ii,:) = sum(dt2,2)/sum(err1inv);   
    % xc(ii,:) = sum(dt1,2)/sum(H);
    
    dt1 = dt.*repmat(reshape(H,1,[]),nFrame,1);
    xc(ii,:) = sum(dt1,2)/sum(H);

    % t1 = tmpIdx(H>0.08,:);
    % m1 = H(1,:)/max(H(1,:));
%     m1 = err1inv./max(err1inv);
    %     m1a = H(2,:)/max(H(2,:));
    % x0(sub2ind(size(dd.bgMean),t1(:,1),t1(:,2))) = 1;

    m0 = H0./max(H0);
    m0a = H./max(H);
    x0idx{ii} = m0;
    x0aidx{ii} = m0a;
    
%     x0(beadIdxEle) = m0;
%     x0a(beadIdxEle) = m0a;
    %     x0a(beadIdxEle) = m1a;
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

if 0
    figure;plot(xc.');
end





