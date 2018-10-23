%% simulation setup
rng('default');
rng(88);
nFrame = 33;
trg = 0:(nFrame-1);
a = 0.7;
tau = 3;
xc = a*(1-exp(-trg/tau));
F0 = 0.05;  % 0.05
K = 20;
Nx = 41;
Ny = 41;
pattenType = 0;

%% init
A = zeros(Nx,Ny);
D = A;  % distance to the image center
E = A;  % distance to the light source
[ax,ay] = find(A==0);

distBead = sqrt((ax-21).^2 + (ay-21).^2);
D(sub2ind(size(A),ax,ay)) = distBead;
A(D<=9) = 1;
[rx,ry] = find(A>0);
resBead = {[rx,ry]};

if pattenType==0
    distDye = sqrt((ax-8).^2 + (ay-8).^2) + 1;
    E(sub2ind(size(A),ax,ay)) = distDye;
    A = A./abs(E-12.01)*2;
    A(A>1) = 1;
%     A(A<0.4) = 0.05;
else
    A = A.*(D.^2./60);
    A(D>8 & D<=9) = 0.5;
%     A(D<7) = 0.01;
    A(A>1) = 1;
%     A(D<7) = 0.6;
%     A(D<5) = 0.3;
%     A(D<2) = 0.1;
end
figure;imshow(A)

X = zeros(Nx,Ny,nFrame);
for ii=1:nFrame
%     X0 = A*xc(ii);
    X0 = A*xc(ii)*F0*K + F0;
    X0(A==0) = 0;
    X(:,:,ii) = X0;
end
d0 = A(A>0) * xc + F0;
% figure;plot(d0');ylim([0,1])

%% detection
serg = flip([0.01,0.02,0.05,0.1,0.2,0.3,0.5]);
% serg = 0.001;
Nsnr = length(serg);
err_snr = zeros(Nsnr,5);
Nsim = 1000;
for ss=1:Nsnr
    se0 = serg(ss);
    fprintf('SE is %f\n',se0);
    err = zeros(Nsim,5);
%     for ii=1:Nsim
    parfor ii=1:Nsim
        %     fprintf('%d\n',ii);
        Y = X + randn(size(X))*se0;
%         Y(Y>1) = 1;
%         Y(Y<0) = 0;
        % figure;plot(reshape(Y(21,21,:),[],1));
        
        % FA -----
        m1 = beadExtractFANorm(Y, resBead, './output/', 1, 0, 0);
%         m1 = m1-min(m1);
        m1norm = m1/mean(m1)*mean(xc);
        err1 = mean((m1norm - xc).^2);
        
%         if err1>0.02
%             keyboard
%         end
        
        % Mean --
        tmpIdx = resBead{1};
        beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
        sig1 = reshape(Y,[],nFrame);
        sig1Sel = sig1(beadIdxEle,:);
        m2 = mean(sig1Sel-F0,1);
%         m2 = m2 - min(m2);
        m2norm = m2/mean(m2)*mean(xc);
        err2 = mean((m2norm - xc).^2);
        
        % Oracle --
        wt1 = reshape(A,[],1);
        wt1 = wt1(beadIdxEle);
        coef = 1/(wt1'*wt1)*wt1';
        m3 = coef*(sig1Sel-F0);        
        
%         sig1Wt = bsxfun(@times,sig1Sel-F0,wt1);
%         m3 = mean(sig1Wt,1);
% %         m3 = m3 - F0;
% %         m3 = m3 - min(m3);
        m3norm = m3/mean(m3)*mean(xc);
        err3 = mean((m3norm - xc).^2);
        
        % Ring --
        m4 = beadExtract(Y, [21 21], 9, 5, 8, './output/', 1, 0);
        m4 = m4 - F0;
%         m4 = m4-min(m4);
        m4norm = m4/mean(m4)*mean(xc);
        err4 = mean((m4norm - xc).^2);
        
        % NMF --
        [W,m5] = nnmf(sig1Sel-F0,1);
%         m5 = m5-F0;
%         m5 = m5-min(m5);
        m5norm = m5/mean(m5)*mean(xc);
        err5 = mean((m5norm - xc).^2);
        
%         W1 = W./mean(W).*mean(wt1);
%         var(wt1-W1)
        
        err(ii,:) = [err1,err2,err3,err4,err5];
    end
    err_snr(ss,:) = mean(err,1);
%     err_snr(ss,:) = median(err,1);
end

figure;
semilogy(-10*log10(serg.^2),err_snr);
legend('factor analysis','mean','oracle','ring','NMF','Location','northeast');
xlabel('-10*log10(noise variance)');
ylabel('MSE');

% erra = err(:,1);
% erra = erra(erra<1e-3);
% fprintf('%f\n',mean(erra))
% fprintf('%f\n',mean(err(:,2)))
% fprintf('%f\n',mean(err(:,3)))
% fprintf('%f\n',mean(err(:,4)))

% % display --
% figure;plot(sig1')
% figure;plot(xc);hold on
% plot(m1norm);
% plot(m2norm);
% plot(m3norm);
% plot(m4norm);
% legend('ground truth','factor analysis','mean','oracle','ring','Location','southeast');
% hold off







