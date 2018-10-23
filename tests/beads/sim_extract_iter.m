%% beads
% load('D:\neuro_WORK\peptide_ucdavis_WORK\dump\20150810\stp3.session.mat');
sig = dd.sig;
[Nx,Ny,nFrame] = size(sig);
sig1 = reshape(sig,[],nFrame);
beadIdx = dd.resBead;
nBead = length(beadIdx);

x0 = zeros(size(dd.bgMean));
x0a = zeros(size(dd.bgMean));
xc = zeros(length(beadIdx),nFrame);

plotMe = 1;

% iirg = 2143;
% iirg = 157;
iirg = 11397;
% iirg = 1:nBead;
for ii=iirg
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = beadIdx{ii};
    beadIdxEle = sub2ind(size(dd.bgMean),tmpIdx(:,1),tmpIdx(:,2));
    tmp = sig1(beadIdxEle,:);
    d = double(tmp)/pp.maxVal;
    %     d = sqrt(d);
    
    W = mean(d,1);
    if plotMe
        figure
    end
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
        if plotMe
            plot(W);hold on
        end
    end
    if plotMe
        hold off
    end     

    xc(ii,:) = W';   
    m1 = err1Weight./max(err1Weight);
    x0(beadIdxEle) = m1;
end

if 1
    K0 = cat(3,x0,dd.bgMean*0.25,x0*0);
    figure;imshow(K0);
%     figure;plot(xc.');
end






