%% beads
% load('D:\neuro_WORK\peptide_ucdavis_WORK\dump\20150811_ring\stp3.session.mat');
sig = double(dd.sig)/pp.maxVal;
sigMean = mean(sig,3);
[Nx,Ny,nFrame] = size(sig);
sig1 = reshape(sig,[],nFrame);
beadIdx = dd.resBead;

% mc = dd.myCurves;
% tv = dd.time_vect;
% a = dd.dff_fit;
% tau = dd.tau_fit;
% idxSel = a<=2;
% % idxSel = a > 0.66 & a <= 0.7;
% idxSelNum = find(idxSel>0);
% beadIdx = beadIdx(idxSel);
% mc = mc(idxSel,:);
% plot(transpose(mc(1:1000,:)))

nBead = length(beadIdx);

% idx0 = 21;
% xcSel = mc(idx0:(idx0+10-1),:);
% plot(xcSel');
% tmpCell = num2cell( (1:10) + idx0 - 1);
% mystr = cellfun(@num2str,tmpCell,'UniformOutput',0);
% legend(mystr);

%% extraction
x0 = zeros(Nx,Ny);
x0a = zeros(Nx,Ny);
xc = zeros(length(beadIdx),nFrame);
xcAll = zeros(length(beadIdx),nFrame);

% iirg = 1:nBead;
iirg = 1481;
for ii=iirg
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = beadIdx{ii};
    beadIdxEle = sub2ind(size(dd.bgMean),tmpIdx(:,1),tmpIdx(:,2));
    sigEle = sig1(beadIdxEle,:);
    
    %     sigEleMax = max(sigEle,[],2);
    %     sigEleSel = sigEleMax<0.99;
    %     if sum(sigEleSel)<10
    %         continue
    %     end
    %     sigEle = sigEle(sigEleSel,:);
    %     beadIdxEle = beadIdxEle(sigEleSel);
    
    sigMeanEle = sigMean(beadIdxEle);
    level = graythresh(sigMeanEle);
    beadIdxEleSel = beadIdxEle(sigMeanEle>=level);
    x0(beadIdxEleSel) = 1;
    
    sigEleSel = sigEle(sigMeanEle>=level,:);
    F0 = sigEleSel(:,1);
    sigEleSelDf = bsxfun(@minus,sigEleSel,F0);
    sigEleSelDff = bsxfun(@rdivide,sigEleSelDf,max(F0,0.005));
    xc(ii,:) = mean(sigEleSelDff,1);
    %     xc(ii,:) = median(sigEleSelDff,1);
    
    F0 = sigEle(:,1);
    sigEleSelDf = bsxfun(@minus,sigEle,F0);
    sigEleSelDff1 = bsxfun(@rdivide,sigEleSelDf,max(F0,0.005));
    xcAll(ii,:) = mean(sigEleSelDff1,1);
    %     xcAll(ii,:) = median(sigEleSelDff,1);
    
    t1 = mean(sigEleSelDff1,2);
    x0a(beadIdxEle) = t1/max(t1);
    
end

% plot(sigEleSelDff');xlabel('Time point');ylabel('\DeltaF/F');title('Otsu');
% plot(sigEleSelDff1');xlabel('Time point');ylabel('\DeltaF/F');title('All pixels');

% figure;plot(xcAll(ii,:));hold on;plot(xc(ii,:));
% xlabel('Time point');ylabel('\DeltaF/F');
% legend('All pixels','Otsu method','Location','southeast');hold off

% plot(xc(ii,:));xlabel('Time point');ylabel('\DeltaF/F')

%% check bad curves
idx0 = 1101;
xcSel = xcAll(idx0:(idx0+100-1),:);
plot(xcSel');
% tmpCell = num2cell( (1:10) + idx0 - 1);
% mystr = cellfun(@num2str,tmpCell,'UniformOutput',0);
% legend(mystr);

%%
[B I] = sort(abs(xcAll(:,27) - 19.84));
% figure;plot(xcAll(ii,:))
% figure;plot(xc(ii,:))

figure;plot(xcAll(1481,:));
% figure;plot(transpose(sigEleSel(1:10,:)))

%% plot
K0 = cat(3,x0*0,zeros(Nx,Ny),sigMean);
figure;imshow(K0);

K1 = cat(3,x0a,zeros(Nx,Ny),sigMean);
figure;imshow(K1);

bg0 = double(dd.sig(:,:,1))/pp.maxVal;
K2 = cat(3,bg0*0,bg0*0,bg0*4);
figure;imshow(K2);

figure;imshow(dd.bgMean);

xcMax = max(xcAll,[],2);
% xcTop = xc(xcMax>10,:);
xcAllTop = xcAll(xcMax>20,:);
figure;plot(xcAllTop')

idxBad = find(xcAll(:,28)>133);

% [Nbead,Ntps] = size(xcTop);
% Ntest = 20;
% idxx = randperm(Nbead,Ntest);
% for ii=1:Ntest
%     figure
%     plot(xcTop(idxx(ii),:))
%     hold on
%     plot(xcAllTop(idxx(ii),:))
%     hold off
%     legend('Otsu','Mean')
% end





