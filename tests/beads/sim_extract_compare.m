% compare ring and fa
% with DFF

fdump = 'D:\neuro_WORK\peptide_ucdavis_WORK\dump\';
resRing = load([fdump,'\20150815_ring_dff\stp2.session.mat']);
resFA = load([fdump,'\20150818_otsu_dff_p1\stp2.session.mat']);
% resFA = load([fdump,'\20150818_otsu_dff\stp3.session.mat']);

curRing0 = resRing.dd.dff_mat;
curFA0 = resFA.dd.dff_mat;

%% plot 
sel0 = max(curRing0,[],1)>10;
curRing = transpose(curRing0(:,sel0));
curFA = transpose(curFA0(:,sel0));

[Nbead,Ntps] = size(curRing);
Ntest = 20;
idxx = randperm(Nbead,Ntest);
for ii=1:Ntest
    figure
    plot(curRing(idxx(ii),:))
    hold on
    plot(curFA(idxx(ii),:))
    hold off
    legend('Ring','Otsu')
end

dif0 = mean(curFA,2) - mean(curRing,2);
figure;hist(dif0,50)
