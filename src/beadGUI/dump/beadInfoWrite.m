function beadInfoWrite( dd, outPath )
%BEADINFOWRITE Get information of beads and write them to a Excel file

goodBeadIdx = dd.goodBeadIdx;
controlBeadIdx = dd.controlBeadIdx;

nSel = length(goodBeadIdx);
nCtrl = length(controlBeadIdx);
resBead = dd.resBead;
nTot = length(resBead);

% last row is whole image mean
Idx = cell(nSel+nCtrl+1,1);
Row = cell(nSel+nCtrl+1,1);
Col = cell(nSel+nCtrl+1,1);
DeltaF = cell(nSel+nCtrl+1,1);
Tau = cell(nSel+nCtrl+1,1);
Typ0 = cell(nSel+nCtrl+1,1);

h = msgbox('Saving results');
allIdx = [goodBeadIdx,controlBeadIdx];
for ii=1:(nSel+nCtrl)
    idx0 = allIdx(ii);
    Idx{ii} = idx0;
    DeltaF{ii} = dd.dff_fit(idx0);
    Tau{ii} = dd.tau_fit(idx0);
    pix0 = resBead{idx0};
    minx = min(pix0(:,1));
    maxx = max(pix0(:,1));
    miny = min(pix0(:,2));
    maxy = max(pix0(:,2));    
    Row{ii} = round(mean([minx,maxx]));
    Col{ii} = round(mean([miny,maxy]));
    if ii<=nSel
        Typ0{ii} = 'Selected';
    else
        Typ0{ii} = 'Control';
    end
end

Idx{nSel+nCtrl+1} = 'Image Mean';
DeltaF{nSel+nCtrl+1} = dd.dff_fit_mean;
Tau{nSel+nCtrl+1} = dd.tau_fit_mean;
A = [Idx,Row,Col,DeltaF,Tau,Typ0];
cN0 = {'Idx','Row','Col','Delta F','Tau','Type'};
A = [cN0;A];

% T = table(Idx,Row,Col,DeltaF,Tau);
fname = [outPath,filesep,'beadsInfo.xls'];
% writetable(T,fname,'sheet',1);
% xlswrite(fname,A);
xlwrite(fname,A,'Beads');

ncnt = {'Summary'};
xlwrite(fname,ncnt,'Beads','H2');
ncnt = {'Control beads',nCtrl};
xlwrite(fname,ncnt,'Beads','H4');
ncnt = {'Beads Selected',nSel};
xlwrite(fname,ncnt,'Beads','H3');
ncnt = {'Beads detected',nTot};
xlwrite(fname,ncnt,'Beads','H5');

close(h);

end

