function [xVar,xVarVar,xVarSd,nPix] = getVarCurveMulti(dats,meanRg)
% TODO: change output data type to table

datMean = mean(dats,2);
N1 = length(meanRg)-1;
xVar = zeros(1,N1);
xVarVar = zeros(1,N1);
xVarSd = zeros(1,N1);
nPix = zeros(1,N1);

for ii=1:N1
    fprintf('%d\n',ii);
    x1 = meanRg(ii);
    x2 = meanRg(ii+1);
    idx = find(datMean>=x1 & datMean<x2);
    nPix(ii) = length(idx);
    datSel = dats(idx,:);
    varSel = var(datSel,0,2);
    xVar(ii) = mean(varSel);
    xVarVar(ii) = var(varSel);
    xVarSd(ii) = std(varSel);    
end

end