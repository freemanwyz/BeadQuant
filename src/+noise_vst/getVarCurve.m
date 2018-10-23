function [xVar,nPix,nPairs,xMean] = getVarCurve(dat0,meanRg,usePara,mthd,divRatio)
% TODO: change output data type to table

N1 = length(meanRg)-1;
xVar = zeros(1,N1);
nPix = zeros(1,N1);
nPairs = zeros(1,N1);
xMean = zeros(1,N1);

if usePara
    parfor ii=1:N1
        fprintf('%d\n',ii);
        [ xVar0, nPix0, nPairs0, xMean0 ] = getVarCurvePoint(dat0,meanRg,mthd,divRatio,ii);
        nPix(ii) = nPix0;
        xVar(ii) = xVar0;
        nPairs(ii) = nPairs0;
        xMean(ii) = xMean0;
    end
else
    for ii=1:N1
        fprintf('%d\n',ii);
        [ xVar0, nPix0, nPairs0, xMean0 ] = getVarCurvePoint(dat0,meanRg,mthd,divRatio,ii);
        nPix(ii) = nPix0;
        xVar(ii) = xVar0;
        nPairs(ii) = nPairs0;
        xMean(ii) = xMean0;
    end
end
end

function [ xVar0, nPix0, nPairs0, xMean0 ] = getVarCurvePoint(dat0,meanRg,mthd,divRatio,ii)
x1 = meanRg(ii);
x2 = meanRg(ii+1);
idx = find(dat0>=x1 & dat0<x2);
nPix0 = length(idx);
if ii==110
%     keyboard
end
[ xVar0, nPairs0, xMean0 ] = getVarNeib1( dat0,idx,mthd,divRatio );
end
