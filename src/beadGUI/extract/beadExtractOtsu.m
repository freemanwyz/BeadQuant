function myCurves = beadExtractOtsu(sigch2, resBead, outPath, maxVal)
% sigch2 is the (combined) signal channel
% USE delta F over F (DFF)

nFrame = size(sigch2,3);
[Nx,Ny] = size(sigch2(:,:,1));
nBead = length(resBead);

sigch2 = double(sigch2)/maxVal;
sig1 = reshape(sigch2,[],nFrame);
sigMean = mean(sigch2,3);

x0 = zeros(Nx,Ny);
myCurves = zeros(length(resBead),nFrame);
iirg = 1:nBead;
USE_OTSU = 0;
USE_SAT = 1;
for ii=iirg
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    tmpIdx = resBead{ii};    
    beadIdxEle = sub2ind([Nx,Ny],tmpIdx(:,1),tmpIdx(:,2));
    sigEle = sig1(beadIdxEle,:);
    
    % !! remove saturated pixels, test only -----
    if ~USE_SAT
        sigEleMax = max(sigEle,[],2);
        idxSel = sigEleMax<0.999;
        if sum(idxSel)<10
            continue
        end
        beadIdxEle = beadIdxEle(idxSel);
        sigEle = sigEle(idxSel,:);
    end
    
    % curve extraction -----
    if USE_OTSU  % mean of high pixels
        sigMeanEle = sigMean(beadIdxEle);  % !!
        %     sigMeanEle = sqrt(sigMean(beadIdxEle));  % !!
        level = graythresh(sigMeanEle);
        beadIdxEleSel = beadIdxEle(sigMeanEle>=level);
        x0(beadIdxEleSel) = 1;
        sigEleSel = sigEle(sigMeanEle>=level,:);
        myCurves(ii,:) = mean(sigEleSel,1);
    else  % mean of all pixels
        myCurves(ii,:) = mean(sigEle,1);
    end

end

%% plot
% K0 = cat(3,x0*0.4,zeros(Nx,Ny),sigMean);
% imwrite(double(K0),[outPath,filesep,'extracted_mean.tif']);

end








