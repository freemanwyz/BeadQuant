function [ varOut, nPairs, meanOut ] = getVarNeib1( G,maskIn,varargin )
% GETVARNEIB1 Compute variance using neighboring pairs
% mean, meanCorrect, meanAvg

numvarargs = length(varargin);
optargs = {'mean',2};
optargs(1:numvarargs) = varargin;
[varMthd,divRatio] = optargs{:};

[Nx,Ny] = size(G);
d = zeros(Nx+2,Ny+2)+Inf;
xrg = 2:(Nx+1);
yrg = 2:(Ny+1);
d(xrg,yrg) = G;

% get pair differences -----
zx = zeros(length(maskIn),8);
xvec = [0,1,1,-1,0,-1,1,-1];
yvec = [1,0,1,1,-1,0,-1,-1];
for ii=1:8
    tmp0 = G - d(xrg+xvec(ii),yrg+yvec(ii));
    zx(:,ii) = tmp0(maskIn);
end

% compute variance -----
if strcmp(varMthd,'mean')
    noise = zx(:);
    noise = noise(~isinf(noise) & ~isnan(noise));
    n0 = length(noise);
    meanOut = mean(noise);
%     if meanOut>0
%         noise = noise(noise<=0);
%     else
%         noise = noise(noise>=0);
%     end
%     meanOut = length(noise)/n0;
    noise = noise - mean(noise);
    noise2 = noise.^2;
    nPairs = length(noise2);
    varOut = mean(noise2)/divRatio;
end
if strcmp(varMthd,'meanSingleDirection')
    noise = zx(:,2);
    noise = noise(~isinf(noise) & ~isnan(noise));
    meanOut = mean(noise);
    noise2 = noise.^2;
    nPairs = length(noise2);
    varOut = mean(noise2)/divRatio;
end
if strcmp(varMthd,'meanAvg')
    zx(isinf(zx) | isnan(zx)) = NaN;
    noise = nanmean(zx,2);
    noise = noise(~isinf(noise) & ~isnan(noise));
    meanOut = mean(noise);
%     if meanOut>0
%         noise = noise(noise<=0);
%     else
%         noise = noise(noise>=0);
%     end
    noise = noise - mean(noise);    
    noise2 = noise.^2;
    nPairs = length(noise2);
%     varOut = mean(noise2);
    varOut = mean(noise2)*8/9;
end
if strcmp(varMthd,'meanCorrect')
    noise = zx(:);
    noise = noise(~isinf(noise) & ~isnan(noise));
    meanOut = mean(noise);
    noise2 = noise.^2;
    NN = length(noise2);
    [~,nI] = sort(noise2,'descend');
    noise2 = noise2( nI((round(NN*0.1)+1):NN) );  % !!
    nPairs = length(noise2);
    varOut = mean(noise2)*1.6076/divRatio;  % !! 0.15:1.9214  0.1:1.6076
end

end




