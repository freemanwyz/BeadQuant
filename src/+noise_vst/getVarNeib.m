function [ varOut, nPairs ] = getVarNeib( G,maskIn,maskOnly,varargin )

numvarargs = length(varargin);
optargs = {'mean',2};
optargs(1:numvarargs) = varargin;
[varMthd,divRatio] = optargs{:};

[Nx,Ny] = size(G);
if maskOnly
    tmp = ones(Nx,Ny);
    tmp(maskIn) = 0;
    maskOut = find(tmp>0);
    if ~isempty(maskOut)
        G(maskOut) = Inf;
    end
end
d = zeros(Nx+2,Ny+2)+Inf;
xrg = 2:(Nx+1);
yrg = 2:(Ny+1);
d(xrg,yrg) = G;

l1 = G - d(xrg,yrg+1);
l2 = G - d(xrg+1,yrg);
l3 = G - d(xrg+1,yrg+1);
l4 = G - d(xrg-1,yrg+1);

l1a = G - d(xrg,yrg-1);
l2a = G - d(xrg-1,yrg);
l3a = G - d(xrg+1,yrg-1);
l4a = G - d(xrg-1,yrg-1);

noise = [l1(maskIn),l2(maskIn),l3(maskIn),l4(maskIn),l1a(maskIn),l2a(maskIn),l3a(maskIn),l4a(maskIn)];
noise = noise(~isinf(noise) & ~isnan(noise));

% noise = noise - mean(noise);  % !!
noise2 = noise.^2;
if strcmp(varMthd,'mean')
    NN = length(noise2);
    [~,nI] = sort(noise2,'descend');
    noise2 = noise2( nI((round(NN*0.1)+1):NN) );  % !!
    nPairs = length(noise2);
    varOut = mean(noise2)*1.6076/divRatio;  % !! 0.15:1.9214  0.1:1.6076
%     varOut = mean(noise2)/divRatio;
else
    nPairs = length(noise2);
    varOut = median(noise2)/divRatio;
end


end

