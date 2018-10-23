G = randn(1024);
Gt = G;

maskOnly = 0;
maskIn = 1:(1024*1024);
[Nx,Ny] = size(G);
d = zeros(Nx+2,Ny+2)+Inf;
xrg = 2:(Nx+1);
yrg = 2:(Ny+1);
d(xrg,yrg) = G;

l1 = G - d(xrg,yrg+1);
l2 = G - d(xrg+1,yrg);
l3 = G - d(xrg+1,yrg+1);
l4 = G - d(xrg-1,yrg+1);

noise = [l1(maskIn),l2(maskIn),l3(maskIn),l4(maskIn)];
noise = noise(~isinf(noise));
noise2 = (noise.^2)/2;
nPairs = length(noise2);

[ns,nsI] = sort(noise2);
% nq1 = ns(round(nPairs*0.02));
nq1 = -1;
nq9 = ns(round(nPairs*0.98));

noiseSel = noise2(noise2>nq1 & noise2<nq9);

% varOut = var(noise);
varOutMean = mean(noiseSel);
varOutMode = median(noiseSel);