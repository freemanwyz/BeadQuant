function [ rga,rgb ] = myCropRg( pixVec, Nx, Ny, varargin )

if nargin==3
    ofst = 5;
else
    ofst = varargin{1};
end

dat1Ele = zeros(Nx,Ny);
dat1Ele(pixVec) = 1;
[ai,bi] = find(dat1Ele>0);
rga = (min(ai)-ofst):(max(ai)+ofst);
rgb = (min(bi)-ofst):(max(bi)+ofst);
rga = rga(rga>0 & rga<Nx);
rgb = rgb(rgb>0 & rgb<Ny);

end

