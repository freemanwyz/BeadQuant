function res = pgL(pCur,yest,sigmaest,ci,di)
% PGL compute numerical integration for Foi's non-clipped likelihood

a = pCur(1);
b = pCur(2);
M = 10000;           % integration step
N = length(yest);   % number of level sets
yrg = 0:1/M:1;
res = 0;
for ii=1:N
    sigma2reg = max(a*yrg+b,1e-6);
    k = 1./(2 * pi * sqrt(ci(ii)*di(ii)) * sigma2reg);     
    tc = (yest(ii)-yrg).^2/ci(ii);
    td = (sigmaest(ii)-sqrt(sigma2reg)).^2/di(ii);
    t = -1./(2*sigma2reg).*(tc+td);
    res = res + trapz(yrg,k.*exp(t));    
end
res = -res;

end