function [ dffFit, tauFit, a0, tau0 ] = expDecayFit( x, y, options )
%EXPDECAYFIT fit time course with exponential decay
% Output is the fit value for last time point and tau
% Unconstraint with good init and little post filtering is good enough.
% I prefer quasi-newton for fminunc
% Constraint method is also good, but slow. We use the interior point, the
% sgp method can not update tau, it may be stuck at local minimum

dffFit = 0;
tauFit = 0;

x = reshape(x,[],1);
y = reshape(y,[],1);
% nTps = length(x);

%% initial value
[~,I] = max(abs(y));
if abs(y(I)) < 1/65535
    return
end
a0 = y(I)*1.01;
% a0 = max(y)*1.1;

% ySel = y > y(I)*0.3;
if a0>0
    ySel = y>0;
else
    ySel = y<0;
end
% ySel = y~=0;
% ySel(I) = 0;
tau0 = mean(x(ySel)./(-log(1-y(ySel)/a0)));
% tau0 = mean(x(2:(nTps-1))./(-log(1-y(2:(nTps-1))/a0)));
% tau0 = x(I)/(-log(1-y(I)/a0));  % use single point to estimate tau

if tau0 > 2*max(x)
    tau0 = 2*max(x);
end
if tau0 < 0.01*x(2)
    tau0 = 0.01*x(2);
%     uu = rand()*0.1 + 0.05;
%     tau0 = uu*x(2);
end

% if tau0>2*x(end)
%     tau0 = 2*x(end);
%     a0 = mean(y(2:(nTps-1))./(1-exp(-x(2:(nTps-1))/tau0)));
% end

v0 = [a0,tau0];  % !!

%% curve fitting
f1 = @(v)expDecayObjWithGrd(v,x,y);
% f1 = @(v)expDecayObj(v,x,y);  % !!

res = fminunc(f1,v0,options);

% lb = [0,0];
% ub = [min(2*a0,1),max(x)*2];
% res = fmincon(f1,v0,[],[],[],[],lb,ub,[],options);

aFit = res(1);
tauFit = res(2);

% if aFit<0 || tauFit<0
%     tauFit = 0;
%     return
% end

dffFit = aFit*(1-exp(-x(end)/tauFit));

if tauFit > 2*max(x)
    tauFit = 2*max(x);
end

% if aFit>a0*1.5 || tauFit>2*x(end) || aFit<0 || tauFit<0
%     aFit = a0;
%     tauFit = tau0;
% end

if 0  % !!
    yFit = aFit*(1-exp(-x/tauFit));
    figure;plot(x,y);hold on;plot(x,yFit);
end

end

