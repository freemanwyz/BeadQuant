%% test the function
options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on','Display','off');
% options = optimoptions('fmincon','GradObj','on','Display','off');
% options = optimoptions('fmincon','Algorithm','sqp','GradObj','off','Display','off');
% options = optimoptions('fmincon','Algorithm','sqp','GradObj','on','Display','off');

dff_mat = dd.dff_mat;
nBeads = size(dff_mat,2);
x = dd.time_vect;
tic
a0x = zeros(nBeads,1);
tau0x = zeros(nBeads,1);
dffx = zeros(nBeads,1);
taux = zeros(nBeads,1);
% for ii=682
parfor ii=1:nBeads
%     fprintf('%d\n',ii);
    [ dffFit, tauFit, a0, tau0 ] = expDecayFit( x, dff_mat(:,ii), options );
    dffx(ii) = dffFit;
    taux(ii) = tauFit;
    a0x(ii) = a0;
    tau0x(ii) = tau0;
end
toc

figure;scatter(dffx,taux);
figure;scatter(a0x,tau0x);

%% test fitting exponential decay
a = 0.9;
x = 0:30;
tau = 5;
% x = 0:5;
% tau = 2;
noiseSd = 0.01;
% y = y - y(1);

nTps = length(x);
% v0 = [0.8,4];

tic
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
for ii=1:1000
    y = a*(1-exp(-x/tau)) + randn(1,length(x))*noiseSd;
    a0 = max(y)*1.05;
    tau0 = mean(x(2:(nTps-1))./(-log(1-y(2:(nTps-1))/a0)));
    v0 = [a0,tau0];
    f = @(v)expDecayObj(v,x,y);
    [res,fval] = fminunc(f,v0,options);
end
toc

tic
options = optimoptions('fminunc','Algorithm','quasi-newton','GradObj',...
    'on','Display','off','TolFun',1e-6,'TolX',1e-6);
for ii=1:1000
    y = a*(1-exp(-x/tau)) + randn(1,length(x))*noiseSd;
    a0 = max(y)*1.05;
    tau0 = mean(x(2:(nTps-1))./(-log(1-y(2:(nTps-1))/a0)));
    v0 = [a0,tau0];
    f1 = @(v)expDecayObjWithGrd(v,x,y);
    [res,fval] = fminunc(f1,v0,options);
end
toc