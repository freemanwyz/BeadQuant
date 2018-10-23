function [W, psi, Ez, mu, llh] = fa(X, p)
% EM algorithm for factor analysis model
% reference:
% Pattern Recognition and Machine Learning by Christopher M. Bishop 

Y = X - 0.05;
[W,m5] = nnmf(Y,1);
dif0 = Y - W * m5;
invpsi = 1./var(dif0,0,2);

[d,n] = size(X);
mu = mean(X,2);
X = bsxfun(@minus,X,mu);

tol = 1e-6;
converged = false;
llh = -inf;

% initialize parameters
% W = rand(d,p); 
% invpsi = 1./rand(d,1);
% invpsi = 5./rand(d,1);
% W = mu;
% invpsi = ones(d,1)*5;

% precompute quantities
I = eye(p);
normX = sum(X.^2,2);

U = bsxfun(@times,W,sqrt(invpsi));
M = U'*U+I;                     % M = W'*inv(Psi)*W+I
R = chol(M);
invM = R\(R'\I);
WinvPsiX = bsxfun(@times,W,invpsi)'*X;       % WinvPsiX = W'*inv(Psi)*X
ii = 0;
while ~converged && ii<100
    ii = ii + 1;
    % E step
    Ez = invM*WinvPsiX;
    Ezz = n*invM+Ez*Ez';
    % end
    
    R = chol(Ezz);
    XEz = X*Ez';
    
    % M step
    W = (XEz/R)/R';
    dif0 = normX-sum(W.*XEz,2);
%     dif0(dif0<0.05) = 0.05;  % lower bound on noise
    invpsi = n./dif0;
    % end

    % compute quantities needed
    U = bsxfun(@times,W,sqrt(invpsi));
    M = U'*U+I;                     % M = W'*inv(Psi)*W+I
    R = chol(M);
    invM = R\(R'\I);
    WinvPsiX = bsxfun(@times,W,invpsi)'*X;       % WinvPsiX = W'*inv(Psi)*X
    % end
    
    % likelihood
    last = llh;
    logdetC = 2*sum(log(diag(R)))-sum(log(invpsi));              % log(det(C))
    trinvCS = (normX'*invpsi-sum(sum((R'\WinvPsiX).^2)))/n;  % trace(inv(C)*S)
    llh = -n*(d*log(2*pi)+logdetC+trinvCS)/2;
    % end
    converged = abs(llh-last) < tol*abs(llh);   % check likelihood for convergence
end
psi = 1./invpsi;
