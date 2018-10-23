function [objf, devf] = expDecayObjWithGrd( v,x,y )
%FITEXPDECAY compute the objective function of exponential decay

a = v(1);
tau = v(2);

b0 = exp(-x/tau);
c0 = 1 - b0;

objf = sum(-2*y*a.*c0 + a^2*c0.^2);

df_a = sum(-2*y.*c0 + 2*a*c0.^2);
df_tau = sum(2*a*x.*b0.*(y-a+a*b0))/tau^2;

devf = [df_a,df_tau];

end