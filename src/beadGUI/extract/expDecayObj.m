function f = expDecayObj( v,x,y )
%FITEXPDECAY compute the objective function of exponential decay

a = v(1);
tau = v(2);

c0 = 1 - exp(-x/tau);
f = sum(-2*y*a.*c0 + a^2*c0.^2);

end

