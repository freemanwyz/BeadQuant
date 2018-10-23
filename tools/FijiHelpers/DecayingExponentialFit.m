function [cf_,gof] = DecayingExponentialFit(time_vect,beadn)
time_vect = time_vect(:);
beadn = beadn(:);

a0 = max(beadn);
t0 = max(time_vect);
amin = 0;
amax = a0*1.5;
tmin = 0;
tmax = t0*2;

% fo_ = fitoptions('method','NonlinearLeastSquares');
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[amin tmin],'Upper',[amax tmax]);
ok_ = isfinite(time_vect) & isfinite(beadn);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
st_ = [beadn(end) time_vect(2)];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*(-exp(-x/tau)+1)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'tau'});

% Fit this model using new data
[cf_,gof] = fit(time_vect(ok_),beadn(ok_),ft_,fo_);



