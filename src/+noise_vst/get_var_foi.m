function [fitparams,cx,cy] = get_var_foi(z)
% ALGORITHM MAIN PARAMETERS
%                                %  the standard-deviation function has the form \sigma=(a*y^polyorder+b*y^(polyorder-1)+c*y^(polyorder-2)+...).^(variance_power/2), where y is the unclipped noise-free signal.
polyorder=1;                     %  order of the polynomial model to be estimated [default 1, i.e. affine/linear]  Note: a large order results in overfitting and difficult and slow convergence of the ML optimization.
variance_power=1;                %  power of the variance [default 1, i.e. affine/linear variance]
%                                %   The usual Poissonian-Gaussian model has the form \sigma=sqrt(a*y+b), which follows from setting polyorder=1 and variance_power=1.

median_est=0;                    %  0: sample standard deviation;  1: MAD   (1)
LS_median_size=0;                %  size of median filter for outlier removal in LS fitting (enhances robustness for initialization of ML) 0: disabled  [default 1 = auto]
tau_threshold_initial=1;         %  (initial) scaling factor for the tau threshold used to define the set of smoothness   [default 1]

% prior_density=0; %(SET BELOW)  %  type of prior density to use for ML    (0)
%                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
%                                %    1: zero_one uniform prior density [0,1];
%                                %    2: -infty_infty uniform prior density (R);

level_set_density_factor=1;      %   density of the slices in for the expectations   [default 1 ( = 600 slices)]   (if image is really small it should be reduced, say, to 0.5 or less)
integral_resolution_factor=1;    %   integral resolution (sampling) for the finite sums used for evaluatiing the ML-integral   [default 1]
speed_factor=1;                  %   factor controlling simultaneously density and integral resolution factors  [default 1] (use large number, e.g. 1000, for very fast algorithm)

text_verbosity=1;                %  how much info to print to screen 0: none, 1: little, 2: a lot
figure_verbosity=3;              %  show/keep figures?        [default 3]
%                                     0: no figures
%                                     1: only figure with final ML result is shown and kept
%                                     2: few figures are shown during processing but none kept
%                                     3: few figures are shown but only figure with final ML result is kept
%                                     4: show and keep all figures
lambda=0;                        %  [0,1]:  models the data distribution as a mixture of Gaussian and Cauchy PDFs (each with scale parameter sigma), with mixture parameters lambda and (1-lambda): p(x) = (1-lambda)*N(x,y,sigma^2)+lambda*Cauchy(x,y,sigma)     [default 1]
auto_lambda=0;                   %  include the mixture parameter lambda in the maximization of the likelihood   [default 1]
%                                     0: lambda is fixed equal to its input value and not optimized
%                                     1: lambda is optimized (the input value of lambda is used as initial value for the optimization)

clipping_below=0;   %  on/off   [keep off for pure-poissonian (no gaussian terms) noise, since there are no negative errors]
clipping_above=0;   %  on/off
prior_density=1;                 %  type of prior density to use for ML    (0)

% CALL ESTIMATION FUNCTION WITH GIVEN OBSERVATION AND PARAMETERS
fitparams=function_ClipPoisGaus_stdEst2D(z,polyorder,variance_power,clipping_below,...
    clipping_above,prior_density,median_est,LS_median_size,tau_threshold_initial,...
    level_set_density_factor,integral_resolution_factor,speed_factor,text_verbosity,...
    figure_verbosity,lambda,auto_lambda);
% fitparams = [ \hat{a} and \hat{b} ]  contains the estimates of the a,b parameters.

% Get plots
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
cx = get(dataObjs,'XData');  % cx{4} is clippled curve, cx{5} is linear line
cy = get(dataObjs,'YData');

end

