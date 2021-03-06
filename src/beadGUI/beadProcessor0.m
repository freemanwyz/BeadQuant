function handles = beadProcessor0(handles)
%Bead ProcessingScript
% !! OLD
%It will then automatically fit time series for each bead with a decaying
%exponential, returning the asymptotic value as a max DFF estimate and tau
%as an estimation of speed with which bead changes fluorescence.
%
%User then selects beads based on desirable DFF and tau combinations.
%Output is series of graphs with fits for all slected beads.

dd = handles.dd;
pp = handles.pp;

M = [transpose(dd.myCurves),reshape(dd.meanCurve,[],1)];

% Making the time vector using size of the M matrix
time_prep = 0:(size(M,1)-1);
time_vect = double(time_prep*pp.timeStep);
% tt = reshape(time_vect,[],1);

% Calculate and plot DFF for all beads
dff_mat = double(calcDF(M,1,pp.useDff,1));

if size(dff_mat,2) > 10000
    pixVal = max(dff_mat,[],1);
    [~,I] = sort(pixVal,'descend');
    idxSel1 = I(1:1000);
    idxSel = [randperm(size(dff_mat,2),10000),idxSel1];
    idxSel = unique(idxSel);
    dff_mat_show = dff_mat(:,idxSel);
else
    dff_mat_show = dff_mat;
end

if 1
    plot(time_vect,dff_mat_show);
    xlabel('time (sec)');
    if pp.useDff
        ylabel('\DeltaF/F');
    else
        ylabel('\DeltaF');
    end
end

% save figure
if 1
    f0 = figure('visible','off');
    plot(time_vect,dff_mat_show);
    xlabel('time (sec)');
    if pp.useDff
        ylabel('\DeltaF/F');
    else
        ylabel('\DeltaF');
    end
    saveas(f0,[pp.outPath,filesep,'all_beads_curves.tif']);
    saveas(f0,[pp.outPath,filesep,'all_beads_curves'],'epsc');
    close(f0);
end
set(handles.myinfo,'String','Curves of all beads');
pause(0.1);

% keyboard

% Fit all of the data. Return maxDFF, tau, and R^2

nBeadsAll = size(dff_mat,2);
fprintf('%d beads to fit\n',nBeadsAll);

dff_fit = zeros(1,size(dff_mat,2));
tau_fit = zeros(1,size(dff_mat,2));
if pp.usePara
    h = msgbox('Fitting curves');
    tic
    parfor i = 1:size(dff_mat,2)
        try
            if pp.fitModel==1
                [cf,~] = DecayingExponentialFit(time_vect,dff_mat(:,i));
                % Inputting the fit evaluated at the last time point for max activation
                fitVal = feval(cf,time_vect(end));
                dff_fit(i) = fitVal;
                % Using the tau fit for speed of activation
                tau_fit(i) = cf.tau;
            else
                xx = time_vect;
                yy = reshape(dff_mat(:,i),1,[]);
                coeffs = polyfit(xx, yy, 1);
                dff_fit(i) = coeffs(1);  % slope
                tau_fit(i) = 0;  % intercept
            end
        catch
            dff_fit(i) = 0;
            tau_fit(i) = 0;
        end
    end
    tt = toc;
    fprintf('%f\n',tt);
    if isvalid(h)
        close(h)
    end
else
    h = waitbar(0,'Fitting curves...');
    tic
%     for i = 1:size(dff_mat,2)
    for i = 1:1000
        try
            if pp.fitModel==1
                [cf,~] = DecayingExponentialFit(time_vect,dff_mat(:,i));
                fitVal = feval(cf,time_vect(end));
                dff_fit(i) = fitVal;
                tau_fit(i) = cf.tau;                
            else
                xx = time_vect;
                yy = reshape(dff_mat(:,i),1,[]);
                coeffs = polyfit(xx, yy, 1);
                dff_fit(i) = coeffs(1);  % slope
                tau_fit(i) = 0;  % intercept
            end
        catch
            dff_fit(i) = 0;
            tau_fit(i) = 0;
        end
        h = waitbar(i/size(dff_mat,2));
    end
    tt = toc;
    fprintf('%f\n',tt);
    if isvalid(h)
        close(h)
    end
end

dff_fit_mean = dff_fit(end);
tau_fit_mean = tau_fit(end);
dff_fit = dff_fit(1:(end-1));
tau_fit = tau_fit(1:(end-1));

handles.dd.dff_fit = dff_fit;
handles.dd.tau_fit = tau_fit;
handles.dd.dff_fit_mean = dff_fit_mean;
handles.dd.tau_fit_mean = tau_fit_mean;
handles.dd.dff_mat = dff_mat;
handles.dd.time_vect = time_vect;

end




