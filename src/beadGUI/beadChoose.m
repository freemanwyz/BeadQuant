function goodBeadIdx = beadChoose(dff_fit,tau_fit,...
    dff_mat,time_vect,controlBeadIdx,datMode,fitModel,outPath,handles)
%% Visualize Data
%Plot ROIs with bad fits
badROI = dff_fit==0 & tau_fit==0;

if sum(badROI)>0
    f = figure('visible','off');
    plot(time_vect,dff_mat(:,badROI));
    if fitModel==1
        xlabel('time (sec)');
        if datMode
            ylabel('\DeltaF/F');
        else
            ylabel('\DeltaF');
        end
    else
        xlabel('k');
    end
    title('Plots of ROIs with failed fits');
    saveas(f,[outPath,filesep,'bad_beads_curves.tif']);
    saveas(f,[outPath,filesep,'bad_beads_curves'],'epsc');
    delete(f);
end

%% Select region of points for analysis
%Select the region of plot where interesting beads are
f = figure('visible','off');
if fitModel==1
    scatter(dff_fit,tau_fit);
    ylabel('\tau')
    if datMode
        xlabel('\DeltaF/F');
    else
        xlabel('\DeltaF');
    end
else
    scatter(dff_fit,tau_fit*0);
    xlabel('k');
end
saveas(gcf,[outPath,filesep,'beads_dff_tau_scatter.tif']);
saveas(gcf,[outPath,filesep,'beads_dff_tau_scatter'],'epsc');
close(f)
if fitModel==1
    scatter(dff_fit,tau_fit);
    ylabel('\tau')
    if datMode
        xlabel('\DeltaF/F');
    else
        xlabel('\DeltaF');
    end
else
    scatter(dff_fit,tau_fit*0);
    xlabel('k');
end
h = uicontrol('Position',[20 20 60 25],'String','Continue',...
    'Callback','uiresume(gcbf)');
% title('Zoom and center plot on points of interest, then hit "Continue"')
set(handles.myinfo,'String','Zoom and center plot on points of interest, then hit "Continue"');
uiwait;
delete(h);
axLim = [get(gca,'xlim'),get(gca,'ylim')];
% close(f);

%Select beads of interest by enclosing within a rectangle
set(handles.myinfo,'String','Draw a box and type 0 to continue');
FIGDisp = 1;%loop until user happy with rectangle
while FIGDisp==1
%     f = figure;
    scatter(dff_fit,tau_fit);
    axis(axLim)
    ylabel('\tau')
    if datMode
        xlabel('\DeltaF/F');
    else
        xlabel('\DeltaF');
    end
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    plot(x,y,'r','linewidth',2)          % draw box around selected region
    prompt = 'Input 1 to select a different box, 0 to continue';
    FIGDisp = str2double(cell2mat(inputdlg(prompt)));
    hold off
%     close(f);
end

%% Get and make plots for selected beads
set(handles.myinfo,'String','');
%Get selected beads
goodDFF = dff_fit>min(x) & dff_fit<max(x);
goodTau = tau_fit>min(y) & tau_fit<max(y);
goodBeads = goodDFF & goodTau;
goodBeadIdx = find(goodBeads);

h = msgbox('Saving curves...');
%make plots for selected beads
for ii = 1:sum(goodBeads)
    f0 = figure('visible','off');
    if fitModel==1
        DecayingExponentialFitwGraph1(time_vect,dff_mat(:,goodBeadIdx(ii)));
    else
        xx = time_vect;
        yy = reshape(dff_mat(:,goodBeadIdx(ii)),1,[]);
        scatter(xx, yy, 10, [0 0.5 0], 'filled');
        coeffs = polyfit(xx, yy, 1);
        fittedX = linspace(min(xx), max(xx), 200);
        fittedY = polyval(coeffs, fittedX);
        hold on;
        plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
        legend('Fl vs time','Linear fit','Location','southeast');
    end
    xlabel('time (sec)')
    if datMode
        ylabel('\DeltaF/F');
    else
        ylabel('\DeltaF');
    end
    titStr = sprintf('Plot for Bead %d',goodBeadIdx(ii));
    title(titStr)
    %     set(f0, 'Visible', 'off')
    saveas(f0,[outPath,filesep,'bead_',num2str(goodBeadIdx(ii)),'_curves.tif']);
    saveas(f0,[outPath,filesep,'bead_',num2str(goodBeadIdx(ii)),'_curves'],'epsc');
    close(f0);
end

%plot control beads curves
for ii = 1:length(controlBeadIdx)
%     f0 = figure('name','Control beads','visible','off');
    f0 = figure('visible','off');
    if fitModel==1
        DecayingExponentialFitwGraph1(time_vect,dff_mat(:,controlBeadIdx(ii)));
    else
        xx = time_vect;
        yy = reshape(dff_mat(:,controlBeadIdx(ii)),1,[]);
        scatter(xx, yy, 10, [0 0.5 0], 'filled');
        coeffs = polyfit(xx, yy, 1);
        fittedX = linspace(min(xx), max(xx), 200);
        fittedY = polyval(coeffs, fittedX);
        hold on;
        plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
        legend('Fl vs time','Linear fit','Location','southeast');
    end
    xlabel('time (sec)')
    if datMode
        ylabel('\DeltaF/F');
    else
        ylabel('\DeltaF');
    end
    titStr = sprintf('Plot for Control Bead %d',controlBeadIdx(ii));
    title(titStr)
%     set(f0, 'Visible', 'off')
    saveas(gcf,[outPath,filesep,'bead_control_',num2str(controlBeadIdx(ii)),'_curves.tif']);
    saveas(gcf,[outPath,filesep,'bead_control_',num2str(controlBeadIdx(ii)),'_curves'],'epsc');
    close(f0);
end
close(h);

cla reset
pause(0.25)

end




