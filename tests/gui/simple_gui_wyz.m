function simple_gui_wyz
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

% components -----
set(0,'defaultfigureunits','characters');
set(0,'defaultuicontrolunits','characters');
set(0,'defaultaxesunits','characters');
f = figure('Visible','off','Position',[1,1,100,30],...
    'MenuBar','none','ToolBar','none','NumberTitle','off');
set(f,'Color', get(0,'DefaultUicontrolBackgroundColor'));
hsurf = uicontrol('Style','pushbutton','String','Surf','Position',[80,22,10,2],'Callback',{@surfbutton_Callback});
hmesh = uicontrol('Style','pushbutton','String','Mesh','Position',[80,18,10,2],'Callback',@meshbutton_Callback);
hcontour = uicontrol('Style','pushbutton','String','Countour','Position',[80,14,10,2],'Callback',@contourbutton_Callback);
htext = uicontrol('Style','text','String','Select Data','Position',[75,10,20,1]);
hpopup = uicontrol('Style','popupmenu','String',{'Peaks','Membrane','Sinc'}, ...
    'Position',[75,6,20,2],'Callback',@popup_menu_Callback);
align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
% ha = axes('Position',[10,5,60,20]);
ha = axes('Parent',f,'Units','normalized','Position',[.1 .15 .6 .7]);

% init -----
set(f,'Units','normalized');
set(ha,'Units','normalized');
set(hsurf,'Units','normalized');
set(hmesh,'Units','normalized');
set(hcontour,'Units','normalized');
set(htext,'Units','normalized');
set(hpopup,'Units','normalized');

peaks_data = peaks(35);
membrane_data = membrane;
[x,y] = meshgrid(-8:.5:8);
r = sqrt(x.^2+y.^2) + eps;
sinc_data = sin(r)./r;
current_data = peaks_data;
surf(current_data);
set(f,'Name','Simple GUI');
movegui(f,'center')
set(f,'Visible','on')

% call backs -----
    function popup_menu_Callback(source,~)
        % Determine the selected data set.
        str = get(source,'String');
        val = get(source,'Value');
        switch str{val};
            case 'Peaks' % User selects Peaks.
                current_data = peaks_data;
            case 'Membrane' % User selects Membrane.
                current_data = membrane_data;
            case 'Sinc' % User selects Sinc.
                current_data = sinc_data;
        end
    end

    function surfbutton_Callback(source,~)
        surf(current_data);
    end

    function meshbutton_Callback(source,~)
        mesh(current_data);
    end

    function contourbutton_Callback(source,~)
        contour(current_data);
    end


end